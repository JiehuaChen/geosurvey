#!/usr/bin/Rscript



line <- commandArgs(trailingOnly=TRUE)
line <- line[1]

# load library
library(raster)
library(rgdal)
library(h2o)

# fetch tif files from s3
covnames <-  c("AFPOPs", "BSANs",  "BSASs",  "BSAVs", "CTIs",   "ELEVs", "EVIs",  "FPAR", "GROADs",
               "GWCRP", "IUCN2", "IUCN4", "IUCNNR",  "LSTDs",  "LSTNs", "REF1s",  "REF2s", "REF3s", "REF7s",  "RELIs", "TARDs", "TMAPs",  "TMFIs")


tiffiles <- paste(paste(covnames, line, sep="_"),".tif", sep="")
for(tiffile in tiffiles){
	tiffile_download <- 1
	while (tiffile_download != 0){
		tiffile_download <- system(paste("s3cmd get -f s3://afsis.geosurvey/data/dataForEachNode/", tiffile, sep=""))
		Sys.sleep(5)
	}

}


grid <- stack(tiffiles)

fpar_file <-  paste(paste("FPAR_mask", line,sep="_"), ".tif", sep="")

system(paste("s3cmd get --skip-existing s3://afsis.geosurvey/data/dataForEachNode/fpar_mask/", fpar_file, sep=""), intern=TRUE)

tif.info <- GDALinfo(fpar_file, silent=TRUE)
totalrows <- tif.info[1]
totalcols <- tif.info[2]
bandsnum <- tif.info[3]
origin.x <- tif.info[4]
res.x <- tif.info[6]
res.y <- tif.info[7]
if(attr(tif.info, "ysign")==-1){
	#origin of data array should be the left upper corner
	origin.y <- tif.info[5] + totalrows*res.y
}else{
	origin.y <- tif.info[5]	
}

locs <- coordinates(grid)
locs <- locs[locs[,2]==median(locs[,2]), ]


fpar_values <- readGDAL(fpar_file, offset=c(1,0), region.dim=c(1,totalcols))

locs <- locs[!is.na(fpar_values@data$band1), ]
locs <- locs[(fpar_values@data$band1)[!is.na(fpar_values@data$band1)]==1, ]
if(!is.matrix(locs)){
    locs <- t(as.matrix(locs))
}
pixeln.x <- ceiling((locs[,1]-origin.x)/res.x)-1
pixeln.y <- ceiling((origin.y-locs[,2])/res.y)-1
neighbor.x <- c(-1, 0, 1, -1,0,1, -1, 0, 1)
neighbor.y <- c(-1, -1, -1, 0, 0, 0, 1, 1, 1)
neighbor_locations_gid  <-  matrix(NA, dim(locs)[1], 18)

for(i in 1:9){
    neighbor_locations_gid[, (2*i-1)] <- (pixeln.x+neighbor.x[i])*res.x + origin.x + res.x/2
    neighbor_locations_gid[, (2*i)] <- origin.y - ((pixeln.y+neighbor.y[i])*res.y + res.y/2)
}

names_neighbors <- c("sw" ,"s", "se", "w","loc", "e", "nw", "n", "ne")

cov_neighbors <- NULL
for(i in 1:9){
    extract_data <- extract(grid, neighbor_locations_gid[, (2*i-1):(2*i)])
    cov_neighbors <- cbind(cov_neighbors, extract_data)
}
cov_neighbors <- as.data.frame(cov_neighbors)
names(cov_neighbors) <- paste(covnames, rep(names_neighbors, each=length(covnames)), sep="_")

cov_neighbors_withloc <- cbind(locs, cov_neighbors)
cov_neighbors_withloc <- na.omit(cov_neighbors_withloc)
if(dim(cov_neighbors_withloc)[1]>0){
    cov_neighbors <- cov_neighbors_withloc[, -c(1,2)]
    locs <- cov_neighbors_withloc[, c(1:2)]

    kernel_predicts <- extract(stack("/home/ubuntu/kernel_maps/RSPKDs.tif"), locs) 
    
    # load h2o models
    localH2O = h2o.init()

    cov_neighbors.h2o <- as.h2o(localH2O, cov_neighbors, key="predict_cov")

    glm_ens <- h2o.loadModel(localH2O, "/home/ubuntu/h2omodels/glm_ens_hs_withkernel")
    ens_coef <- glm_ens@models[[glm_ens@best_model]]@model$coefficients

    rf_model <- h2o.loadModel(localH2O, "/home/ubuntu/h2omodels/rf_model_hsp")
    gbm_model <- h2o.loadModel(localH2O, "/home/ubuntu/h2omodels/gbm_model_hsp")
    nnet_model <- h2o.loadModel(localH2O, "/home/ubuntu/h2omodels/nnet_model_hsp")


    gbm_predicts <- as.data.frame.H2OParsedData(h2o.predict(gbm_model, cov_neighbors.h2o))$Yes

    rf_predicts <- as.data.frame.H2OParsedData(h2o.predict(rf_model, cov_neighbors.h2o))$Yes

    nnet_predicts <- as.data.frame.H2OParsedData(h2o.predict(nnet_model, cov_neighbors.h2o))$Yes

    ens_cov <- cbind(kernel=kernel_predicts, gbm = gbm_predicts, rf=rf_predicts, nnet=nnet_predicts, 1)
    ens_predict <- ens_cov %*% as.matrix(ens_coef)
    ens_predict <- exp(ens_predict)/(1+exp(ens_predict))
    
    gbm_predicts <- cbind(locs, gbm_predicts)
    colnames(gbm_predicts) <- c("x", "y", "prob")
    predict_file <-  paste("predicted_values_gbm", line, "csv", sep=".")
    write.csv(gbm_predicts, predict_file, row.names=FALSE)        

	tiffile_upload <- 1
	while (tiffile_upload != 0){
        tiffile_upload <- system(paste("s3cmd put -f", predict_file, "s3://afsis.geosurvey/results_h2o_hs_withkernel/"))
    }

    system(paste("rm", predict_file))

    
    rf_predicts <- cbind(locs, rf_predicts)
    colnames(rf_predicts) <- c("x", "y", "prob")
    predict_file <-  paste("predicted_values_rf", line, "csv", sep=".")
    write.csv(rf_predicts, predict_file, row.names=FALSE)        
	tiffile_upload <- 1

	while (tiffile_upload != 0){
        tiffile_upload <- system(paste("s3cmd put -f", predict_file, "s3://afsis.geosurvey/results_h2o_hs_withkernel/"))
    }
    
    system(paste("rm", predict_file))
    
    nnet_predicts <- cbind(locs, nnet_predicts)
    colnames(nnet_predicts) <- c("x", "y", "prob")
    predict_file <-  paste("predicted_values_nnet", line, "csv", sep=".")
    write.csv(nnet_predicts, predict_file, row.names=FALSE)        
	tiffile_upload <- 1

	while (tiffile_upload != 0){
        tiffile_upload <- system(paste("s3cmd put -f", predict_file, "s3://afsis.geosurvey/results_h2o_hs_withkernel/"))
    }
    
    system(paste("rm", predict_file))
        
    ens_predict <- cbind(locs, ens_predict)
    colnames(ens_predict) <- c("x", "y", "prob")
    predict_file <-  paste("predicted_values", line, "csv", sep=".")
    write.csv(ens_predict, predict_file, row.names=FALSE)        
	tiffile_upload <- 1
	while (tiffile_upload != 0){
        tiffile_upload <- system(paste("s3cmd put -f", predict_file, "s3://afsis.geosurvey/results_h2o_hs_withkernel/"))
    }
    
    system(paste("rm", predict_file))
    
    cat(line, predict_file, "uploaded.")
}else{
    cat(line, ": no prediction file generated")
}
system("rm *.tif")



