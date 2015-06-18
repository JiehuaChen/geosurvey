library(raster)
library(rgdal)
library(h2o)


geosv <- read.table("../../data/1MQ_validation_data.csv", header=T, sep=",")

glist <- list.files(path="../../../remotesensing_data/Af_grids_std", pattern="tif", full.names=T)
grid <- stack(glist)


tif.info <- GDALinfo(glist[1], silent=TRUE)
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


geosv.proj <- as.data.frame(project(cbind(geosv$Lon, geosv$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geosv.proj) <- c("x","y")
geosv <- cbind(geosv, geosv.proj)
coordinates(geosv) <- ~x+y
projection(geosv) <- projection(grid)


geosv <- as.data.frame(geosv)

locs <- cbind(geosv$x, geosv$y)

names_geosurvey <- names(geosv)

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

for(i in 1:9){
    extract_data <- extract(grid, neighbor_locations_gid[, (2*i-1):(2*i)])
    geosv <- cbind(geosv, extract_data)
    print(i)
}
colnames(geosv) <- c(names_geosurvey, paste(names(grid), rep(names_neighbors, each=length(names(grid))), sep="_"))

geosv_cov <- geosv[, (11:199)]

localH2O = h2o.init(max_mem_size='6g')
geosv.h2o <- as.h2o(localH2O, geosv, key="COV")


    rf_model <- h2o.loadModel(localH2O, "../modelfiles/rf_model/rf_model")
    gbm_model <- h2o.loadModel(localH2O, "../modelfiles/gbm_model/gbm_model")
    nnet_model <- h2o.loadModel(localH2O, "../modelfiles/nnet_model/nnet_model")


    gbm_predicts <- as.data.frame.H2OParsedData(h2o.predict(gbm_model, geosv.h2o))$Yes
    rf_predicts <- as.data.frame.H2OParsedData(h2o.predict(rf_model, geosv.h2o))$Yes
    nnet_predicts <- as.data.frame.H2OParsedData(h2o.predict(nnet_model, geosv.h2o))$Yes

    kernel_predicts <- extract(stack("../data/GeoSurveyKD/CRPKDs.tif"), locs) 
 
geosv_ens <- cbind(geosv$CRPc, kernel_predicts, gbm_predicts, rf_predicts, nnet_predicts)
geosv_ens <- na.omit(geosv_ens)
geosv_ens[,1] <- ifelse(geosv_ens[,1]==1, 0, 1)
geosv_ens <- as.data.frame(geosv_ens)
names(geosv_ens) <- c("CRP","kernel","gbm", "rf", "nnet")

geosv_ens_crp.h2o <- as.h2o(localH2O, geosv_ens, key="geosv_env_crp")

#elastic net ensemble
h2o_glm_ens <- h2o.glm(x=2:5, y=1, family="binomial", lambda_search=TRUE, data=geosv_ens_crp.h2o)

h2o.saveModel(h2o_glm_ens, "glm_ens")

