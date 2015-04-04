#!/usr/bin/Rscript

### extracting covariates in parallel
# stop when there is an warning
ptm <- proc.time()

options(warn=2)

# load library
library(rgdal)
# read in command arguments
args <- commandArgs(trailingOnly=TRUE)


# get list of cross validation parameters

cvpar = readLines(args[[1]])

library(raster)
glist <- list.files(path="../../remotesensing_data/Af_grids_std", pattern="tif", full.names=T)
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


geosurvey_data <- function(indata) {
   library(raster)
   print("In cv_bart")
   print(indata)
   if(indata!=" "){
       filename <- indata
       geosurvey_data <- read.csv(filename)
       names_geosurvey <- names(geosurvey_data)
       locs <- cbind(geosurvey_data$x, geosurvey_data$y)
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
           geosurvey_data <- cbind(geosurvey_data, extract_data)
           print(i)
       }
       colnames(geosurvey_data) <- c(names_geosurvey, paste(names(grid), rep(names_neighbors, each=length(names(grid))), sep="_"))
       write.csv(geosurvey_data, paste(filename, "output", sep=""), row.names=FALSE)
       myres = paste("finished", filename, sep=":")
   }
}

library(doParallel)
registerDoParallel(cores=2)

extract_results <- foreach(filenames=cvpar, .packages = c("raster"), .verbose=TRUE) %dopar% geosurvey_data(indata = filenames)

system(paste("awk 'FNR > 1'",  "geosurveydata_[0-9]*.csv", "> geosurveydata_withcov.csv"))
system(paste("awk 'FNR == 1'", "geosurveydata_1.csvoutput", "> header.csv"))
system("cat header.csv geosurveydata_withcov.csv > geosurveydata_withcov_withheader.csv")

