#!/usr/bin/Rscript

### CV
# stop when there is an warning
ptm <- proc.time()

options(warn=2)

# load library
## install SparkR
#library(devtools)
#install_github("amplab-extras/SparkR-pkg", subdir="pkg")

library(SparkR)

# read in command arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
   print("Usage : ML_parallel.R <master> <cv_file> <nslices=1>")
   print("Example : Rscript ML_parallel.R local[10] cv_small.txt")
   q("no")
}


# get list of cross validation parameters

cvpar = readLines(args[[2]])

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


#Initialize Spark context 
sc <- sparkR.init(args[[1]], "geosurvey_predict", sparkEnvir=list(spark.executor.memory="1g"))

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

nslices = 1
if(length(args)==3) {nslices = as.numeric(args[[3]])}

rdd <- parallelize(sc, coll=cvpar, numSlices = nslices)
myerr <- flatMap(rdd,geosurvey_data)
output <- collect(myerr)

# write.csv(results, "results.csv", row.names=FALSE)



