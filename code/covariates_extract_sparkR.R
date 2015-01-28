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

#Initialize Spark context 
sc <- sparkR.init(args[[1]], "geosurvey_predict", sparkEnvir=list(spark.executor.memory="1g"))

geosurvey_data <- function(indata) {
   library(raster)
   print("In cv_bart")
   print(indata)
   if(indata!=" "){
       filename <- indata
       locs <- read.csv(filename)
       locs <- cbind(locs$x, locs$y)
       extract_data <- extract(grid, locs)
       extract_data <- cbind(locs, extract_data)
       colnames(extract_data) <- c("Longitude", "Latitude",  names(grid)) 
       write.csv(extract_data, paste(filename, "output", sep=""), row.names=FALSE)
       myres = paste("finished", filename, sep=":")
   }
}

nslices = 1
if(length(args)==3) {nslices = as.numeric(args[[3]])}

rdd <- parallelize(sc, coll=cvpar, numSlices = nslices)
myerr <- flatMap(rdd,geosurvey_data)
output <- collect(myerr)

# write.csv(results, "results.csv", row.names=FALSE)



