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
glist <- list.files(path="../../../remotesensing_data/Af_Gtif_std_1k", pattern="tif", full.names=T)
grid <- stack(glist)

#Initialize Spark context 
sc <- sparkR.init(args[[1]], "geosurvey_predict")

cv_bart <- function(indata) {
    library(raster)
   print("In cv_bart")
   print(indata)
   params <- strsplit(indata, split=",")[[1]]
   extract_data <- extract(grid, cbind(as.numeric(params[1]), as.numeric(params[2])))
   myres = paste(params[1],params[2], extract_data, collapse=",")
}

nslices = 1
if(length(args)==3) {nslices = as.numeric(args[[3]])}

rdd <- parallelize(sc, coll=cvpar, numSlices = nslices)
myerr <- flatMap(rdd,cv_bart)
output <- collect(myerr)
results <- do.call("rbind", strsplit(do.call("rbind", output), split=" "))

write.csv(results, "results.csv", row.names=FALSE)

