library(rgdal)
library(downloader)

download("https://www.dropbox.com/s/6wk1dmz8xxo6vdw/fAPAR.zip?dl=0", "fAPAR.zip", mode="wb")
unzip("fAPAR.zip", overwrite=T)

fpar_ave <- readGDAL("fPAR_ave.tif", silent=TRUE)
coordinates_fpar <- coordinates(fpar_ave)[!is.na(fpar_ave@data$band1), ]
coordinates_fpar_selected <- sample(1:dim(coordinates_fpar)[1], 1000000)
coordinates_fpar_selected <- coordinates_fpar[coordinates_fpar_selected, ]

coordinates_fpar_doublecheck <- sample(1:1000000, 500000)
coordinates_fpar_doublecheck <- coordinates_fpar_selected[coordinates_fpar_doublecheck, ]

