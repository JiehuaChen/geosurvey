library(rgdal)
library(downloader)

download("https://www.dropbox.com/s/6wk1dmz8xxo6vdw/fAPAR.zip?dl=0", "fAPAR.zip", mode="wb")
unzip("fAPAR.zip", overwrite=T)

fpar_ave <- readGDAL("fPAR_ave.tif", silent=TRUE)
coordinates_fpar <- coordinates(fpar_ave)[!is.na(fpar_ave@data$band1), ]
coordinates_fpar_selected <- sample(1:dim(coordinates_fpar)[1], 1000000)
coordinates_fpar_selected <- coordinates_fpar[coordinates_fpar_selected, ]
coordinates_fpar_selected <- project(coordinates_fpar_selected, "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs", inv = TRUE)


coordinates_fpar_doublecheck <- sample(1:1000000, 500000)
coordinates_fpar_doublecheck <- coordinates_fpar_selected[coordinates_fpar_doublecheck, ]


write.csv(data.frame(longitude=coordinates_fpar_selected[, 1], latitude=coordinates_fpar_selected[,2]), "1m_fparmasked_sampled.csv", row.names=FALSE, quote=FALSE)
write.csv(data.frame(longitude=coordinates_fpar_doublecheck[, 1], latitude=coordinates_fpar_doublecheck[,2]), "1m_fparmasked_doublecheck.csv", row.names=FALSE, quote=FALSE)




