# Generate GeoSurvey sampling points in Africa on a 1 km grid
# J. Chen, W. Wu & M. Walsh, November, 2014

require(downloader)
require(rgdal)

# Download sampling mask grid ---------------------------------------------
download("https://www.dropbox.com/s/brikq1ceek16w3s/Af_mask.tif.zip?dl=0", "AF_mask.tif.zip", mode="wb")
unzip("Af_mask.tif.zip", overwrite=T)
mask <- readGDAL("AF_mask.tif", silent=TRUE)

# Generate (n) sample and (v) validation locations ------------------------
n <- 1000000
v <- n/4
set.seed(1385321)
coord_mask <- coordinates(mask)[!is.na(mask@data$band1), ]
mask_n <- sample(1:dim(coord_mask)[1], n)
coord_mask_n <- coord_mask[mask_n, ]
coord_mask_n <- project(coord_mask_n, "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs", inv = TRUE)
coord_mask_v <- sample(1:n, v)
coord_mask_v <- coord_mask_n[coord_mask_v, ]

# Write csv's -------------------------------------------------------------
k <- 4
nk <- n/k
for(i in 1:k){
	write.csv(data.frame(Lon=coord_mask_n[(((i-1)*nk+1):(i*nk)),1], Lat=coord_mask_n[(((i-1)*nk+1):(i*nk)),2]), paste("AF_GS_sample_4p_", i, ".csv", sep=""), row.names=FALSE, quote=FALSE)
}

k <- 2
nk <- n/k
for(i in 1:k){
	write.csv(data.frame(Lon=coord_mask_n[(((i-1)*nk+1):(i*nk)),1], Lat=coord_mask_n[(((i-1)*nk+1):(i*nk)),2]), paste("AF_GS_sample_2p_", i,".csv", sep=""), row.names=FALSE, quote=FALSE)
}

write.csv(data.frame(Lon=coord_mask_v[,1], Lat=coord_mask_v[,2]), "AF_GS_valpts.csv", row.names=FALSE, quote=FALSE)



