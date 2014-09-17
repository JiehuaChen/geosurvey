require(downloader)
require(proj4)
require(raster)

# Data downloads ----------------------------------------------------------

# Geosurvey data

download("https://www.dropbox.com/s/8dgzphvwn0ls4sq/TZ_geos_0914.csv?dl=0", "../data/TZ_geos_0914.csv", mode="wb")


# Tanzania Gtifs (~35.5 Mb)
download("https://www.dropbox.com/s/0ucx0oqx8c4lpej/TZ_grids.zip?dl=0", "../data/TZ_grids.zip", mode="wb")
unzip("../data/TZ_grids.zip", exdir="../data", overwrite=T)

