require(rgdal)

geosurvey_a1 <- read.csv("../data/2015-01-20_18-00-02_248261/6_1m-point-survey-a1.csv", stringsAsFactors=FALSE)
geosurvey_a2 <- read.csv("../data/2015-01-20_18-00-02_248261/7_1m-point-survey-a2.csv", stringsAsFactors=FALSE)

geosurvey_total <- rbind(geosurvey_a1, geosurvey_a2)

geosurvey_proj <- as.data.frame(project(cbind(geosurvey_total$Longitude, geosurvey_total$Latitude), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geosurvey_proj) <- c("x","y")
geosurvey_total <- cbind(geosurvey_total, geosurvey_proj)

write.csv(geosurvey_total, "geosurveydata_total.csv")


n <- dim(geosurvey_total)[1]
rownum <- round(seq(1, n, length.out=100)[-1])

file.create("geofiles.txt")

write.csv(geosurvey_total[(1:rownum[1]), ], "geosurveydata_1.csv", row.names=FALSE) 
write.table("geosurveydata_1.csv", "geofiles.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
for(i in 2:(length(rownum))){
    write.csv(geosurvey_total[((rownum[i-1]+1):rownum[i]),], paste( "geosurveydata_",i, ".csv", sep=""), row.names=FALSE) 
    write.table(paste( "geosurveydata_",i, ".csv", sep=""), "geofiles.txt", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
}


