#########
setwd("C:/Users/Kelly/Desktop/Frontiers_GitHub")
rm(list=ls())


output <- read.table("Rescue_GLOBAL_RCP85_SN_50_without_DISP_8.5_zooregions.txt",h=T)
PA.orig <-read.table("PAM_filtered.txt",h=T)
coords <-PA.orig[,1:2]
PA <-PA.orig[,3:7195]
remove(PA.orig) #release RAM memory
riq <-rowSums(PA)


vetorID_rescueP <- output$ID [which (output$rescuep_save_095==1)]
rescue.P_PAM <- PA[, vetorID_rescueP]

PAM_rescueP <- cbind(coords,rescue.P_PAM)
ncol(PAM_rescueP)
str(PAM_rescueP)

write.table(PAM_rescueP, "PAM_rescueP.txt")

library(raster)
X11()
raster_riq <- rasterFromXYZ(xyz =  cbind(coords, riq), res = 1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
plot(raster_riq)
writeRaster(raster_riq,"raster_riq.asc")

x11()
raster_rescueP <- rasterFromXYZ(xyz =  cbind(coords, rowSums(rescue.P_PAM)), res = 1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
plot(raster_rescueP)
writeRaster(raster_rescueP,"raster_riq_rescueP.95.asc")

library(maptools)
data("wrld_simpl")
plot(wrld_simpl, add = T, border = "grey80")

riqueza <- rowSums(PA)
riqueza_rescuep <- rowSums(rescue.P_PAM)
riqueza_rescue_P <- ifelse(test = riqueza > 0, yes = riqueza_rescuep, no = NA)

riqueza <- ifelse(test = riqueza > 0, yes = riqueza, no = NA)


prop_perda_Rich_Rescue_P <- 1-(riqueza_rescue_P/riqueza)
raster_prop_perda_Rich_Rescue_P <- rasterFromXYZ(xyz =  cbind(coords, prop_perda_Rich_Rescue_P), res = 1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
x11()
plot(raster_prop_perda_Rich_Rescue_P)
writeRaster(raster_prop_perda_Rich_Rescue_P,"raster_prop_perda_Rich_Rescue_P.asc")


#########################
library(raster)

paleartic<- output$P_RescueP [which (output$zooregions=="PALEARTIC")]
neartic <- output$Rescue [which (output$zooregions=="NEARTIC")]
neotropical <- output$Rescue [which (output$zooregions=="NEOTROPICAL")]
australia <- output$Rescue [which (output$zooregions=="AUSTRALIA")]
afrotropical <- output$Rescue [which (output$zooregions=="AFROTROPICAL")]
malasia <- output$Rescue [which (output$zooregions=="MALASIA")]
oceania <- output$Rescue [which (output$zooregions=="OCEANIA")]

X11()
par(mfrow=c(4,2))
hist(output$Rescue,nclass=100,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Global")
hist(neotropical,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Neotropic")
hist(afrotropical,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Afrotropic")
hist(malasia,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Indo-Malay")
hist(australia,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Australasia")
hist(neartic,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Nearctic")
hist(paleartic,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Palearctic")
hist(oceania,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Oceania")
