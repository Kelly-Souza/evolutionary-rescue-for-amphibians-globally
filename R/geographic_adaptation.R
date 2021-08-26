#########
setwd("C:/Users/Kelly/Desktop/Frontiers_GitHub")
rm(list=ls())


output <- read.table("Rescue_GLOBAL_RCP85_SN_50_without_DISP_8.5.txt",h=T)
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


prop_perda_Rich_Rescue_P <- (riqueza_rescue_P)/riqueza
raster_prop_perda_Rich_Rescue_P <- rasterFromXYZ(xyz =  cbind(coords, prop_perda_Rich_Rescue_P), res = 1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
x11()
plot(raster_prop_perda_Rich_Rescue_P)
writeRaster(raster_prop_perda_Rich_Rescue_P,"raster_prop_perda_Rich_Rescue_P.asc")

rangeSize <- colSums(PA)

class(PA.orig)
library(raster)

zooregions<- readShapePoly(fn = "zooreg_diss1.shp" )
plot(zooregions, add=T )


dados <- read.table("DADOS_ZOO.txt", h=T)
coords <- c(dados$X,dados$Y)
zooreg <- dados$ZooReg

x11()
mapas <- rasterFromXYZ(xyz =  cbind(coords, riqueza_lost), res = 1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
plot(mapas)
writeRaster(mapas,"raster_riqueza.asc")
writeRaster(mapas,"raster_rescue95.asc")
writeRaster(mapas,"raster_lost.asc")


regions<- readShapePoly("wwf_terr_ecos.shp")
regions@data$BIOME
plot(raster_prop_perda_Rich_Rescue)
plot(regions[regions@data$BIOME == 1,], add = T, border = "red")
View(regions@data)


rowSums(PA)
any(rowSums(rescue_PAM) > rowSums(rescue.P_PAM))

any(rowSums(rescue_PAM) > rowSums(PA))

any(rowSums(rescue.P_PAM) > rowSums(PA))



paleartic<- output$Rescue [which (output$zooregions=="PALEARTIC")]
neartic <- output$Rescue [which (output$zooregions=="NEARTIC")]
neotropical <- output$Rescue [which (output$zooregions=="NEOTROPICAL")]
australia <- output$Rescue [which (output$zooregions=="AUSTRALIA")]
afrotropical <- output$Rescue [which (output$zooregions=="AFROTROPICAL")]
malasia <- output$Rescue [which (output$zooregions=="MALASIA")]
oceania <- output$Rescue [which (output$zooregions=="OCEANIA")]

length(paleartic)
length(neartic)
length(neotropical)
length(australia)
length(afrotropical)
length(oceania)

head(rescue)
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

head(rescue)

#3x3
X11()
par(family="serif", cex.lab=1.8, cex.axis=1.8,cex.main = 3)
tiff("global.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(output$Rescue,nclass=100,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Global", ylim=c(0,2500))

dev.off()

tiff("global.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(neotropical,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Neotropic", ylim=c(0,1000))


tiff("global.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(afrotropical,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Afrotropic", ylim=c(0,1000))



hist(malasia,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Indo-Malay", ylim=c(0,1000))
hist(australia,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Australasia", ylim=c(0,200))
hist(neartic,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Nearctic", ylim=c(0,200))
hist(paleartic,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Palearctic", ylim=c(0,200))
hist(oceania,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Oceania", ylim=c(0,200))


#4X2
X11()
par(family="serif", cex.lab=1.2, cex.axis=1.2,cex.main = 1.8, mar=c(5,4,3,1))

tiff("Global.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(output$Rescue,nclass=100,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Global", ylim=c(0,2500))
dev.off()


tiff("Neotropic.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(neotropical,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Neotropic", ylim=c(0,1000))
dev.off()

tiff("Afrotropic.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(afrotropical,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Afrotropic", ylim=c(0,1000))
dev.off()

tiff("Indo-Malay.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(malasia,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Indo-Malay", ylim=c(0,1000))
dev.off()

tiff("Australasia.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(australia,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Australasia", ylim=c(0,200))
dev.off()


tiff("Nearctic.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(neartic,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Nearctic", ylim=c(0,200))
dev.off()

tiff("Palearctic.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(paleartic,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Palearctic", ylim=c(0,200))
dev.off()

tiff("Oceania.tif",  res = 300, width = 300/72*300, height = 300/72*300)
hist(oceania,nclass=50,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="",main="Oceania", ylim=c(0,200))
dev.off()
