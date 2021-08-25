#########
vetorID_rescue <- rescue$ID [which (rescue$Rescue_Save_95==1)]

rescue_PAM <- PA[, vetorID_rescue]
str(rescue_PAM)

vetorID_rescueP <- rescue$ID [which (rescue$RescueP_Save_95==1)]
rescue.P_PAM <- PA[, vetorID_rescueP]

PAM_rescue <- cbind(coords,rescue_PAM)
ncol(PAM_rescue)
PAM_rescueP <- cbind(coords,rescue.P_PAM)
ncol(PAM_rescueP)
str(PAM_rescueP)

write.table(PAM_rescue, "PAM_rescue.txt")
write.table(PAM_rescueP, "PAM_rescueP.txt")

library(raster)
x11()
par(mfrow=c(1,1))
#)
raster_riq <- rasterFromXYZ(xyz =  cbind(coords, riq), res = 1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

plot(raster_riq)
writeRaster(raster_riq,"raster_riq.asc")

raster_rescue<- rasterFromXYZ(xyz =  cbind(coords, rowSums(rescue_PAM)), res = 1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
plot(raster_rescue)
writeRaster(raster_rescue,"raster_riq_rescue.95.asc")
nrow(coords)
nrow(rescue_PAM)

x11()
raster_rescueP <- rasterFromXYZ(xyz =  cbind(coords, rowSums(rescue.P_PAM)), res = 1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
plot(raster_rescueP)
writeRaster(raster_rescueP,"raster_riq_rescueP.95.asc")

library(maptools)
data("wrld_simpl")
plot(wrld_simpl, add = T, border = "black")

riqueza <- rowSums(PA)
riqueza_rescue <- rowSums(rescue_PAM)
riqueza_rescue <- ifelse(test = riqueza > 0, yes = riqueza_rescue, no = NA)

riqueza_rescue_P <- rowSums(rescue.P_PAM)
riqueza_rescue_P <- ifelse(test = riqueza > 0, yes = riqueza_rescue_P, no = NA)

riqueza <- ifelse(test = riqueza > 0, yes = riqueza, no = NA)

prop_perda_Rich_Rescue <- (riqueza_rescue)/riqueza
raster_prop_perda_Rich_Rescue <- rasterFromXYZ(xyz =  cbind(coords, prop_perda_Rich_Rescue), res = 1, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
x11()
plot(raster_prop_perda_Rich_Rescue)
writeRaster(raster_prop_perda_Rich_Rescue,"raster_prop_perda_Rich_Rescue.asc")

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



paleartic<- rescue$Rescue [which (rescue$zooregions=="PALEARTIC")]
neartic <- rescue$Rescue [which (rescue$zooregions=="NEARTIC")]
neotropical <- rescue$Rescue [which (rescue$zooregions=="NEOTROPICAL")]
australia <- rescue$Rescue [which (rescue$zooregions=="AUSTRALIA")]
afrotropical <- rescue$Rescue [which (rescue$zooregions=="AFROTROPICAL")]
malasia <- rescue$Rescue [which (rescue$zooregions=="MALASIA")]
oceania <- rescue$Rescue [which (rescue$zooregions=="OCEANIA")]

length(paleartic)
length(neartic)
length(neotropical)
length(australia)
length(afrotropical)
length(oceania)

head(rescue)
X11()
par(mfrow=c(4,2))
hist(rescue$Rescue,nclass=100,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Global")
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
hist(rescue$Rescue,nclass=100,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Global", ylim=c(0,2500))

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
hist(rescue$Rescue,nclass=100,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="Global", ylim=c(0,2500))
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