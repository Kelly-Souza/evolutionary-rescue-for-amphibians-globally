setwd("C:/Users/Kelly/Desktop/Frontiers_GitHub")
rm(list=ls())
library(adegenet)
library(spdep)


PA.orig <-read.table("PAM_filtered.txt",h=T)

str(PA.orig)

head(PA.orig)
coords <-PA.orig[,1:2]
PA <-PA.orig[,3:7195]
remove(PA.orig) #release RAM memory

#write.table(coords,"Glob.Coords.txt")

ranges <-apply(PA,2,sum)
hist(ranges,col="grey",nclass=50,main="",xlab="Geographic ranges",ylab="number of species")
hist(log(ranges),col="grey",nclass=35,main="",xlab="Geographic ranges",ylab="number of species")
which(ranges==0) #ok

sum(ifelse(ranges == 1,1,0))/length(ranges)
sum(ifelse(ranges >= 10,1,0))/length(ranges)

mean.temp.pres <-read.table("mean_temp_present.txt",h=T)
mean.temp.fut <-read.table("mean_temp_future.txt",h=T)

head(mean.temp.pres)
head(mean.temp.fut)

riq <-rowSums(PA)
sum(ifelse(riq > 0,1,0))/length(riq)
riqP <-which(riq > 0)

#richness
s.value(coords,riq,xax = 1, yax = 2, method = c("greylevel"),csize = 0.1, cpoint = 0, pch = 20, clegend = 1, cneig = 1, xlim = NULL,
        ylim = NULL, grid = FALSE, addaxes = FALSE, cgrid = 0, include.origin = FALSE,
        origin = c(0,0),sub = "MAXENT", csub = 1.0, possub = "topleft", pixmap = NULL, contour = NULL,
        area = NULL, add.plot = FALSE)

#mapping spp1
#s.value(coords,PA[,1],xax = 1, yax = 2, method = c("greylevel"),csize = 0.1, cpoint = 0, pch = 20, clegend = 1, cneig = 1, xlim = NULL,
#ylim = NULL, grid = FALSE, addaxes = FALSE, cgrid = 0, include.origin = FALSE,
#origin = c(0,0),sub = "MAXENT", csub = 1.0, possub = "topleft", pixmap = NULL, contour = NULL,
#area = NULL, add.plot = FALSE)

#SET EMISSION SCENARIO
t0 <-mean.temp.pres[,1] # col1=CCSM
t1a <-mean.temp.fut[,1] # col1=CCSM Rcp26
t1b <-mean.temp.fut[,2] # col1=CCSM Rcp45
t1c <-mean.temp.fut[,3] # col1=CCSM Rcp60
t1d <-mean.temp.fut[,4] # col1=CCSM Rcp85

t1 <-t1d#for shifting the scenarios
#t1 <-t1a #for shifting the scenarios


#sd for spp with range > 10 (for bootstrapping)
sd.rg10 <-as.matrix(read.table("sd_rg10.txt",h=T))
median(sd.rg10)



#SIMULATIONS on
#Peak Shifts and Rescue across species

out.spp <-matrix(0,ncol(PA),10)
pb2 <- txtProgressBar(min = 0, max = ncol(PA), style = 3)


for(i in 1:ncol(PA)){
  
  setTxtProgressBar(pb2, i)
  
  range <-PA[,i]
  loc <-which(PA[,i]==1)
  mlat <-(sum(coords[,2]*range))/sum(range)
  mlong <-sum(coords[,1]*range)/sum(range)
  
  t0rg <-t0[which(range==1)]
  t1rg <-t1[which(range==1)]
  
  pk0 <-mean(t0[which(range==1)])
  pk1 <-mean(t1[which(range==1)])
  pk.shift <-pk1 - pk0
  
  #deviation in geographical space
  sd.rg <-sd(t0rg)
  sd.rg <-ifelse(sum(range) < 10,sd.rg <-sample(sd.rg10,size=1,replace=F),sd.rg)
  
  gen <-round(runif(1,40,50)) #generations
  
  
  #randomizations within each species
  sFST <-numeric()
  sd.imp <-numeric()
  hald <-numeric()
  rescue <-numeric()
  rescue.P <-numeric()
  freq.rg.genet <-numeric()
  freq.rg.plast <-numeric()
  
  for(j in 1:1000){
    
    #defining stochastic values
    h2 <-runif(1,0.20,0.4)
    sd.geo <-sd.rg
    
    sd <- sd.geo
    vT = sd*sd
    vA = vT*h2
    vE <-vT-vA
    
    #ss <-runif(1,80,120) #variation in strenght of selection; weak
    #ss <-runif(1,8,12) #variation in strenght of selection; medium
    #ss <-runif(1,1,1.2) #variation in strenght of selection; medium
    ss <-50 #more or less like Schneideri...
    w2 <-(sd.rg*sd.rg)*ss*h2
    Vs <- w2 + vE
    
    B <-runif(1,1.1,1.25)
    #sd <-sqrt((1 - fst) * (sd.geo*sd.geo))# for means the FST is not necessary...
    
    
    if(B*sqrt(w2/(vA+Vs)) < 1){
      TB=1.0001
    } else {
      TB <-B*sqrt(w2/(vA+Vs))
    }
    
    #Genetic adapt
    hald[j] <- (abs(pk.shift))/sd/gen
    kc <-(vA*(sqrt(2*log(TB)/(vA+Vs))))/sqrt(vT)
    rescue[j] <-ifelse(hald[j] > kc,0,1)
    
    #Chevin plasticity
    b <-runif(1,0,0.25)
    GT <- runif(1,1.5,3)
    kp <-sqrt(((2*log(B)*(1/Vs))/GT)*((h2*vT)/(1 - b)))
    rescue.P[j] <-ifelse(hald[j] > kp,0,1)
    
    
  }
  
  #Outputs
  out.spp[i,1] <-sum(range)
  out.spp[i,2] <-sd.rg
  out.spp[i,3] <-sd.rg
  out.spp[i,4] <-mlat
  out.spp[i,5] <-mlong
  out.spp[i,6] <-mean(t0[which(range==1)])
  out.spp[i,7] <-mean(t1[which(range==1)])
  out.spp[i,8] <-mean(hald)
  out.spp[i,9] <-mean(rescue)
  out.spp[i,10] <-mean(rescue.P)
 
}
head(out.spp)
colnames(out.spp) <-c("Range","sd.phy","sd.rg","Latitude", "Longitude","meanT0","meanT1",
                      "Haldane","P_Rescue","P_RescueP")
write.table(out.spp,"Rescue_GLOBAL_RCP85_SN_50_without_DISP_8.5.txt")
coorsdm <- out.spp[,4:5]
colnames(coorsdm) <- c("lat", "long")
write.table(coorsdm,"coordsm.txt")
#out.spp <-read.table("Rescue_GLOBAL_RCP85_SN_50",h=T)


sum(ifelse(out.spp[,1] >= 10,1,0))
sum(ifelse(out.spp[,1] >= 10,1,0))/nrow(out.spp)


hist(out.spp[,2],nclass=100,col="grey",xlab="Standard deviations (phylogeny)",main="")

x11()
par(family="serif", cex.lab=1.55, cex.axis=1.55)
hist(out.spp[,3],nclass=100,col="grey",xlab="Standard deviations ºC (range)",main="",ylab="Frequency of species" )

hist(out.spp[,7]-out.spp[,6],nclass=100,col="grey",xlab="Peak shift",main="",ylab="Frequency of species")

x11()
hist(out.spp[,7]-out.spp[,6],nclass=100,col="grey",xlab="Peak shift (ºC)",main="",ylab="Frequency of species", xlim = c(1,10))

hist(out.spp[,8],nclass=100,col="grey",xlab="Haldanes",ylab="Frequency of species",main="",xlim=c(0,1.0))
hist(out.spp[,9],nclass=100,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="")
hist(out.spp[,10],nclass=50,col="grey",xlab="Probability of Rescue (H > MSER-P)",ylab="Frequency of species",main="")

#rescue genetic vs. plasticity
plot(out.spp[,9],out.spp[,10],cex=0.5,pch=16,xlab="Probability of Rescue (Genetic) - MSER",ylab="Probabilit of rescue (Plasticity) - MSER-P")
abline(a=0,b=1,lwd=2,lty=2,col="black")

median(out.spp[,3])
median(out.spp[,7]-out.spp[,6])
median(out.spp[,8]) #haldMedian

plot(out.spp[,3],out.spp[,2]) #plot sds...
par(mfrow=c(1,1))
#for large range spp only...
plot(out.spp[out.spp[,1] >10,3],out.spp[out.spp[,1] >10,2]) #plot sds...
hist(out.spp[out.spp[,1] >10,9],nclass=100,col="grey",xlab="Probability of Rescue (H > MSER)",ylab="Frequency of species",main="")
hist(out.spp[out.spp[,1] >10,10],nclass=50,col="grey",xlab="Probability of Rescue (H > MSER-P)",ylab="Frequency of species",main="")


#persistence...
sum(ifelse(out.spp[,9] < 0.05,1,0)) / nrow(out.spp)
sum(ifelse(out.spp[,9] > 0.95,1,0)) / nrow(out.spp)

sum(ifelse(out.spp[,10] < 0.05,1,0)) / nrow(out.spp)
sum(ifelse(out.spp[,10] > 0.95,1,0)) / nrow(out.spp)

median(out.spp[,9])
median(out.spp[,10])



sum(ifelse(out.spp[,9] <= 0.025,1,0)) / nrow(out.spp)
sum(ifelse(out.spp[,9] >= 0.975,1,0)) / nrow(out.spp)

sum(ifelse(out.spp[,10] < 0.2,1,0)) / nrow(out.spp)
sum(ifelse(out.spp[,10] > 0.8,1,0)) / nrow(out.spp)

#rescue and range latitudinal position
plot(out.spp[,1],out.spp[,9],cex=1,pch=16,xlab="Geographic range",ylab="Probabilit of rescued")
plot(out.spp[,4],out.spp[,9],cex=1,pch=16,xlab="Latitudinal Midpoint",ylab="Probability of range rescued")
plot(out.spp[,5],out.spp[,9],cex=1,pch=16,xlab="Longitudinal Midpoint",ylab="Probability of trailing edge rescued")

#rescue and range latitudinal position
plot(out.spp[,1],out.spp[ ,10],cex=1,pch=16,xlab="Geographic range",ylab="Probabilit of rescued")
plot(out.spp[,4],out.spp[,10],cex=1,pch=16,xlab="Latitudinal Midpoint",ylab="Probability of range rescued")
plot(out.spp[,5],out.spp[,10],cex=1,pch=16,xlab="Longitudinal Midpoint",ylab="Probability of trailing edge rescued")


summary(lm(out.spp[,9]~log(out.spp[,1])))

summary(lm(out.spp[,9]~log(out.spp[,1])+out.spp[,2]+out.spp[,4]+out.spp[,5]))

max(out.spp[,7]-out.spp[,6])
min(out.spp[,7]-out.spp[,6])
median(out.spp[,7]-out.spp[,6]) # 3.75
median(out.spp[,8])#0.046
sum(ifelse(out.spp[,9] > 0.05,1,0)) / nrow(out.spp) # hald>0.05 - 55%
sum(ifelse(out.spp[,9] > 0.95,1,0)) / nrow(out.spp) #12%
sum(ifelse(out.spp[,10] > 0.05,1,0)) / nrow(out.spp) #79%
sum(ifelse(out.spp[,10] > 0.95,1,0)) / nrow(out.spp) #44%
sum(ifelse(out.spp[,9] > 0.95,1,0)) / nrow(out.spp)  #12%


rescue <- read.table("Rescue_GLOBAL_RCP85_SN_50_without_DISP_8.5_zooregions.txt",h=T)

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
