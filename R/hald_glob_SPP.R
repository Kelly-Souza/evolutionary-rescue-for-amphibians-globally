setwd("C:/Users/Kelly/Desktop/Frontiers_GitHub")

PA.orig <-read.table("PAM_filtered.txt",h=T)
coords <-PA.orig[,1:2]
PA <-PA.orig[,3:7195]
remove(PA.orig) #release RAM memory

write.table(coords,"Glob.Coords.txt")

ranges <-apply(PA,2,sum)
hist(ranges,col="grey",nclass=50,main="",xlab="Geographic ranges",ylab="number of species")
hist(log(ranges),col="grey",nclass=35,main="",xlab="Geographic ranges",ylab="number of species")

sum(ifelse(ranges == 1,1,0))/length(ranges)
sum(ifelse(ranges >= 10,1,0))/length(ranges)

mean.temp.pres <-read.table("mean_temp_present.txt",h=T)
mean.temp.fut <-read.table("mean_temp_future.txt",h=T)

riq <-rowSums(PA)
sum(ifelse(riq > 0,1,0))/length(riq)
riqP <-which(riq > 0)


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
