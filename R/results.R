setwd("C:/Users/Kelly/Desktop/Frontiers_GitHub")

out.spp <- read.table("Rescue_GLOBAL_RCP85_SN_50_without_DISP_8.5.txt",h=T)

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
