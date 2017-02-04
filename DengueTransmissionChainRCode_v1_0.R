##########################################################
## Code for manuscript entitled "Dengue diversity across spatial and temporal scales: local structure and the impact of host population size"
## Author: Henrik Salje
## Date: 2 February 2014
##########################################################

##########################################################
## Libraries
##########################################################
library(spatstat)


##########################################################
## Bring in datasets
##########################################################

load("gen.dat") ## R object for genetic data called gen.dat that consists of a list of four sublists, one for each serotype. Each sublist has the following components:
# $N 		- How many sequences of that serotype there are
# $dmat 	- An N x N matrix  - each cell represents the spatial distance between case i and case j (where there are n sequences available for that serotype) (diagonal is NA)
# $tmat  	- An N x N matrix - each cell represents the time (in days) between case i and case j (diagonal is NA)
# $gmat  	- An N x N matrix - each cell represents the total evolutionary time (in days) between that pair of viruses (diagonal is NA)
# $ttotmrca  	- An N x N matrix - each cell represents the time (in days) to the most recent common ancestor for the earlier of the two cases (diagonal is NA)
# dat 		- A dataframe with N rows comprising:
	# $date 		- Date (in decimal years - e.g., 2000.09 for 1 Feb 2000)
	# $serotype 	- The serotype of each case
	# $X 			- X coordinate of each case (UTM)
	# $Y 			- Y coordinate of each case (UTM)

load("ser.dat") ## R object for serotype data called ser.dat that consists of a list with the following components:
# $N 		- How many cases there are
# $dmat - An N x N matrix  - each cell represents the spatial distance between case i and case j (where there are n cases available for that serotype) (diagonal is NA)
# $tmat  - An N x N matrix - each cell represents the time (in days) between case i and case j (diagonal is NA)
# $smat  - An N x N matrix - each cell represents whether the serotype is the same (coded as 1) or not (coded as 0) (diagonal is NA)
# dat 		- A dataframe with N rows comprising:
	# $date 		- Date (in decimal years)
	# $serotype 	- The serotype of each case
	# $X 			- X coordinate of each case (UTM)
	# $Y 			- Y coordinate of each case (UTM)

load("regional.gen.dat") ## R object for genetic data called regional.gen.dat that consists of a list of four sublists, one for each serotype. Each sublist has the following components:
# $N 		- How many sequences of that serotype there are
# $dmat 	- An N x N matrix  - each cell represents the spatial distance between case i and case j (where there are n sequences available for that serotype) (diagonal is NA) - in this example, pairs where one is outside Bangkok is coded NA
# $tmat  	- An N x N matrix - each cell represents the time (in days) between case i and case j (diagonal is NA)
# $ttotmrca  	- An N x N matrix - each cell represents the time (in days) to the most recent common ancestor for the earlier of the two cases (diagonal is NA)
# dat 		- A dataframe with N rows comprising:
	# $DistrictID 	- District ID (e.g., 1 if Bangkok, 2 if Lampang etc., NA if not Thailand)
	# $CountryID 	- Country ID (e.g., 1 if Thailand, 2 if not Thailand)



##########################################################
## Calculate median spatial distance between pairs based on genetic separation (Basis of Figure 1G)
##########################################################

gdist.max<-seq(0,20*365,30)  	## Max genetic distances
gdist.window<-365					## Size of genetic distance window
gdist.min<-gdist.max-gdist.window ## Min genetic distances
gdist.min[which(gdist.min<0)]<-0
gdist.mid<-(gdist.min+gdist.max)/2 ## Midpoint of genetic distances

med.dist<-list()
for (i in 1:4){med.dist[[i]]<-list()}
for (ii in 1:4){
	for (i in 1:length(gdist.max)){
		tmp<-as.numeric((gen.dat[[ii]]$gmat<gdist.max[i])*(gen.dat[[ii]]$gmat>=gdist.min[i])*gen.dat[[ii]]$dmat)
		med.dist[[ii]][[i]]<-tmp[which(tmp>0)]
	}
}
med.dists.out<-rep(NaN,length(gdist.max))
for (i in 1:length(gdist.seq)){
	med.dists.out[i]<-median(c(med.dist[[1]][[i]],med.dist[[2]][[i]],med.dist[[3]][[i]],med.dist[[4]][[i]]))
}
plot(gdist.mid,med.dists.out/1000,pch=20,xlab="",ylab="",las=1,frame.plot=F,cex=0.75,cex.lab=0.75,cex.axis=0.75)



##########################################################
## Probability same chain and number of chains by distance based on genetic data (Basis of Figures 2 and 3)
##########################################################

max.tdist<-365/2		# Max time between cases (days)
max.mrca.dist<-365/2	# Max evolutionary time to MRCA (days)
sdists<-seq(0,5000,200) # Spatial distances at which to calculate probability

numerator<-denominator<-array(NaN,c(4,length(sdists)))
for (ii in 1:4){
	for (j in 1:length(sdists)){
		numerator[ii,j]<-sum((gen.dat[[ii]]$dmat<sdists[j])*(gen.dat[[ii]]$ttotmrca<max.mrca.dist)*(gen.dat[[ii]]$tmat<max.tdist),na.rm=T)
		denominator[ii,j]<-sum((gen.dat[[ii]]$dmat<sdists[j])*(gen.dat[[ii]]$tmat<max.tdist),na.rm=T)
	}
}
probSameChainHomotypic<-apply(numerator,c(2),sum)/apply(denominator,c(2),sum)  # Sum over serotypes - gives probability same chain for homotypic cases

## Adjust for underlying distribution of serotypes
tmp.num<-tmp.den<-rep(0,length(sdists))
for (ii in 1:4){
	tmat<-(abs(outer(gen.dat[[ii]]$dat$date,ser.dat$dat$date,"-"))*365)<max.tdist
	tmat[which(tmat==0)]<-NaN
	smat<-outer(gen.dat[[ii]]$dat$serotype,ser.dat$dat$serotype,"==")
	smat[which(smat==0)]<-NaN
	dmat<-crossdist(gen.dat[[ii]]$dat$X,gen.dat[[ii]]$dat$Y,ser.dat$dat$X,ser.dat$dat$Y)

	dt<-cumsum(hist(dmat*tmat,breaks=c(0,sdists,1e10),plot=F)$counts)-gen.dat[[ii]]$N #Subtract gen.dat[[ii]]$N as those in Gen database also appear in Ser database (avoids self comparison)
	dst<-cumsum(hist(dmat*tmat*smat,breaks=c(0,sdists,1e10),plot=F)$counts)-gen.dat[[ii]]$N 

	tmp.num<-tmp.num+dst[-(length(sdists)+1)]
	tmp.den<-tmp.den+dt[-(length(sdists)+1)]
}
probSameChain<-probSameChainHomotypic*tmp.num/tmp.den

plot(sdists,probSameChain,pch=20,xlab="Distance (m)",ylab="Probability same chain",las=1,frame.plot=F,cex=0.75,cex.lab=0.75,cex.axis=0.75,ylim=c(0,1))

plot(sdists,1/probSameChain,pch=20,xlab="Distance (m)",ylab="# Chains",las=1,frame.plot=F,cex=0.75,cex.lab=0.75,cex.axis=0.75,log="xy")




##########################################################
## Probability same chain and number of chains by distance based on serotype data only (Basis of Figures 2 and 3)
##########################################################

max.tdist<-365/2		# Max time between cases (days)
max.mrca.dist<-365/2	# Max evolutionary time to MRCA (days)
sdists<-seq(0,5000,200) # Spatial distances at which to calculate probability
dist.no.spat.dep<-10000 #Distance where no spatial dependence between cases

p0<-sum((ser.dat$dmat>dist.no.spat.dep)*ser.dat$smat*(ser.dat$tmat<max.tdist),na.rm=T)/sum((ser.dat$dmat>dist.no.spat.dep)*(ser.dat$tmat<max.tdist),na.rm=T)  # Estimate of underlying probability of unrelated cases nevertheless being homotypic

tmat<-ser.dat$tmat<max.tdist
tmat[which(tmat==0)]<-NA
smat<-ser.dat$smat
smat[which(smat==0)]<-NA

a<-cumsum(hist(tmat*ser.dat$dmat*smat,breaks=c(0,sdists,1e10),plot=F)$counts)[-(length(sdists)+1)]
b<-cumsum(hist(tmat*ser.dat$dmat,breaks=c(0,sdists,1e10),plot=F)$counts)[-(length(sdists)+1)]

pX<-a/b
probSameChainSerOnly<-(pX-p0)/(1-p0)

plot(sdists,probSameChainSerOnly,pch=20,xlab="Distance (m)",ylab="Probability same chain (ser only)",las=1,frame.plot=F,cex=0.75,cex.lab=0.75,cex.axis=0.75,ylim=c(0,1),col=2)

plot(sdists,1/probSameChainSerOnly,pch=20,xlab="Distance (m)",ylab="# Chains (ser only)",las=1,frame.plot=F,cex=0.75,cex.lab=0.75,cex.axis=0.75,log="xy",col=2)






##########################################################
## Relative risk of having MRCA within specified range for pairs of cases a certain distance apart relative to reference population (Basis of Figure 4)
##########################################################

max.tdist<-365/2		# Max time between cases (days)
max.mrca.dists<-365*c(0.5,2,5,10)	# Max evolutionary times to MRCA (days)
min.mrca.dists<-c(0,max.mrca.dists[-length(max.mrca.dists)])	# Min evolutionary times to MRCA (days)
max.sdists<-c(500,1000,5000,10000,1e10) # Max spatial distances at which to calculate RR within BKK
min.sdists<-c(0,max.sdists[-length(max.sdists)]) # Max spatial distances at which to calculate RR within BKK
ddist.comp<-10000 # Reference group (pairs separated by greater than this distance within Bangkok)

withinBKK.num<-withinBKK.den<-array(NaN,c(length(max.mrca.dists),length(max.sdists),4)) # Within Bangkok
withinBKK.ref.num<-withinBKK.ref.den<-array(NaN,c(length(max.mrca.dists),4)) # Within Bangkok (same ref for all dists)
intraProv.num<-intraProv.den<-intraProv.ref.num<-intraProv.ref.den<-array(NaN,c(length(max.mrca.dists),4)) # Intra-province (within Thailand)
interProv.num<-interProv.den<-interProv.ref.num<-interProv.ref.den<-array(NaN,c(length(max.mrca.dists),4)) # Inter-province (within Thailand)
SEAsia.num<-SEAsia.den<-SEAsia.ref.num<-SEAsia.ref.den<-array(NaN,c(length(max.mrca.dists),4)) # Across SE Asia
for (i in 1:length(max.mrca.dists)){
	for (j in 1:4){
		withinT<-(regional.gen.dat[[j]]$tmat<=max.tdist)
		ref.withinT<-(regional.gen.dat[[j]]$dmat>ddist.comp)*withinT # Pairs sick within 6 months and distal Bangkok loacted (reference)
		withinBKK.ref.num[i,j]<-sum(ref.withinT*(regional.gen.dat[[j]]$ttotmrca<max.mrca.dists[i])*(regional.gen.dat[[j]]$ttotmrca>=min.mrca.dists[i]),na.rm=T)
		withinBKK.ref.den[i,j]<-sum(ref.withinT,na.rm=T)
		
		# Within Bangkok
		for (k in 1:length(max.sdists)){
			withinDT<-(regional.gen.dat[[j]]$dmat<=max.sdists[k])*(regional.gen.dat[[j]]$dmat>min.sdists[k])*withinT
			
			withinBKK.num[i,k,j]<-sum(withinDT*(regional.gen.dat[[j]]$ttotmrca<(max.mrca.dists[i]))*(regional.gen.dat[[j]]$ttotmrca>=(min.mrca.dists[i])),na.rm=T)
			withinBKK.den[i,k,j]<-sum(withinDT,na.rm=T)
		}

		# Intra-province	
		intra.prov.mat<-outer(regional.gen.dat[[j]]$dat$DistrictID,regional.gen.dat[[j]]$dat$DistrictID,"==")
		diag(intra.prov.mat)<-NA
		ind.BKK<-which(regional.gen.dat[[j]]$dat$DistrictID==1)  # BKK cases
		ind<-which(regional.gen.dat[[j]]$dat$DistrictID!=1)  # Non-BKK cases
		intraProv.num[i,j]<-sum((withinT*intra.prov.mat*(regional.gen.dat[[j]]$ttotmrca<max.mrca.dists[i])*(regional.gen.dat[[j]]$ttotmrca>=min.mrca.dists[i]))[ind,ind],na.rm=T)
		intraProv.den[i,j]<-sum((withinT*intra.prov.mat)[ind,ind],na.rm=T)

		ref.ind<-which(rowSums(withinT[ind.BKK,ind],na.rm=T)>0) # In reference, only include viruses with at least one case within 6 months of case in numerator
		intraProv.ref.num[i,j]<-sum((ref.withinT*(regional.gen.dat[[j]]$ttotmrca<max.mrca.dists[i])*(regional.gen.dat[[j]]$ttotmrca>=min.mrca.dists[i]))[ref.ind,ref.ind],na.rm=T)
		intraProv.ref.den[i,j]<-sum(ref.withinT[ref.ind,ref.ind],na.rm=T)

		# Inter-province
		ind<-which(regional.gen.dat[[j]]$dat$DistrictID!=1)  # Not in BKK
		interProv.num[i,j]<-sum((withinT*(regional.gen.dat[[j]]$ttotmrca<max.mrca.dists[i])*(regional.gen.dat[[j]]$ttotmrca>=min.mrca.dists[i]))[ind.BKK,ind],na.rm=T)
		interProv.den[i,j]<-sum(withinT[ind.BKK,ind],na.rm=T)
		interProv.ref.num[i,j]<-sum((ref.withinT*(regional.gen.dat[[j]]$ttotmrca<max.mrca.dists[i])*(regional.gen.dat[[j]]$ttotmrca>=min.mrca.dists[i]))[ref.ind,ref.ind],na.rm=T)
		interProv.ref.den[i,j]<-sum(ref.withinT[ref.ind,ref.ind],na.rm=T)

		# SE-Asia - BKK as source
		ind<-which(regional.gen.dat[[j]]$dat$CountryID==2)  # SE Asia (Not Thailand)
		SEAsia.num[i,j]<-sum((withinT*(regional.gen.dat[[j]]$ttotmrca<max.mrca.dists[i])*(regional.gen.dat[[j]]$ttotmrca>=min.mrca.dists[i]))[ind.BKK,ind],na.rm=T)
		SEAsia.den[i,j]<-sum(withinT[ind.BKK,ind],na.rm=T)
		ref.ind<-which(rowSums(withinT[ind.BKK,ind],na.rm=T)>0) 
		SEAsia.ref.num[i,j]<-sum((ref.withinT*(regional.gen.dat[[j]]$ttotmrca<max.mrca.dists[i])*(regional.gen.dat[[j]]$ttotmrca>=min.mrca.dists[i]))[ref.ind,ref.ind],na.rm=T)
		SEAsia.ref.den[i,j]<-sum(ref.withinT[ref.ind,ref.ind],na.rm=T)
	}
}

#Sum over serotypes
# Within Bangkok
withinBKK.prop.withMRCA<-apply(withinBKK.num,c(1,2),sum)/apply(withinBKK.den,c(1,2),sum)
withinBKK.prop.withMRCA.ref<-apply(withinBKK.ref.num,c(1),sum)/apply(withinBKK.ref.den,c(1),sum)
withinBKK.RR<-sweep(withinBKK.prop.withMRCA,1,withinBKK.prop.withMRCA.ref,"/")

# Intra-province
intraProv.prop.withMRCA<-apply(intraProv.num,c(1),sum)/apply(intraProv.den,c(1),sum)
intraProv.prop.withMRCA.ref<-apply(intraProv.ref.num,c(1),sum)/apply(intraProv.ref.den,c(1),sum)
intraProv.RR<-intraProv.prop.withMRCA/intraProv.prop.withMRCA.ref

# Inter-province
interProv.prop.withMRCA<-apply(interProv.num,c(1),sum)/apply(interProv.den,c(1),sum)
interProv.prop.withMRCA.ref<-apply(interProv.ref.num,c(1),sum)/apply(interProv.ref.den,c(1),sum)
interProv.RR<-interProv.prop.withMRCA/interProv.prop.withMRCA.ref

# SE-Asia 
SEAsia.prop.withMRCA<-apply(SEAsia.num,c(1),sum)/apply(SEAsia.den,c(1),sum)
SEAsia.prop.withMRCA.ref<-apply(SEAsia.ref.num,c(1),sum)/apply(SEAsia.ref.den,c(1),sum)
SEAsia.RR<-SEAsia.prop.withMRCA/SEAsia.prop.withMRCA.ref

ref=3
plot(c(withinBKK.RR[ref,],intraProv.RR[ref],interProv.RR[ref],SEAsia.RR[ref]),log="y")










