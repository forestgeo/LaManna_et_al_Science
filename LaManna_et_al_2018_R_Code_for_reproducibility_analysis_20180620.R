### ForestGEO analysis of CNDD and latitudinal diversity relationships - LaManna et al. 2017 & 2018, Science
### R Code for Technical Response analyses using data formatted for reproducibility analysis
### Code by: J. A. LaManna, updated: 06-20-2018


##########################################################################################################
###
### DATA AND CODE USE DISCLAIMER:
###   This code and accompanying data are intended only for the reproduction of analyses presented in 
###   LaManna et al. (2018, Science 360, technical responses eaar3824 & eaar5245).  Anyone seeking to use 
###   this code and/or data for any modified and/or new analysis and/or any publications must first obtain 
###   permission from PIs of all 24 Smithsonian ForestGEO plots included in this dataset.  Failure to do so 
###   constitutes violation of data sharing agreements put in place by the CTFS-Smithsonian ForestGEO network.
###
##########################################################################################################


##########################################################################################################
### Load data
### Data list contains quadrat location and DBH for all main stems in each plot, as well as GX GY locations 
### for adults of all species with sample size large enough to be included in CNDD analyses.  
### GX GY locations of saplings are not required to reproduce analyses, only GX GY locations of quadrat 
### centers, which are provided in the data list as well.  
### Code is included to classify individual stems as adults and saplings and to calculate maximum height
### for each species based on CTFS allometric equations (from the CTFSRPackage and Chave 2005 Oecologia). 


load("CTFS_LaManna_et_al_2018_data_for_reproducibility_analysis_20180620.RData")	
load("CTFS_LaManna_et_al_2018_plotdata_for_reproducibility_analysis_20180620.RData")
load("CTFS_LaManna_et_al_spplist_for_reproducibility_analysis_20180620.RData")


##########################################################################################################
### Load/install required R packages

if (!require("doBy", character.only=T, quietly=T)) {
    install.packages("doBy", repos="http://cran.us.r-project.org")}
    library("doBy", character.only=T)
if (!require("reshape", character.only=T, quietly=T)) {
    install.packages("reshape", repos="http://cran.us.r-project.org")}
    library("reshape", character.only=T)
if (!require("vegan", character.only=T, quietly=T)) {
    install.packages("vegan", repos="http://cran.us.r-project.org")}
    library("vegan", character.only=T)
if (!require("spatstat", character.only=T, quietly=T)) {
    install.packages("spatstat", repos="http://cran.us.r-project.org")}
    library("spatstat", character.only=T)
if (!require("gnm", character.only=T, quietly=T)) {
    install.packages("gnm", repos="http://cran.us.r-project.org")}
    library("gnm", character.only=T)


##########################################################################################################
### Load functions

### Funciton to classify saplings and adults
adult.sapling.classify = function(test) {
test$sapling=c(0)
dbhclass5=c()
dbhclass2=c()
test$sapling[which(test$dbh<10)]=1
test$sapling=as.factor(test$sapling)
sapclass=tapply(test$sapling,test$sp,table)
for(i in 1:length(sapclass)) {if(sapclass[[i]][1]/sum(sapclass[[i]])<0.20) {
	test$sapling[which(test$sp==names(sapclass[i]) & test$dbh>5)]=0
	dbhclass5=append(dbhclass5,names(sapclass[i]))
}}
sapclass=tapply(test$sapling,test$sp,table)
for(i in 1:length(sapclass)) {if(sapclass[[i]][1]/sum(sapclass[[i]])<0.20) {
	test$sapling[which(test$sp==names(sapclass[i]) & test$dbh>2)]=0
	dbhclass2=append(dbhclass2,names(sapclass[i]))
}}
dbhclass5=dbhclass5[-which(dbhclass5 %in% dbhclass2==T)]
return(test$sapling)
}

### Function to calculate maximum height for each species
### This function is the 'predht.asym' function from the CTFSRPackage (from 'Chave.AGB' function used in Chave 2005 Oecologia)
### Default height parameters are from Chave et al 2003 on BCI biomass (see http://ctfs.si.edu/Public/CTFSRPackage/index.php/web/topics/biomass~slash~biomass.CTFSdb.r/Chave.AGB)
predht.asym=function(dbh,param)
{
 if(is.null(dim(param)))
  {
   ymax=param[1]
   a=param[2]
   b=param[3]
  }
 else
  {
   ymax=param[,1]
   a=param[,2]
   b=param[,3]
  }

 return(ymax*(1-exp(-a*dbh^b)))
}
ctfs.predict.asym = function(test) {
htparam=c(41.7, 0.057, 0.748) 	# Default CTFS height parameters from BCI data (see http://ctfs.si.edu/Public/CTFSRPackage/index.php/web/topics/biomass~slash~biomass.CTFSdb.r/predht.asym)
d=tapply(test$dbh, test$sp, max) 
max.height = predht.asym(dbh = d, param = htparam)
max.heights = data.frame("sp" = names(max.height), "max.height" = max.height)
return(max.heights)
}

### Clark's 2Dt seed dispersal function with p = 1
seed = function(r, alpha=100) {(1/(pi*alpha*((1+((r^2)/alpha))^2)))}	

### Function to calculate distance-weighted adult abundance (using Clark's 2dT dispersal kernel with alpha = exp(6), or mean dispersal distance of ~30 m)
dist.weighted.abund.function = function(test, quads.mat) {
adultlocs = test[which(test$sapling == 0),c("gx","gy")]
xdist = crossdist(adultlocs[,1], adultlocs[,2], quads.mat[,3], quads.mat[,4])
xdist.weighted.alphae6 = seed(xdist, alpha = exp(6)) * (1/seed(0.1, alpha = exp(6)))
adults.dist.weighted = data.frame(quadrat = quads.mat[,2], sp = test$sp[1], adult6 = apply(xdist.weighted.alphae6, 2, sum))
return(adults.dist.weighted)
}

### Ricker model
Ricker <- function(x,Had,Hsap){
	list(predictors =	list(r = 1, CNDD = 1, HNDDad = 1, HNDDsap = 1),
	variables = list(substitute(x), substitute(Had), substitute(Hsap)),
      term = function(predictors, variables) {
          pred <- paste("(", variables[1],")*exp(", predictors[1], 
				")*exp(", predictors[2], "*", variables[1], 
				")*exp(", predictors[3], "*", variables[2],
				")*exp(", predictors[4], "*", variables[3],")+0.0001", sep = "")
         })
}
class(Ricker ) <- "nonlin"

### Function to fit model with nonlinear Ricker model (distance-weighted approach)
fit.ricker.cndd = function(data){
x=data$adult; y = data$sap; Hsap = data$Hsap; Had = data$Had
return(tryCatch(gnm(y~-1+Ricker(x,Had,Hsap),family = quasipoisson(link="identity")), error=function(e) NULL))
}

### Function to fit model with nonlinear Ricker model (original offset approach)
fit.ricker.cndd.offset = function(data, offset = 0.1){
x=data$adult; y = data$sap; Hsap = data$Hsap; Had = data$Had
x[which(x == 0 & y > 0)] = x[which(x == 0 & y > 0)] + offset	# Add offset only to quads with saplings but zero adults so they are not excluded
return(tryCatch(gnm(y~-1+Ricker(x,Had,Hsap),family = quasipoisson(link="identity")), error=function(e) NULL))
}

### Function to fit model with nonlinear Ricker model (applying offset to all quadrats)
fit.ricker.cndd.offset.all = function(data, offset = 0.1){
x=data$adult; y = data$sap; Hsap = data$Hsap; Had = data$Had
x = x + offset							# Add offset to all quads so quads with saplings but zero adults are not excluded
return(tryCatch(gnm(y~-1+Ricker(x,Had,Hsap),family = quasipoisson(link="identity")), error=function(e) NULL))
}

### Function to calculate relationship between CNDD and species abundance for a forest plot
CNDD.across.spp=function(test) {
model=lm(CNDD~log(ba),weights=(1/(CNDD.se)),data=test)
int=summary(model)$coef[1,1]
int.se=summary(model)$coef[1,2]
beta=summary(model)$coef[2,1]
beta.se=summary(model)$coef[2,2]
ests=c(int,int.se,beta,beta.se)
names(ests)=c("int","int.se","beta","beta.se")
return(ests)
}

### Function to prepare data from each forest plot for CNDD analysis
CNDD.analysis.data.prep = function(test) {

# Name objects from data list and rename species and plot columns
data2 = test[[1]]
quads.mat = test[[2]]
data2$sp = data2$spp.code
data2$plot = data2$plot.code


# Classify individuals as adults or saplings
data2$sapling = adult.sapling.classify(data2)

# Calculate maximum height for each species
max.heights = ctfs.predict.asym(data2)
data2 = data2[,-which(colnames(data2) == "max.height")]
data2 = merge(data2, max.heights, by = "sp", all.x = T)

## Create species by site matrices for CNDD analysis
data2$value=c(1)
adult.spp.forcast=summaryBy(value~quadrat+sp,FUN=sum,data=data2[which(data2$sapling==0),],na.rm=T)
adult.spp=cast(adult.spp.forcast,quadrat~sp,value='value.sum',fun.aggregate=sum,fill=0)
saps.spp.forcast=summaryBy(value~quadrat+sp,FUN=sum,data=data2[which(data2$sapling==1),],na.rm=T)
saps.spp=cast(saps.spp.forcast,quadrat~sp,value='value.sum',fun.aggregate=sum,fill=0)

## Restrict data to quadrats that share saplings and adults
adult.spp.adndd=adult.spp
goodplots=intersect(adult.spp.adndd$quadrat,saps.spp$quadrat)
adult.spp.adndd=adult.spp.adndd[which(adult.spp.adndd$quadrat %in% goodplots==T),]
saps.adndd=saps.spp[which(saps.spp$quadrat %in% goodplots==T),]
adult.spp.adndd=adult.spp.adndd[order(adult.spp.adndd$quadrat),]
saps.adndd=saps.adndd[order(saps.adndd$quadrat),]

## Remove species with no individuals in common quadrats
if(0 %in% colSums(adult.spp.adndd[,-1])==T) {adult.spp.adndd=adult.spp.adndd[,-c(which(colSums(adult.spp.adndd[,-1])==0)+1)]}
if(0 %in% colSums(saps.adndd[,-1])==T) {saps.adndd=saps.adndd[,-c(which(colSums(saps.adndd[,-1])==0)+1)]}

## Restrict data to species with both sapling and adult data
goodspp=intersect(colnames(adult.spp.adndd),colnames(saps.adndd))
adult.spp.adndd=adult.spp.adndd[,which(colnames(adult.spp.adndd) %in% goodspp==T)]
saps.adndd=saps.adndd[,which(colnames(saps.adndd) %in% goodspp==T)]

## Double-check that all quadrats have data
goodplots=intersect(adult.spp.adndd$quadrat,saps.spp$quadrat)
adult.spp.adndd=adult.spp.adndd[which(adult.spp.adndd$quadrat %in% goodplots==T),]
saps.adndd=saps.adndd[which(saps.adndd$quadrat %in% goodplots==T),]
adult.spp.adndd=adult.spp.adndd[order(adult.spp.adndd$quadrat),]
saps.adndd=saps.adndd[order(saps.adndd$quadrat),]

# Arrange dataset for CNDD analysis (1 row = number conspecific and heterospecific individuals of one species in one quadrat)
adults.ndd=adult.spp.adndd[,-1]
saps.ndd=saps.adndd[,-1]
conspp.adult=adults.ndd[,1]
for(i in 2:dim(adults.ndd)[2]) {conspp.adult=c(conspp.adult,adults.ndd[,i])}
conspp.sap=saps.ndd[,1]
for(i in 2:dim(saps.ndd)[2]) {conspp.sap=c(conspp.sap,saps.ndd[,i])}
sp=rep(colnames(saps.ndd)[1],times=dim(saps.ndd)[1])
for(i in 2:dim(saps.ndd)[2]) {sp=c(sp,rep(colnames(saps.ndd)[i],times=dim(saps.ndd)[1]))}
quadrat=as.character(adult.spp.adndd$quadrat)
for(i in 2:dim(adults.ndd)[2]) {quadrat=c(quadrat,as.character(adult.spp.adndd$quadrat))}
heterospp.adult=rowSums(adults.ndd)-adults.ndd[,1]
for(i in 2:dim(adults.ndd)[2]) {heterospp.adult=c(heterospp.adult,rowSums(adults.ndd)-adults.ndd[,i])}
heterospp.sap=rowSums(saps.ndd)-saps.ndd[,1]
for(i in 2:dim(saps.ndd)[2]) {heterospp.sap=c(heterospp.sap,rowSums(saps.ndd)-saps.ndd[,i])}
full=data.frame(conspp.adult,conspp.sap,heterospp.adult,heterospp.sap)
full$quadrat=quadrat
full$sp=sp

# Calculate distance-weighted adult abundance for each species in each quadrat
goodspp = unique(data2$sp[!is.na(data2$order)])
data3 = data2[which(data2$sp %in% goodspp == T),]
data3$sp = as.character(data3$sp)
splist = split(data3, data3$sp)
adults.dist.weighted = lapply(splist, dist.weighted.abund.function, quads.mat = quads.mat)
adults.dist.weighted = data.frame(do.call('rbind', adults.dist.weighted))
adults.dist.weighted$quadrat = as.character(adults.dist.weighted$quadrat)
names(adults.dist.weighted) = c("quadrat", "sp", "adult.dist.wght.alphae6")
full2 = merge(full, adults.dist.weighted, by = c("sp", "quadrat"), all.x = T)
full2$plot = c(data2$plot.code[1])
full = full2[,c("plot", "sp", "quadrat", "conspp.adult", "conspp.sap", "heterospp.adult", "heterospp.sap", "adult.dist.wght.alphae6")]
return(full)
}


basal.area.function = function(test) {
data2 = test[[1]]
data2$sp = data2$spp.code
data2$plot = data2$plot.code
data2$ba=(((data2$dbh/100)/2)^2)*pi  		# Calculate basal area for each individual
batest = summaryBy(ba ~ sp, FUN = sum, data = data2)
names(batest)=c("sp", "uncorrected.ba")
plotsize=(length(unique(data2$quadrat))*400)/10000
batest$ba = batest$uncorrected.ba / plotsize
batest$plot = data2$plot[1]
return(batest)
}



##########################################################################################################
### Analysis (data preparation)

allfull.list = lapply(ctfs.data.list, CNDD.analysis.data.prep)
allfull = data.frame(do.call("rbind", allfull.list))

ba.list = lapply(ctfs.data.list, basal.area.function)
spp.ba = data.frame(do.call("rbind", ba.list))

names(spplist)[which(names(spplist) == "spp.code")] = "sp"
names(spplist)[which(names(spplist) == "plot.code")] = "plot"
spplist2 = merge(spplist, spp.ba, by = c("plot", "sp"), all.x = T)
spplist = spplist2[,c("plot", "sp", "order", "included.spp", "ba")]

allfull2 = merge(allfull, spplist, by = c("plot", "sp"), all.x = T)
allfull = allfull2[,c("plot", "sp", "quadrat", "conspp.adult", "conspp.sap", "adult.dist.wght.alphae6", "heterospp.adult", "heterospp.sap", "order", "included.spp", "ba")]

save(allfull, file = "CTFS_LaManna_et_al_Full_Quad_dataset_20180620.RData")


##########################################################################################################
### Analysis (CNDD estimation with Ricker model)

### Prepare data for Ricker analyses
data = allfull
data.summary = summaryBy(as.numeric(conspp.adult > 0) + as.numeric(conspp.sap > 0) ~ sp, FUN = sum, data = data)
names(data.summary)[c(2:3)]=c("adultquads", "sapquads")
data1 = merge(data, data.summary, by = "sp", all.x = T)
data1 = data1[order(data1$order),]
data.backup = data1





### Ricker CNDD Analysis
data1 = data.backup

set.seed(15)


# Remove species with saplings or adults occupying fewer than 9 quadrats
testq = data1[which(data1$adultquads>9),]
testq = testq[which(testq$sapquads>9),]

# Remove liana and bamboo species
testq = testq[!is.na(testq$order),]


#### Loop through plots and species 
testq = testq[order(testq$order),]
#testq = testq[order(testq$plot, testq$sp),]


est.list = list()
studyplots = unique(testq$plot)
residlist = list()
length(unique(testq$sp))

for(i in 1:length(studyplots)) {
testq2 = testq[which(testq$plot == studyplots[i]),]
sap = testq2$conspp.sap
Had = testq2$heterospp.adult
Hsap = testq2$heterospp.sap

# ***IMPORTANT*** Only run one of the two lines below either for the distance-weighted adult abundances (LaManna et al. Science 2018)
# 		  or the quadrat-based adult abundances (original method) 
adult = testq2$adult.dist.wght.alphae6		### Distance-weighted conspecific adult abundances
# adult = testq2$conspp.adult			### Conspecific adult abundances in each quadrat

spp = as.numeric(as.factor(testq2$sp))
n.groups = length(unique(spp))
n = length(sap)
spnames = data.frame(levels(as.factor(testq2$sp)),c(1:length(unique(testq2$sp))))
names(spnames) = c("sp","code")
data = data.frame("adult" = adult, "sap" = sap, "Had" = Had, "Hsap" = Hsap, "spp" = spp)
qlist = split(data, spp)

# ***IMPORTANT*** Only run one of the three lines below either for the distance-weighted adult abundance approach (LaManna et al. Science 2018), 
#		  the original offset approach (LaManna et al. 2017), or the offset applied to all quadrats approach (LaManna et al. 2018)
fits3 = lapply(qlist, fit.ricker.cndd)					### Distance-weighted conspecific adult abundance approach
# fits3 = lapply(qlist, fit.ricker.cndd.offset, offset = 0.1)		### Original offset approach (Add offset only to quads with saplings but zero adults so they are not excluded)
# fits3 = lapply(qlist, fit.ricker.cndd.offset.all, offset = 0.1)	### Offset applied to all quadrats so quads with saplings but zero adults are not excluded


fits2=list()
splist=c()
r.cor.test=c()
for(j in 1:length(fits3)) {
  possibleError <- tryCatch(summary(fits3[[j]])$coef, error=function(e) e  )
  if(!inherits(possibleError, "error")){best.pow.optim=c(t(summary(fits3[[j]])$coef[,1:2]))
names(best.pow.optim)=paste(rep(rownames(summary(fits3[[j]])$coef[,1:2]),each=2),rep(c("",".se"),times=3),sep="")
splist=append(splist,as.character(spnames$sp[j]))
r.cor.test=append(r.cor.test,cor.test(fits3[[j]]$data$y,fitted(fits3[[j]]))$estimate)
fits2[[j]]=best.pow.optim
}}

fits=data.frame(do.call("rbind", fits2))
fits$plot=c(studyplots[i])
fits$sp=splist
fits$r.cor.test=r.cor.test
fits4 = merge(fits, spplist, by = c("plot", "sp"), all.x = T)
est.list[[i]]=fits4
}






##########################################################################################################
### Examine CNDD estimates

est.list.backup=est.list

est.list=est.list.backup

# Function to remove poorly estimated speces
fReliableEsts=function(test) {
test2=test[!is.na(test$CNDD.se),]
test2=test2[which(test2$CNDD.se<100),]
return(test2)}
est.list=lapply(est.list,fReliableEsts)
est.list.test=data.frame(do.call("rbind", est.list))
dim(est.list.test)
est.list.test = est.list.test[order(est.list.test$plot, est.list.test$sp),]
est.list = split(est.list.test, est.list.test$plot)
save(est.list.test, file = "CTFS_LaManna_et_al_Species_CNDD_Estimates.RData")


# Calculate weighted mean and median estimates for each forest plot
fmedianCNDD=function(test) {
return(median(test$CNDD))}
fmedianHNDDad=function(test) {
return(median(test$HNDDad))}
fmedianHNDDsap=function(test) {
return(median(test$HNDDsap))}
fmedianR=function(test) {
return(median(test$r))}
fweighted.mean.R=function(test) {
return(summary(lm(test$r~1,weights=(1/(test$r.se))))$coef[1])}
fweighted.mean=function(test) {
return(summary(lm(test$CNDD~1,weights=(1/(test$CNDD.se))))$coef[1])}
fweighted.mean.HNDDad=function(test) {
return(summary(lm(test$HNDDad~1,weights=(1/(test$HNDDad.se))))$coef[1])}
fweighted.mean.HNDDsap=function(test) {
return(summary(lm(test$HNDDsap~1,weights=(1/(test$HNDDsap.se))))$coef[1])}

medianR=unlist(lapply(est.list,fmedianR))
medianCNDD=unlist(lapply(est.list,fmedianCNDD))
medianHNDDad=unlist(lapply(est.list,fmedianHNDDad))
medianHNDDsap=unlist(lapply(est.list,fmedianHNDDsap))
meanR=unlist(lapply(est.list,fweighted.mean.R))
meanCNDD=unlist(lapply(est.list,fweighted.mean))
meanHNDDad=unlist(lapply(est.list,fweighted.mean.HNDDad))
meanHNDDsap=unlist(lapply(est.list,fweighted.mean.HNDDsap))




### Plot species rarefied richness against strength of CNDD across forest plots
plot(medianCNDD,plotdata$plotrarefy)
cor.test(medianCNDD,plotdata$plotrarefy,method="spearman")



### Calculation of relationship between CNDD and species abundance for each forest plot
sppcnddlist=lapply(est.list, CNDD.across.spp)
sppcnddlist=data.frame(do.call("rbind", sppcnddlist))
plot(abs(plotdata$lat), sppcnddlist$beta, xlab = "Absolute latitude", ylab = "Slope btw. CNDD and species abundance")
cor.test(abs(plotdata$lat), sppcnddlist$beta, method = "spearman")










####################################################################################################################################
####################################################################################################################################
### Null model analysis (maintain adult locations to preserve habitat specificity, distribute recruits with dispersal kernel)
### Null model recommended by Wiegand & Moloney (2014 book) for analyses involving adults and their recruits
####################################################################################################################################
####################################################################################################################################


###############################################################
### Load data

load("CTFS_LaManna_et_al_2018_data_for_reproducibility_analysis_20180620.RData")
load("CTFS_LaManna_et_al_2018_plotdata_for_reproducibility_analysis_20180620.RData")
load("CTFS_LaManna_et_al_spplist_for_reproducibility_analysis_20180620.RData")
load("CTFS_LaManna_et_al_Species_CNDD_Estimates.RData")
load("CTFS_LaManna_et_al_Full_Quad_dataset_20180620.RData")

###############################################################
### Load Functions

# Clark's 2Dt seed dispersal kernel
dispersal.fun = function(x, y, alpha, n) {
seedp = function(r, alpha=100) {(1/(pi*alpha*((1+((r^2)/alpha))^2)))*2*pi*r}	# Clark's 2Dt p = 1 (integrated in 2 polar dimensions)
disp.dists = seq(0, 10000, by = 1)
rho <- sample(disp.dists, size = n, replace = TRUE, prob = seedp(r = disp.dists, alpha = alpha))
theta <- runif(n, 0, 2*pi)
x1 <- (rho * cos(theta)) + x
y1 <- (rho * sin(theta)) + y
result = data.frame(x1, y1)
return(result)
}

dispersal.fun2 = function(loc, xlim = plotwidth, ylim = plotheight, alpha = 100) {	# Recruits disperse across plot edges in a torus
test = dispersal.fun(loc[1], loc[2], alpha = alpha, n = 1)
test[,1] = test[,1] %% xlim		# Torus
test[,2] = test[,2] %% ylim
return(data.frame(x1 = test[,1], y1 = test[,2]))
}

# Function from Thomson et al 2011 J. of Ecol. Relationship between max tree height and mean dispersal distance (plus observed error around regression fit)
mean.disp.dist = function(max.height) {10 ^ ((0.1975 + (0.936 * log10(max.height))) + rnorm(length(max.height), mean = 0, sd = 0.245))}		# Relationship between maximum tree height and mean dispersal distance from Thomson et al 2011 J. of Ecol.

# Function to calculate alpha from Clark's 2Dt model given mean dispersal distance
clark.alpha = function(mean.disp, shape = 1) {exp(log(mean.disp) - (0.5*log(pi) + log(gamma(shape - 1/2)) - log(2) - log(gamma(shape))))^2}	

# Function to assemble null communities with Clark's 2Dt dispersal kernel with 30 m mean dispersal for all species
null.assembly.30m.mean.dispersal = function(test, alpha = 100) {
alpha = clark.alpha(30)							# 30 m mean dispersal
adults = test[which(test$sapling == 0),]
saps = test[which(test$sapling == 1),]
sgx <- adults$gx					# Use to fix adult locs		
sgy <- adults$gy					# Use to fix adult locs	
adultlocs <- data.frame(x1 = sgx, y1 = sgy, quadrat = adults$quadrat)
adultlocs$sapling = 0
repro.adults <- sample(1:nrow(adults), size = nrow(saps), replace = T)
saplocs = apply(data.frame(sgx, sgy)[repro.adults,], 1, dispersal.fun2, xlim = plotwidth, ylim = plotheight, alpha = alpha)
saplocs = data.frame(do.call("rbind", saplocs))
saplocs$sapling = 1
saplocs$quadrat = c(999999)
simlocs = rbind(adultlocs, saplocs)
simlocs$sp = test$sp[1]
return(simlocs)
}

# Function to assemble null communities with Clark's 2Dt dispersal kernel with allometrically-scaled mean dispersal with error (from Thomson et al. 2011)
null.assembly.allometric.mean.dispersal = function(test, alpha = 100) {
# The following two lines of code use allometrically scaled dispersal kernel (based on a species' maximum height, see Chave 2005 Oecologia for allometric equations relating DBH to height)
mean.disp.distance = mean.disp.dist(test$max.height[1])		# For Thomson et al 2011 tree height-dispersal relationship
alpha = clark.alpha(mean.disp.distance)				# For Thomson et al 2011 tree height-dispersal relationship 
adults = test[which(test$sapling == 0),]
saps = test[which(test$sapling == 1),]
sgx <- adults$gx					# Use to fix adult locs		
sgy <- adults$gy					# Use to fix adult locs	
adultlocs <- data.frame(x1 = sgx, y1 = sgy, quadrat = adults$quadrat)
adultlocs$sapling = 0
repro.adults <- sample(1:nrow(adults), size = nrow(saps), replace = T)
saplocs = apply(data.frame(sgx, sgy)[repro.adults,], 1, dispersal.fun2, xlim = plotwidth, ylim = plotheight, alpha = alpha)
saplocs = data.frame(do.call("rbind", saplocs))
saplocs$sapling = 1
saplocs$quadrat = c(999999)
simlocs = rbind(adultlocs, saplocs)
simlocs$sp = test$sp[1]
return(simlocs)
}

# Function to fit Ricker model for null simulated data
Ricker <- function(x,Had,Hsap){
	list(predictors =	list(r = 1, CNDD = 1),
	variables = list(substitute(x)),
      term = function(predictors, variables) {
          pred <- paste("(", variables[1],")*exp(", predictors[1], 
				")*exp(", predictors[2], "*", variables[1],")+0.0001", sep = "")
         })
}
class(Ricker ) <- "nonlin"

# Function to fit Ricker model
fit.ricker.cndd = function(data){
x = data$adult; y = data$saps
return(tryCatch(gnm(y~-1+Ricker(x),family = quasipoisson(link="identity")), error=function(e) NULL))
}

quad.ID.cast <- function(simlocs, height = plotheight, width = plotwidth, by.meters = 20, sp = "XXXX") {
wdvec = as.numeric(cut(simlocs$x1, breaks = seq(0, width, by = by.meters)))
htvec = as.numeric(cut(simlocs$y1, breaks = seq(0, height, by = by.meters)))
wdvec = sprintf("%03d",wdvec)
htvec = sprintf("%03d",htvec)
simlocs$quad = paste(wdvec, htvec, sep = "")
simlocs$value = c(1)
adultcast = summaryBy(value ~ quad, keep.names = T, FUN = sum, data = simlocs[which(simlocs$sapling == 0),])
names(adultcast)[2] = "adults"
sapcast = summaryBy(value ~ quad, keep.names = T, FUN = sum, data = simlocs[which(simlocs$sapling == 1),])
names(sapcast)[2] = "saps"
simcast = merge(adultcast, sapcast, by = "quad", all = T)
simcast = merge(simcast, totquad, by = "quad", all = T)
simcast[is.na(simcast)] <- 0
adults.dist.weighted.sp = adults.dist.weighted3[which(adults.dist.weighted3$sp == sp),c("quad", "adult.dist")]
simcast2 = merge(simcast, adults.dist.weighted.sp, by = "quad")
simcast = data.frame(quad = simcast2$quad, adults = simcast2$adult.dist, saps = simcast2$saps)
return(simcast)
}


fits.summary = function(Rickerfits, sp = "XXXX") {
fits2 = list()
for(j in 1:length(Rickerfits)) {
  possibleError <- tryCatch(summary(Rickerfits[[j]])$coef, error=function(e) e  )
  if(!inherits(possibleError, "error")) {best.pow.optim = c(t(coef(summary(Rickerfits[[j]]))[,1:2]))
names(best.pow.optim) = paste(rep(rownames(coef(summary(Rickerfits[[j]]))[,1:2]), each = 2), rep(c("", ".se"), times = 1), sep = "")
fits2[[j]] = best.pow.optim
}}
fits = data.frame(do.call("rbind", fits2))
fits$sp = sp
return(fits)
}


null.model.30m.mean.dispersal = function(test, k = 2000) {
null.assemblies.list = replicate(k, null.assembly.30m.mean.dispersal(test), simplify = F)
null.assemblies.quad.20x20 = lapply(null.assemblies.list, quad.ID.cast, by = 20, sp = test$sp[1])
Rickerfits.20x20.dist.weight = lapply(null.assemblies.quad.20x20, fit.ricker.cndd)
mod.list = list(Rickerfits.20x20.dist.weight)
mod.list.results = lapply(mod.list, fits.summary, sp = test$sp[1])
names(mod.list.results) = c("Rickerfits.20x20.dist.weight")
return(mod.list.results)
}

null.model.allometric.mean.dispersal = function(test, k = 2000) {
null.assemblies.list = replicate(k, null.assembly.allometric.mean.dispersal(test), simplify = F)
null.assemblies.quad.20x20 = lapply(null.assemblies.list, quad.ID.cast, by = 20, sp = test$sp[1])
Rickerfits.20x20.dist.weight = lapply(null.assemblies.quad.20x20, fit.ricker.cndd)
mod.list = list(Rickerfits.20x20.dist.weight)
mod.list.results = lapply(mod.list, fits.summary, sp = test$sp[1])
names(mod.list.results) = c("Rickerfits.20x20.dist.weight")
return(mod.list.results)
}






#############################################################################################
### Null model analysis (30m mean dispersal parameter for all species)
### Note: This null model analysis can take 72-120 hours to complete (or more, depending on computing speed)
### Note: This null model analysis also requires completely running all code above this point
### For allometric-scaled dispersal null model, skip to next section below from this point


CNDD.null.model.results.30m.mean.dispersal = list()
for(z in 1:length(ctfs.data.list)){
test = ctfs.data.list[[z]]

# Name objects from data list and rename species and plot columns
data2 = test[[1]]
quads.mat = test[[2]]
data2$sp = data2$spp.code
data2$plot = data2$plot.code
plotname = data2$plot[1]

plotwidth =  test[[3]]	# Plot width
plotheight = test[[4]]	# Plot height

forestplot.20x20.dist.weight = est.list.test[which(est.list.test$plot == plotname),]
goodspp = unique(forestplot.20x20.dist.weight$sp)
data3 = data2[which(data2$sp %in% goodspp == T),]
data3$sp = as.character(data3$sp)
splist = split(data3, data3$sp)

wdvec = c(1:(plotwidth/20))
htvec = c(1:(plotheight/20))
grid = expand.grid(wdvec = wdvec, htvec = htvec)
grid$wdvec = sprintf("%03d", grid$wdvec)
grid$htvec = sprintf("%03d", grid$htvec)
totquad = data.frame(quad = paste(grid$wdvec, grid$htvec, sep = ""))

wdvecm = as.numeric(cut(quads.mat$qx, breaks = seq(0, plotwidth, by = 20)))
htvecm = as.numeric(cut(quads.mat$qy, breaks = seq(0, plotheight, by = 20)))
wdvecm = sprintf("%03d",wdvecm)
htvecm = sprintf("%03d",htvecm)
quads.mat$quad = paste(wdvecm, htvecm, sep = "")
adults.dist.weighted = allfull[which(allfull$plot == plotname),]
adults.dist.weighted2 = data.frame(quadrat = adults.dist.weighted$quadrat, sp = adults.dist.weighted$sp, adult.dist = adults.dist.weighted$adult.dist.wght.alphae6)
adults.dist.weighted3 = merge(adults.dist.weighted2, quads.mat, by = "quadrat")
adults.dist.weighted3 = adults.dist.weighted3[,c("quad", "sp", "adult.dist")]

# Run analysis and save results
set.seed(5064)
Ricker.cluster.null.models = lapply(splist, null.model.30m.mean.dispersal, k = 100)

# Summarize analysis
Rickersims.20x20.dist.weight = list()
for(j in 1:length(Ricker.cluster.null.models)) {
Rickersims.20x20.dist.weight[[j]] = Ricker.cluster.null.models[[j]][[1]]}
Rickersims.20x20.dist.weight = data.frame(do.call('rbind', Rickersims.20x20.dist.weight))

# Merge simulated results with observed results
forestplot.20x20.dist.weight = est.list.test[which(est.list.test$plot == plotname),]
forestplot.20x20.dist.weight = forestplot.20x20.dist.weight[,c("sp", "CNDD", "CNDD.se", "ba")]
names(Rickersims.20x20.dist.weight)[1:4] = c("r.sim", "r.se.sim", "CNDD.sim", "CNDD.se.sim")
Rickersims.20x20.dist.weight = merge(forestplot.20x20.dist.weight, Rickersims.20x20.dist.weight, by = "sp", all.x = T)

# turn these into list
Rickersims = list(Rickersims.20x20.dist.weight)
Model.estimates = list(forestplot.20x20.dist.weight)

obs.medianCNDD = c()
medianCNDD.es = c() 
medianCNDD.ses = c()
obs.rare.medianCNDD = c()
rare.medianCNDD.es = c()
rare.medianCNDD.ses = c()
obs.cnddspp.slope = c()
CNDDSPP.es = c()
CNDDSPP.ses = c()

for(j in 1:length(Rickersims)) {
its = as.numeric(table(Rickersims[[j]]$sp))
its2 = c()
for(i in 1:length(its)) {
its2 = append(its2, 1:its[i])
}
Rickersims[[j]]$it = its2
Rickersims[[j]] = Rickersims[[j]][!is.na(Rickersims[[j]]$CNDD.se.sim),]
if(max(Rickersims[[j]]$CNDD.se.sim) > 100) {Rickersims[[j]] = Rickersims[[j]][-which(Rickersims[[j]]$CNDD.se.sim > 100),]}
test.its = split(Rickersims[[j]], Rickersims[[j]]$it)
test.sp = split(Rickersims[[j]], Rickersims[[j]]$sp)

# Median CNDD
medianCNDD.sim = unlist(lapply(test.its, function(x) {median(x$CNDD.sim)}))
obs.medianCNDD[j] = median(Model.estimates[[j]]$CNDD)
medianCNDD.es[j] = (obs.medianCNDD[j] - mean(medianCNDD.sim)) 
medianCNDD.ses[j] = (obs.medianCNDD[j] - mean(medianCNDD.sim)) / sd(medianCNDD.sim) 

# Rare species median CNDD
rarecutoff = 0.1
rare.medianCNDD.sim = unlist(lapply(test.its, function(x) {median(x$CNDD.sim[which(x$ba <= rarecutoff)])}))
rare.medianCNDD.sim.mean = mean(rare.medianCNDD.sim)
rare.medianCNDD.sim.sd = sd(rare.medianCNDD.sim)
obs.rare.medianCNDD[j] = median(Model.estimates[[j]]$CNDD[which(Model.estimates[[j]]$ba <= rarecutoff)])
rare.medianCNDD.es[j] = (obs.rare.medianCNDD[j] - rare.medianCNDD.sim.mean) 
rare.medianCNDD.ses[j] = (obs.rare.medianCNDD[j] - rare.medianCNDD.sim.mean) / rare.medianCNDD.sim.sd 

# CNDDSPP
CNDD.sp.abund.null.function = function(test) {
model = lm(CNDD.sim ~ log(ba), weights = (1/CNDD.se.sim), data = test)
int=summary(model)$coef[1,1]
int.se=summary(model)$coef[1,2]
beta=summary(model)$coef[2,1]
beta.se=summary(model)$coef[2,2]
ests=c(int,int.se,beta,beta.se)
names(ests)=c("int","int.se","beta","beta.se")
return(ests)}
CNDD.sim.sppcndd = lapply(test.its, CNDD.sp.abund.null.function)
CNDD.sim.sppcndd = data.frame(do.call("rbind", CNDD.sim.sppcndd))
CNDD.sim.sppcndd.mean = apply(CNDD.sim.sppcndd, 2, mean)
CNDD.sim.sppcndd.sd = apply(CNDD.sim.sppcndd, 2, sd)
obs.cnddspp = lm(CNDD ~ log(ba), weights = (1/CNDD.se), data = Model.estimates[[j]])
obs.cnddspp.slope[j] = summary(obs.cnddspp)$coef[2,1]
CNDDSPP.es[j] = (obs.cnddspp.slope[j] - CNDD.sim.sppcndd.mean[3]) 
CNDDSPP.ses[j] = (obs.cnddspp.slope[j] - CNDD.sim.sppcndd.mean[3]) / CNDD.sim.sppcndd.sd[3]
}

# Create summary object
forestplot.wiegand.summary = data.frame(obs.medianCNDD, medianCNDD.es, medianCNDD.ses, 
	obs.cnddspp.slope, CNDDSPP.es, CNDDSPP.ses, obs.rare.medianCNDD, rare.medianCNDD.es, rare.medianCNDD.ses)
forestplot.wiegand.summary$plot = plotname

# return summary object
null.results.file = list(forestplot.wiegand.summary, Ricker.cluster.null.models, Rickersims)
CNDD.null.model.results.30m.mean.dispersal[[z]] = null.results.file
}

save(CNDD.null.model.results.30m.mean.dispersal, file = "CNDD_null_model_results_30m_mean_dispersal.RData")




# Examine null model results

null.results = list()
for(i in 1:length(CNDD.null.model.results.30m.mean.dispersal)) {
null.results[[i]] = CNDD.null.model.results.30m.mean.dispersal[[i]][[1]]}
null.results = data.frame(do.call("rbind", null.results))
null.results$medianCNDD.sim = null.results$obs.medianCNDD - null.results$medianCNDD.es
null.results$CNDDSPP.sim = null.results$obs.cnddspp.slope - null.results$CNDDSPP.es
null.results$rare.medianCNDD.sim = null.results$obs.rare.medianCNDD - null.results$rare.medianCNDD.es

summary((null.results$obs.medianCNDD - null.results$medianCNDD.sim) - null.results$medianCNDD.es)			# Should be all zeros
summary((null.results$obs.cnddspp.slope - null.results$CNDDSPP.sim) - null.results$CNDDSPP.es)				# Should be all zeros
summary((null.results$obs.rare.medianCNDD - null.results$rare.medianCNDD.sim) - null.results$rare.medianCNDD.es)	# Should be all zeros

save(null.results, file = "CTFS_LaManna_et_al_Null_Model_results_30m_mean_dispersal.RData")



### Figures

load("CTFS_LaManna_et_al_Null_Model_results_30m_mean_dispersal.RData")

# Plot observed, simulated, effect size, and std. effect size relationship between median CNDD and species rarefied richness across plots
par(mfrow = c(2,2))
plot(null.results$obs.medianCNDD, plotdata$plotrarefy, xlim = c(min(null.results$obs.medianCNDD), max(null.results$obs.medianCNDD)), 
	xlab = "Observed median CNDD", ylab = "Plot rarefied species richness", las = 1)
abline(v = 0, lty = 2)
legend("bottomleft", legend = c(paste("r = ",round(cor.test(null.results$obs.medianCNDD, plotdata$plotrarefy,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(null.results$obs.medianCNDD, plotdata$plotrarefy, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(null.results$medianCNDD.sim, plotdata$plotrarefy, xlim = c(min(null.results$obs.medianCNDD), max(null.results$obs.medianCNDD)), 
	xlab = "Simulated median CNDD", ylab = "Plot rarefied species richness", las = 1)
abline(v = 0, lty = 2)
plot(null.results$medianCNDD.es, plotdata$plotrarefy, 
	xlab = "Median CNDD ES", ylab = "Plot rarefied species richness", las = 1)
abline(v = 0, lty = 2)
legend("bottomleft", legend = c(paste("r = ",round(cor.test(null.results$medianCNDD.es, plotdata$plotrarefy,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(null.results$medianCNDD.es, plotdata$plotrarefy, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(null.results$medianCNDD.ses, plotdata$plotrarefy, 
	xlab = "Median CNDD SES", ylab = "Plot rarefied species richness", las = 1)
abline(v = 0, lty = 2)
legend("bottomleft", legend = c(paste("r = ",round(cor.test(null.results$medianCNDD.ses, plotdata$plotrarefy,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(null.results$medianCNDD.ses, plotdata$plotrarefy, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
title("Null model results for median CNDD (30 m mean dispersal across species)", outer = T, line = -2)



# Plot observed, simulated, effect size, and std. effect size for relationship between CNDD-species-abundance slopes and absolute latitude across plots
par(mfrow = c(2,2))
plot(abs(plotdata$lat),null.results$obs.cnddspp.slope, ylim = c(min(null.results$CNDDSPP.sim), max(null.results$obs.cnddspp.slope)), 
	xlab = "Absolute latitude", ylab = "Obs. slope btw. CNDD & spp. abund.", las = 1)
abline(h = 0, lty = 2)
legend("topright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$obs.cnddspp.slope,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$obs.cnddspp.slope, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(abs(plotdata$lat),null.results$CNDDSPP.sim, ylim = c(min(null.results$CNDDSPP.sim), max(null.results$obs.cnddspp.slope)), 
	xlab = "Absolute latitude", ylab = "Sim. slope btw. CNDD & spp. abund.", las = 1)
abline(h = 0, lty = 2)
plot(abs(plotdata$lat),null.results$CNDDSPP.es, 
	xlab = "Absolute latitude", ylab = "Slope btw. CNDD & spp. abund. ES", las = 1)
abline(h = 0, lty = 2)
legend("topright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$CNDDSPP.es,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$CNDDSPP.es, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(abs(plotdata$lat),null.results$CNDDSPP.ses, 
	xlab = "Absolute latitude", ylab = "Slope btw. CNDD & spp. abund. SES", las = 1)
abline(h = 0, lty = 2)
legend("topright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$CNDDSPP.ses,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$CNDDSPP.ses, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
title("Null model results for CNDD-spp. abund. (30 m mean dispersal across spp.)", outer = T, line = -2)



# Plot observed, simulated, effect size, and std. effect size relationship between median CNDD of rare species and absolute latitude across plots
par(mfrow = c(2,2))
plot(abs(plotdata$lat),null.results$obs.rare.medianCNDD, ylim = c(min(null.results$obs.rare.medianCNDD), max(null.results$obs.rare.medianCNDD)), 
	xlab = "Observed median CNDD", ylab = "Observed median for rare spp. CNDD", las = 1)
abline(h = 0, lty = 2)
legend("bottomright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$obs.rare.medianCNDD,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$obs.rare.medianCNDD, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(abs(plotdata$lat),null.results$rare.medianCNDD.sim, ylim = c(min(null.results$obs.rare.medianCNDD), max(null.results$obs.cnddspp.slope)), 
	xlab = "Simulated median CNDD", ylab = "Simulated median for rare spp. CNDD", las = 1)
abline(h = 0, lty = 2)
plot(abs(plotdata$lat),null.results$rare.medianCNDD.es, 
	xlab = "Median CNDD ES", ylab = "Median CNDD for rare spp. ES", las = 1)
abline(h = 0, lty = 2)
legend("bottomright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$rare.medianCNDD.es,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$rare.medianCNDD.es, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(abs(plotdata$lat),null.results$rare.medianCNDD.ses, 
	xlab = "Median CNDD SES", ylab = "Median CNDD for rare spp. SES", las = 1)
abline(h = 0, lty = 2)
legend("bottomright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$rare.medianCNDD.ses,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$rare.medianCNDD.ses, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
title("Null model results for median CNDD for rare spp. (30 m mean dispersal across spp.)", outer = T, line = -2)










#############################################################################################
### Null model analysis (allometrically-scaled mean dispersal parameter with stochastic error)
### Note: This null model analysis can take 72-120 hours to complete (or more, depending on computing speed)
### Note: This null model analysis also requires completely running all code above the 30-m null model analysis, including all functions required for null model analysis (above)


CNDD.null.model.results.allometric.mean.dispersal = list()
for(z in 1:length(ctfs.data.list)){
test = ctfs.data.list[[z]]

# Name objects from data list and rename species and plot columns
data2 = test[[1]]
quads.mat = test[[2]]
data2$sp = data2$spp.code
data2$plot = data2$plot.code
plotname = data2$plot[1]

plotwidth =  test[[3]]	# Plot width
plotheight = test[[4]]	# Plot height

forestplot.20x20.dist.weight = est.list.test[which(est.list.test$plot == plotname),]
goodspp = unique(forestplot.20x20.dist.weight$sp)
data3 = data2[which(data2$sp %in% goodspp == T),]
data3$sp = as.character(data3$sp)
splist = split(data3, data3$sp)

wdvec = c(1:(plotwidth/20))
htvec = c(1:(plotheight/20))
grid = expand.grid(wdvec = wdvec, htvec = htvec)
grid$wdvec = sprintf("%03d", grid$wdvec)
grid$htvec = sprintf("%03d", grid$htvec)
totquad = data.frame(quad = paste(grid$wdvec, grid$htvec, sep = ""))

wdvecm = as.numeric(cut(quads.mat$qx, breaks = seq(0, plotwidth, by = 20)))
htvecm = as.numeric(cut(quads.mat$qy, breaks = seq(0, plotheight, by = 20)))
wdvecm = sprintf("%03d",wdvecm)
htvecm = sprintf("%03d",htvecm)
quads.mat$quad = paste(wdvecm, htvecm, sep = "")
adults.dist.weighted = allfull[which(allfull$plot == plotname),]
adults.dist.weighted2 = data.frame(quadrat = adults.dist.weighted$quadrat, sp = adults.dist.weighted$sp, adult.dist = adults.dist.weighted$adult.dist.wght.alphae6)
adults.dist.weighted3 = merge(adults.dist.weighted2, quads.mat, by = "quadrat")
adults.dist.weighted3 = adults.dist.weighted3[,c("quad", "sp", "adult.dist")]

# Run analysis and save results
set.seed(5064)
Ricker.cluster.null.models = lapply(splist, null.model.allometric.mean.dispersal, k = 100)

# Summarize analysis
Rickersims.20x20.dist.weight = list()
for(j in 1:length(Ricker.cluster.null.models)) {
Rickersims.20x20.dist.weight[[j]] = Ricker.cluster.null.models[[j]][[1]]}
Rickersims.20x20.dist.weight = data.frame(do.call('rbind', Rickersims.20x20.dist.weight))

# Merge 20x20 offset = 0.1 sims with observed results
forestplot.20x20.dist.weight = est.list.test[which(est.list.test$plot == plotname),]
forestplot.20x20.dist.weight = forestplot.20x20.dist.weight[,c("sp", "CNDD", "CNDD.se", "ba")]
names(Rickersims.20x20.dist.weight)[1:4] = c("r.sim", "r.se.sim", "CNDD.sim", "CNDD.se.sim")
Rickersims.20x20.dist.weight = merge(forestplot.20x20.dist.weight, Rickersims.20x20.dist.weight, by = "sp", all.x = T)

# turn these into list
Rickersims = list(Rickersims.20x20.dist.weight)
Model.estimates = list(forestplot.20x20.dist.weight)

obs.medianCNDD = c()
medianCNDD.es = c() 
medianCNDD.ses = c()
obs.rare.medianCNDD = c()
rare.medianCNDD.es = c()
rare.medianCNDD.ses = c()
obs.cnddspp.slope = c()
CNDDSPP.es = c()
CNDDSPP.ses = c()

for(j in 1:length(Rickersims)) {
its = as.numeric(table(Rickersims[[j]]$sp))
its2 = c()
for(i in 1:length(its)) {
its2 = append(its2, 1:its[i])
}
Rickersims[[j]]$it = its2
Rickersims[[j]] = Rickersims[[j]][!is.na(Rickersims[[j]]$CNDD.se.sim),]
if(max(Rickersims[[j]]$CNDD.se.sim) > 100) {Rickersims[[j]] = Rickersims[[j]][-which(Rickersims[[j]]$CNDD.se.sim > 100),]}
test.its = split(Rickersims[[j]], Rickersims[[j]]$it)
test.sp = split(Rickersims[[j]], Rickersims[[j]]$sp)

# Median CNDD
medianCNDD.sim = unlist(lapply(test.its, function(x) {median(x$CNDD.sim)}))
obs.medianCNDD[j] = median(Model.estimates[[j]]$CNDD)
medianCNDD.es[j] = (obs.medianCNDD[j] - mean(medianCNDD.sim)) 
medianCNDD.ses[j] = (obs.medianCNDD[j] - mean(medianCNDD.sim)) / sd(medianCNDD.sim) 

# Rare species median CNDD
rarecutoff = 0.1
rare.medianCNDD.sim = unlist(lapply(test.its, function(x) {median(x$CNDD.sim[which(x$ba <= rarecutoff)])}))
rare.medianCNDD.sim.mean = mean(rare.medianCNDD.sim)
rare.medianCNDD.sim.sd = sd(rare.medianCNDD.sim)
obs.rare.medianCNDD[j] = median(Model.estimates[[j]]$CNDD[which(Model.estimates[[j]]$ba <= rarecutoff)])
rare.medianCNDD.es[j] = (obs.rare.medianCNDD[j] - rare.medianCNDD.sim.mean) 
rare.medianCNDD.ses[j] = (obs.rare.medianCNDD[j] - rare.medianCNDD.sim.mean) / rare.medianCNDD.sim.sd 

# CNDDSPP
CNDD.sp.abund.null.function = function(test) {
model = lm(CNDD.sim ~ log(ba), weights = (1/CNDD.se.sim), data = test)
int=summary(model)$coef[1,1]
int.se=summary(model)$coef[1,2]
beta=summary(model)$coef[2,1]
beta.se=summary(model)$coef[2,2]
ests=c(int,int.se,beta,beta.se)
names(ests)=c("int","int.se","beta","beta.se")
return(ests)}
CNDD.sim.sppcndd = lapply(test.its, CNDD.sp.abund.null.function)
CNDD.sim.sppcndd = data.frame(do.call("rbind", CNDD.sim.sppcndd))
CNDD.sim.sppcndd.mean = apply(CNDD.sim.sppcndd, 2, mean)
CNDD.sim.sppcndd.sd = apply(CNDD.sim.sppcndd, 2, sd)
obs.cnddspp = lm(CNDD ~ log(ba), weights = (1/CNDD.se), data = Model.estimates[[j]])
obs.cnddspp.slope[j] = summary(obs.cnddspp)$coef[2,1]
CNDDSPP.es[j] = (obs.cnddspp.slope[j] - CNDD.sim.sppcndd.mean[3]) 
CNDDSPP.ses[j] = (obs.cnddspp.slope[j] - CNDD.sim.sppcndd.mean[3]) / CNDD.sim.sppcndd.sd[3]
}

# Create summary object
forestplot.wiegand.summary = data.frame(obs.medianCNDD, medianCNDD.es, medianCNDD.ses, 
	obs.cnddspp.slope, CNDDSPP.es, CNDDSPP.ses, obs.rare.medianCNDD, rare.medianCNDD.es, rare.medianCNDD.ses)
forestplot.wiegand.summary$plot = plotname

# return summary object
null.results.file = list(forestplot.wiegand.summary, Ricker.cluster.null.models, Rickersims)
CNDD.null.model.results.allometric.mean.dispersal[[z]] = null.results.file
}


save(CNDD.null.model.results.allometric.mean.dispersal, file = "CNDD_null_model_results_allometric_mean_dispersal.RData")


# Examine null model results

null.results = list()
for(i in 1:length(CNDD.null.model.results.allometric.mean.dispersal)) {
null.results[[i]] = CNDD.null.model.results.allometric.mean.dispersal[[i]][[1]]}
null.results = data.frame(do.call("rbind", null.results))
null.results$medianCNDD.sim = null.results$obs.medianCNDD - null.results$medianCNDD.es
null.results$CNDDSPP.sim = null.results$obs.cnddspp.slope - null.results$CNDDSPP.es
null.results$rare.medianCNDD.sim = null.results$obs.rare.medianCNDD - null.results$rare.medianCNDD.es

summary((null.results$obs.medianCNDD - null.results$medianCNDD.sim) - null.results$medianCNDD.es)			# Should be all zeros
summary((null.results$obs.cnddspp.slope - null.results$CNDDSPP.sim) - null.results$CNDDSPP.es)				# Should be all zeros
summary((null.results$obs.rare.medianCNDD - null.results$rare.medianCNDD.sim) - null.results$rare.medianCNDD.es)	# Should be all zeros

save(null.results, file = "CTFS_LaManna_et_al_Null_Model_results_allometric_mean_dispersal.RData")



### Figures

load("CTFS_LaManna_et_al_Null_Model_results_allometric_mean_dispersal.RData")

# Plot observed, simulated, effect size, and std. effect size relationship between median CNDD and species rarefied richness across plots
par(mfrow = c(2,2))
plot(null.results$obs.medianCNDD, plotdata$plotrarefy, xlim = c(min(null.results$obs.medianCNDD), max(null.results$obs.medianCNDD)), 
	xlab = "Observed median CNDD", ylab = "Plot rarefied species richness", las = 1)
abline(v = 0, lty = 2)
legend("bottomleft", legend = c(paste("r = ",round(cor.test(null.results$obs.medianCNDD, plotdata$plotrarefy,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(null.results$obs.medianCNDD, plotdata$plotrarefy, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(null.results$medianCNDD.sim, plotdata$plotrarefy, xlim = c(min(null.results$obs.medianCNDD), max(null.results$obs.medianCNDD)), 
	xlab = "Simulated median CNDD", ylab = "Plot rarefied species richness", las = 1)
abline(v = 0, lty = 2)
plot(null.results$medianCNDD.es, plotdata$plotrarefy, 
	xlab = "Median CNDD ES", ylab = "Plot rarefied species richness", las = 1)
abline(v = 0, lty = 2)
legend("bottomleft", legend = c(paste("r = ",round(cor.test(null.results$medianCNDD.es, plotdata$plotrarefy,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(null.results$medianCNDD.es, plotdata$plotrarefy, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(null.results$medianCNDD.ses, plotdata$plotrarefy, 
	xlab = "Median CNDD SES", ylab = "Plot rarefied species richness", las = 1)
abline(v = 0, lty = 2)
legend("bottomleft", legend = c(paste("r = ",round(cor.test(null.results$medianCNDD.ses, plotdata$plotrarefy,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(null.results$medianCNDD.ses, plotdata$plotrarefy, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
title("Null model results for median CNDD (Allometric dispersal across species)", outer = T, line = -2)



# Plot observed, simulated, effect size, and std. effect size for relationship between CNDD-species-abundance slopes and absolute latitude across plots
par(mfrow = c(2,2))
plot(abs(plotdata$lat),null.results$obs.cnddspp.slope, ylim = c(min(null.results$CNDDSPP.sim), max(null.results$obs.cnddspp.slope)), 
	xlab = "Absolute latitude", ylab = "Obs. slope btw. CNDD & spp. abund.", las = 1)
abline(h = 0, lty = 2)
legend("topright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$obs.cnddspp.slope,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$obs.cnddspp.slope, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(abs(plotdata$lat),null.results$CNDDSPP.sim, ylim = c(min(null.results$CNDDSPP.sim), max(null.results$obs.cnddspp.slope)), 
	xlab = "Absolute latitude", ylab = "Sim. slope btw. CNDD & spp. abund.", las = 1)
abline(h = 0, lty = 2)
plot(abs(plotdata$lat),null.results$CNDDSPP.es, 
	xlab = "Absolute latitude", ylab = "Slope btw. CNDD & spp. abund. ES", las = 1)
abline(h = 0, lty = 2)
legend("topright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$CNDDSPP.es,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$CNDDSPP.es, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(abs(plotdata$lat),null.results$CNDDSPP.ses, 
	xlab = "Absolute latitude", ylab = "Slope btw. CNDD & spp. abund. SES", las = 1)
abline(h = 0, lty = 2)
legend("topright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$CNDDSPP.ses,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$CNDDSPP.ses, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
title("Null model results for CNDD-spp. abund. (Allometric dispersal across spp.)", outer = T, line = -2)



# Plot observed, simulated, effect size, and std. effect size relationship between median CNDD of rare species and absolute latitude across plots
par(mfrow = c(2,2))
plot(abs(plotdata$lat),null.results$obs.rare.medianCNDD, ylim = c(min(null.results$obs.rare.medianCNDD), max(null.results$obs.rare.medianCNDD)), 
	xlab = "Observed median CNDD", ylab = "Observed median for rare spp. CNDD", las = 1)
abline(h = 0, lty = 2)
legend("bottomright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$obs.rare.medianCNDD,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$obs.rare.medianCNDD, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(abs(plotdata$lat),null.results$rare.medianCNDD.sim, ylim = c(min(null.results$obs.rare.medianCNDD), max(null.results$obs.cnddspp.slope)), 
	xlab = "Simulated median CNDD", ylab = "Simulated median for rare spp. CNDD", las = 1)
abline(h = 0, lty = 2)
plot(abs(plotdata$lat),null.results$rare.medianCNDD.es, 
	xlab = "Median CNDD ES", ylab = "Median CNDD for rare spp. ES", las = 1)
abline(h = 0, lty = 2)
legend("bottomright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$rare.medianCNDD.es,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$rare.medianCNDD.es, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
plot(abs(plotdata$lat),null.results$rare.medianCNDD.ses, 
	xlab = "Median CNDD SES", ylab = "Median CNDD for rare spp. SES", las = 1)
abline(h = 0, lty = 2)
legend("bottomright", legend = c(paste("r = ",round(cor.test(abs(plotdata$lat),null.results$rare.medianCNDD.ses,method="spearman")$estimate,3)), 
	paste("p = ",round(cor.test(abs(plotdata$lat),null.results$rare.medianCNDD.ses, method = "spearman")$p.value,3))), 
	bty = "n", cex = 0.9)
title("Null model results for median CNDD for rare spp. (Allometric dispersal across spp.)", outer = T, line = -2)














