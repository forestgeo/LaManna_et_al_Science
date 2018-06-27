### R script to simulate data with known value of CNDD, simulated error,
### and dispersal of seeds outside the forest plot for CTFS-ForestGEO analysis.
### 
### By: J. A. LaManna, updated 02-23-2017


# Load required packages
library(doBy)
library(reshape)
library(VGAM)
library(emdbook)
library(nlme)
library(ggplot2)
library(RColorBrewer)
library(mvtnorm)
library(doBy)
library(ggplot2)
library(pscl)
library(MASS)
library(boot)
library(spatstat)
library(vegan)
library(gnm)



################################################################################################
# I. Functions
################################################################################################
### Functions to simulate data with known CNDD
# n = number of quadrats in simulated forest plot
# meanTrees = number of mean adult trees per quadrat
# lambda = per capita recruitment rate (in absence of density dependence)
# trueCNDD = conspecific density dependence (0 = no CNDD)
# theta = negative binomial overdispersion parameter
# d = proportion of seeds dispersing outside of the simulated forest plot

################################################################################################
### Fit Functions

### Load functions for analyses

# Nonlinear Ricker function - Distance-weighted approach
Ricker <- function(x,Had,Hsap){
	list(predictors =	list(r = 1, CNDD = 1),
	variables = list(substitute(x)),
      term = function(predictors, variables) {
          pred <- paste("(", variables[1],")*exp(", predictors[1], 
				")*exp(", predictors[2], "*", variables[1],")+0.0001", sep = "")
         })
}
class(Ricker ) <- "nonlin"

fit.ricker.cndd.dist.wgt = function(data){
x = data$distwgtTrees; y = data$recruits
return(tryCatch(gnm(y~-1+Ricker(x),family = quasipoisson(link="identity")), error=function(e) NULL))
}


# Nonlinear Ricker function - Extra-intercept approach (proposed by Chisholm and Fung)
Ricker.Int <- function(x){
	list(predictors =	list(r = 1, CNDD = 1, Int = 1),
	variables = list(substitute(x)),
      term = function(predictors, variables) {
          pred <- paste(predictors[3], "+", "(", variables[1],")*exp(", predictors[1], 
				")*exp(", predictors[2], "*", variables[1],")+0.0001", sep = "")
         })
}
class(Ricker.Int ) <- "nonlin"

fit.ricker.int.cndd = function(data){
x = data$consppTrees; y = data$recruits
return(tryCatch(gnm(y~-1+Ricker.Int(x),family = quasipoisson(link="identity")), error=function(e) NULL))
}


### Other necessary functions

seed = function(r, alpha=100) {(1/(pi*alpha*((1+((r^2)/alpha))^2)))}	# Clark's 2Dt p = 1 

# Function to calculate distance-weighted adult abundances for each quadrat
dist.weighted.abund.function = function(test, quads) {
adultlocs = test[,c("gx", "gy")]
quadlocs = quads[,c("cx", "cy")]
xdist = crossdist(adultlocs[,1], adultlocs[,2], quadlocs[,1], quadlocs[,2])
xdist.weighted = seed(xdist, alpha = exp(6)) * (1/seed(0.1, alpha = exp(6)))
return(apply(xdist.weighted, 2, sum))
}


################################################################################################
### Data simulation functions

# Simple Model (no error)
sim.data.simple.ricker = function(h = 25, w = 50, meanTrees, r, trueCNDD, theta){
n = h * w
trueTrees = rnbinom(n, size = theta, mu = meanTrees)
cy1 = (1:h * 20) - 20
cx1 = (1:w * 20) - 20
quad.centers = expand.grid(cy1, cx1); quad.centers = data.frame(cx = quad.centers[,2], cy = quad.centers[,1])
data1 = data.frame(trueTrees, cx = quad.centers$cx, cy = quad.centers$cy)
data1$quadrat = c(1:nrow(data1))
trees.expanded = data1[rep(seq.int(1,nrow(data1)), data1$trueTrees),]
trees.expanded$gx = trees.expanded$cx + runif(nrow(trees.expanded),0.001,19.999)
trees.expanded$gy = trees.expanded$cy + runif(nrow(trees.expanded),0.001,19.999)
data1$dist.wgt.adults = dist.weighted.abund.function(trees.expanded, data1)
recruits = data1$dist.wgt.adults * exp(r + (trueCNDD * data1$dist.wgt.adults)) 	# Ricker population model
data=data.frame(consppTrees = trueTrees, distwgtTrees = data1$dist.wgt.adults, recruits)}


# Error Model
# same thing but add process error to observed recruits per quadrat, simulates demographic stochasticity
sim.data.error.ricker = function(h = 25, w = 50, meanTrees, r, trueCNDD, theta){
n = h * w
trueTrees = rnbinom(n, size = theta, mu = meanTrees)
cy1 = (1:h * 20) - 20
cx1 = (1:w * 20) - 20
quad.centers = expand.grid(cy1, cx1); quad.centers = data.frame(cx = quad.centers[,2], cy = quad.centers[,1])
data1 = data.frame(trueTrees, cx = quad.centers$cx, cy = quad.centers$cy)
data1$quadrat = c(1:nrow(data1))
trees.expanded = data1[rep(seq.int(1,nrow(data1)), data1$trueTrees),]
trees.expanded$gx = trees.expanded$cx + runif(nrow(trees.expanded),0.001,19.999)
trees.expanded$gy = trees.expanded$cy + runif(nrow(trees.expanded),0.001,19.999)
data1$dist.wgt.adults = dist.weighted.abund.function(trees.expanded, data1)
truerecruits = data1$dist.wgt.adults * exp(r + (trueCNDD * data1$dist.wgt.adults)) 	# Ricker population model
recruits = rnbinom(n, size = theta, mu = truerecruits)
data=data.frame(consppTrees = trueTrees, truerecruits, distwgtTrees = data1$dist.wgt.adults, recruits)}


# Dispersal and error model
# same as error version but some fraction (d) of recruits globally dispersed
sim.data.dispersal.ricker = function(h = 25, w = 50, meanTrees, r, trueCNDD, theta, d){
n = h * w
trueTrees = rnbinom(n, size = theta, mu = meanTrees)
cy1 = (1:h * 20) - 20
cx1 = (1:w * 20) - 20
quad.centers = expand.grid(cy1, cx1); quad.centers = data.frame(cx = quad.centers[,2], cy = quad.centers[,1])
data1 = data.frame(trueTrees, cx = quad.centers$cx, cy = quad.centers$cy)
data1$quadrat = c(1:nrow(data1))
trees.expanded = data1[rep(seq.int(1,nrow(data1)), data1$trueTrees),]
trees.expanded$gx = trees.expanded$cx + runif(nrow(trees.expanded),0.001,19.999)
trees.expanded$gy = trees.expanded$cy + runif(nrow(trees.expanded),0.001,19.999)
data1$dist.wgt.adults = dist.weighted.abund.function(trees.expanded, data1)
totRecruits = data1$dist.wgt.adults * exp(r + (trueCNDD * data1$dist.wgt.adults)) 	# Ricker population model
localRecruits = totRecruits * (1 - d)
recruits = localRecruits + ((sum(totRecruits) * d) / n)
recruits = rnbinom(n, size = theta, mu = recruits)
data=data.frame(consppTrees = trueTrees, truerecruits = totRecruits, distwgtTrees = data1$dist.wgt.adults, recruits)}




################################################################################################
# II. Simple: no observation or measurement error, no dispersal
################################################################################################
# Ricker models that do not use distance-weighted adult abundance are blind to number of 
# adults trees in neighboring quadrats

# Simple Ricker Simulation
h = 25
w = 50
meanTrees = 2
trueCNDD = -0.1	# True conspecific negative density-dependence (0 = no CNDD; lower values = stronger CNDD)
r = 0
theta = 1
data = sim.data.simple.ricker(h=h,w=w,meanTrees,r=r,trueCNDD=trueCNDD,theta=theta)
fit.ricker.cndd.dist.wgt(data) 
fit.ricker.int.cndd(data)



################################################################################################
# III. Add process error to recruits
################################################################################################
# add process error to number of recruits in each quadrat

meanTrees = 2
h = 25
w = 50
trueCNDD = -0.10	# True conspecific negative density-dependence (0 = no CNDD; lower values = stronger CNDD)
r = 0.00
theta = 10
data = sim.data.error.ricker(h = h, w = w, meanTrees = meanTrees, r = r, trueCNDD = trueCNDD, theta = theta)
fit.ricker.cndd.dist.wgt(data) 
fit.ricker.int.cndd(data) 



################################################################################################
# IV. Add dispersal and process error
################################################################################################
# Incorporate dispersal out of the plot, ((1-d)*100) recruits stay put, (d*100) are globally dispersed

meanTrees = 0.1	# Mean number of adult trees per quadrat
h = 25		# Number of quadrats the forest plot is in height
w = 50		# Number of quadrats the forest plot is in width
trueCNDD = -0.5	# True conspecific negative density-dependence (1 = no CNDD; lower values = stronger CNDD)
r = 0.00		# Density-independent population growth rate 
theta = 1		# Error
d = 0.10		# Dispersal Factor
data = sim.data.dispersal.ricker(h = h, w = w, meanTrees = meanTrees, r = r, trueCNDD = trueCNDD, theta = theta, d = d)
fit.ricker.cndd.dist.wgt(data) 
fit.ricker.int.cndd(data) 






########################################################################################################################
########################################################################################################################
### Evaluation of Ricker models over parameter space


### Set parameter values for simulations

h = 25			# Number of quadrats the forest plot is in height
w = 50			# Number of quadrats the forest plot is in width
its = 15 		# number of breaks between extreme values of CNDD in parameter space
k = 1			# number of iterations for each parameter combination

testCNDD = seq(-2.0,0.05,length=its)								# Range of CNDD values observed in data
testR = seq(-0.5,1.25,length=8)									# Based on actual mean recruitment (10% to 90%)
testmeanTrees = c(0.026, 0.043, 0.062, 0.092, 0.145, 0.214, 0.349, 0.625, 1.279)		# Based on deciles of mean adult abundance across all species in dataset (10% to 90%)
testd = seq(0.0,0.2,by=0.05)									# Increasing dispersal out of the plot
testTheta = c(1:10)										# Values of simulated process error (values of shape parameter in neg. binomial distribution)


### Create matrix of all parameter combinations for simulations
rickermatrix = data.frame(matrix(c(NA),nrow=its*length(testR)*length(testmeanTrees)*length(testd)*length(testTheta), ncol = 7))
names(rickermatrix) = c("knownCNDD", "knownR", "meanTrees", "d", "theta", "estCNDD", "estR")
rickermatrix$knownCNDD = rep(testCNDD, each = dim(rickermatrix)[1]/length(testCNDD))
rickermatrix$knownR = rep(testR, each = dim(rickermatrix)[1]/length(testCNDD)/length(testR))
rickermatrix$meanTrees = rep(testmeanTrees, each = dim(rickermatrix)[1]/length(testCNDD)/length(testR)/length(testmeanTrees))
rickermatrix$d = rep(testd, each = dim(rickermatrix)[1]/length(testCNDD)/length(testR)/length(testmeanTrees)/length(testd))
rickermatrix$theta = rep(testTheta, each = dim(rickermatrix)[1]/length(testCNDD)/length(testR)/length(testmeanTrees)/length(testd)/length(testTheta))


### Function to automatically run simulations for all parameter combinations - Distance-weighted appraoch
Ricker.simulation.function.dist.wgt.rows = function(test, k = 1) {
meanTrees = test['meanTrees']; r = test['knownR']; trueCNDD = test['knownCNDD']
theta = test['theta']; d = test['d']
data = replicate(k, sim.data.dispersal.ricker(h = h, w = w, meanTrees = meanTrees, r = r,	trueCNDD = trueCNDD, theta = theta, d = d), simplify = F)
Rickerfits = lapply(data, fit.ricker.cndd.dist.wgt)
fits2 = list()
for(j in 1:length(Rickerfits)) {
  possibleError <- tryCatch(summary(Rickerfits[[j]])$coef, error=function(e) e  )
  if(!inherits(possibleError, "error")) {best.pow.optim = c(t(summary(Rickerfits[[j]])$coef[,1:2]))
names(best.pow.optim) = paste(rep(rownames(summary(Rickerfits[[j]])$coef[,1:2]), each = 2), rep(c("", ".se"), times = 1), sep = "")
fits2[[j]] = best.pow.optim
}}
fits = data.frame(do.call("rbind", fits2))
fits = fits[!is.na(fits$CNDD.se),]
if(max(fits$CNDD.se) > 100) {fits = fits[-which(fits$CNDD.se > 100),]}
fit.summary = c(colMeans(fits), trueCNDD, r, meanTrees, d, theta)
names(fit.summary)[5:9] = c("knownCNDD", "knownR", "meanTrees", "d", "theta")
return(fit.summary)
}

### Function to automatically run simulations for all parameter combinations - Extra-intercept approach (proposed by Chisholm and Fung)
Ricker.simulation.function.int.rows = function(test, k = 10) {
meanTrees = test['meanTrees']; r = test['knownR']; trueCNDD = test['knownCNDD']
theta = test['theta']; d = test['d']
data = replicate(k, sim.data.dispersal.ricker(h = h, w = w, meanTrees = meanTrees, r = r,	trueCNDD = trueCNDD, theta = theta, d = d), simplify = F)
Rickerfits = lapply(data, fit.ricker.int.cndd)
fits2 = list()
for(j in 1:length(Rickerfits)) {
  possibleError <- tryCatch(summary(Rickerfits[[j]])$coef, error=function(e) e  )
  if(!inherits(possibleError, "error")) {best.pow.optim = c(t(summary(Rickerfits[[j]])$coef[,1:2]))
names(best.pow.optim) = paste(rep(rownames(summary(Rickerfits[[j]])$coef[,1:2]), each = 2), rep(c("", ".se"), times = 1), sep = "")
fits2[[j]] = best.pow.optim
}}
fits = data.frame(do.call("rbind", fits2))
fits = fits[!is.na(fits$CNDD.se),]
if(max(fits$CNDD.se) > 100) {fits = fits[-which(fits$CNDD.se > 100),]}
if(nrow(fits) > 0) {
fit.summary = c(colMeans(fits), trueCNDD, r, meanTrees, d, theta)
names(fit.summary)[7:11] = c("knownCNDD", "knownR", "meanTrees", "d", "theta")
return(fit.summary)
} else {return(rep(NA,times=11))}
}




### Run this code for benchmark tests of distance-weighted approach
set.seed(254)
begin.time = Sys.time()
rickersims = apply(rickermatrix, 1, Ricker.simulation.function.dist.wgt.rows, k = 1)
rickersims2 = data.frame(t(rickersims))
end.time = Sys.time()
(duration.ricker = end.time - begin.time)
rickermatrix.dst.wgt = rickersims2
names(rickermatrix.dst.wgt) = c("r", "r.se", "CNDD", "CNDD.se", "knownCNDD", "knownR", "meanTrees", "d", "theta")


### Run this code for benchmark tests of extra-intercept approach (proposed by Chisholm and Fung)
set.seed(254)
begin.time = Sys.time()
rickersims = apply(rickermatrix, 1, Ricker.simulation.function.int.rows, k = 1)
rickersims2 = data.frame(t(rickersims))
end.time = Sys.time()
(duration.ricker = end.time - begin.time)
rickermatrix.ext.int = rickersims2[!is.na(rickersims2[,1]),]
names(rickermatrix.ext.int) = c("r", "r.se", "CNDD", "CNDD.se", "Int", "Int.se", "knownCNDD", "knownR", "meanTrees", "d", "theta")






##########################################################################################
##########################################################################################
### FIGURES 


################################################################
### Known CNDD vs. estimated CNDD - Distance-weighted appraoch
### Expected slope of 1
### This is Fig. 2b in our response to Chisholm & Fung, LaManna et al. 2018.
### The distance-weighted approach reasonably recovers known values of CNDD after process error and 
### immigration out of the plot are simulated (slope and correlation near 1)

rickermatrix = rickermatrix.dst.wgt
rickermatrix2 = rickermatrix[order(rickermatrix$knownCNDD),]
x = rickermatrix2$knownCNDD
y = rickermatrix2$CNDD
plot(x,y,xlab="Known CNDD value",ylab="Estimated CNDD value",las=1, ylim = c(-5.75,3.75), cex=0.5,
	main = "Distance-weighted Adult Abundance", type = "n")
boxplot(y~x,las=3,at=(unique(x)),add=T,xaxt="n",pch=19,cex=0.3,
	pars = list(boxwex = 0.09, staplewex = 0.8, outwex = 0.8),
	col="red",yaxt="n")
abline(0,1,lty=2)
legend("bottomright",legend=c(paste("r = ",round(cor.test(x,y)$estimate,3)),paste("slope = ",round(summary(lm(y~x))$coef[2,1],3))),bty="n",cex=0.8)




################################################################
### Known CNDD vs. estimated CNDD - Extra-intercept approach (proposed by Chisholm and Fung) 
### Expected slope of 1
### This is Fig. 2a in our response to Chisholm & Fung, LaManna et al. 2018.
### The extra-intercept approach cannot accurately recover known values of CNDD after process error and 
### immigration out of the plot are simulated (slope and correlation not close to 1). 
### This approach estimates weak values of CNDD - near zero - even when the known value of 
### CNDD is strong (i.e. the left side of the figure).

rickermatrix = rickermatrix.ext.int
rickermatrix2 = rickermatrix[order(rickermatrix$knownCNDD),]
x = rickermatrix2$knownCNDD
y = rickermatrix2$CNDD
plot(x,y,xlab="Known CNDD value",ylab="Estimated CNDD value",las=1, ylim = c(-5.75,3.75), cex=0.5,
	main = "Distance-weighted Adult Abundance", type = "n")
boxplot(y~x,las=3,at=(unique(x)),add=T,xaxt="n",pch=19,cex=0.3,
	pars = list(boxwex = 0.09, staplewex = 0.8, outwex = 0.8),
	col="red",yaxt="n")
abline(0,1,lty=2)
legend("bottomright",legend=c(paste("r = ",round(cor.test(x,y)$estimate,3)),paste("slope = ",round(summary(lm(y~x))$coef[2,1],3))),bty="n",cex=0.8)




################################################################
### Known CNDD vs. estimated CNDD for rarest 1/3 of species - Distance-weighted appraoch 
### Expected slope of 1

rickermatrix = rickermatrix.dst.wgt
meanTreestestlevel = 0.07
rickermatrix2 = rickermatrix[order(rickermatrix$knownCNDD),]
rickermatrix2 = rickermatrix2[which(rickermatrix2$meanTrees <= meanTreestestlevel),]
x = rickermatrix2$knownCNDD
y = rickermatrix2$CNDD
plot(x,y,xlab="Known CNDD value",ylab="Estimated CNDD value",las=1, ylim = c(-5.75,3.75), cex=0.5,
	main = "Distance-weighted Adult Abundance", type = "n")
boxplot(y~x,las=3,at=(unique(x)),add=T,xaxt="n",pch=19,cex=0.3,
	pars = list(boxwex = 0.09, staplewex = 0.8, outwex = 0.8),
	col="red",yaxt="n")
abline(0,1,lty=2)
legend("bottomright",legend=c(paste("r = ",round(cor.test(x,y)$estimate,3)),paste("slope = ",round(summary(lm(y~x))$coef[2,1],3))),bty="n",cex=0.8)




################################################################
### Known CNDD vs. estimated CNDD for rarest 1/3 of species - Extra-intercept approach (proposed by Chisholm and Fung)
### Expected slope of 1
### The extra-intercept approach cannot accurately recover known values of CNDD, and the 
### extreme positive bias is greatest when CNDD is stronger (i.e. more negative, on the left side of the figure).

rickermatrix = rickermatrix.ext.int
meanTreestestlevel = 0.07
rickermatrix2 = rickermatrix[order(rickermatrix$knownCNDD),]
rickermatrix2 = rickermatrix2[which(rickermatrix2$meanTrees <= meanTreestestlevel),]
x = rickermatrix2$knownCNDD
y = rickermatrix2$CNDD
plot(x,y,xlab="Known CNDD value",ylab="Estimated CNDD value",las=1, ylim = c(-5.75,3.75), cex=0.5,
	main = "Distance-weighted Adult Abundance", type = "n")
boxplot(y~x,las=3,at=(unique(x)),add=T,xaxt="n",pch=19,cex=0.3,
	pars = list(boxwex = 0.09, staplewex = 0.8, outwex = 0.8),
	col="red",yaxt="n")
abline(0,1,lty=2)
legend("bottomright",legend=c(paste("r = ",round(cor.test(x,y)$estimate,3)),paste("slope = ",round(summary(lm(y~x))$coef[2,1],3))),bty="n",cex=0.8)




################################################################
### CNDD bais (estimated CNDD minus known CNDD) across species abundances - Distance-weighted appraoch 
### Expected slope of 0
### Minimal bias for rare species, but only a fraction of the observed values of CNDD for 
### rare species (see Fig. 4 in our response to Chisholm & Fung, LaManna et al. 2018). 
### Also, this minimal bias is corrected by the use of the null model presented in our
### response to Huelsmann & Hartig, LaManna et al. Science 2018.

rickermatrix = rickermatrix.dst.wgt
rickermatrix2 = rickermatrix[order(rickermatrix$meanTrees),]
x = log(rickermatrix2$meanTrees)
y = (rickermatrix2$CNDD - rickermatrix2$knownCNDD)
plot(x,y,xlab="Mean adult abundance",ylab="CNDD bias",las=1, ylim = c(-8, 7), cex=0.5, xaxt = "n", 
	main = "Distance-weighted Adult Abundance", type = "n")
axis(1,at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),labels=c(0.0001,0.001,0.01,0.1,1,10,100))
ticks = c(seq(0.0001,0.001,length=10),seq(0.002,0.01,length=9),seq(0.02,0.1,length=9),seq(0.2,1,length=9),seq(2,10,length=9),seq(20,100,length=9))
axis(1, at = c(log(ticks)), labels = c(rep("", times = length(ticks))))
abline(h=0,lty=2,col="gray55")
boxplot(y~x,las=3,at=(unique(x)),add=T,xaxt="n",pch=19,cex=0.3,
	pars = list(boxwex = 0.18, staplewex = 0.8, outwex = 0.8),
	col="red",yaxt="n")
legend("bottomright",legend=c(paste("slope = ",round(summary(lm(y~x))$coef[2,1],3))),bty="n",cex=0.8)




################################################################
### CNDD bais (estimated CNDD minus known CNDD) across species abundances - Extra-intercept approach (proposed by Chisholm and Fung) 
### Expected slope of 0
### Extreme positive bias in CNDD for rare species, part of the reason their analysis removes the 
### latitudinal pattern in CNDD

rickermatrix = rickermatrix.ext.int
rickermatrix2 = rickermatrix[order(rickermatrix$meanTrees),]
x = log(rickermatrix2$meanTrees)
y = (rickermatrix2$CNDD - rickermatrix2$knownCNDD)
plot(x,y,xlab="Mean adult abundance",ylab="CNDD bias",las=1, ylim = c(-8, 7), cex=0.5, xaxt = "n", 
	main = "Distance-weighted Adult Abundance", type = "n")
axis(1,at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),labels=c(0.0001,0.001,0.01,0.1,1,10,100))
ticks = c(seq(0.0001,0.001,length=10),seq(0.002,0.01,length=9),seq(0.02,0.1,length=9),seq(0.2,1,length=9),seq(2,10,length=9),seq(20,100,length=9))
axis(1, at = c(log(ticks)), labels = c(rep("", times = length(ticks))))
abline(h=0,lty=2,col="gray55")
boxplot(y~x,las=3,at=(unique(x)),add=T,xaxt="n",pch=19,cex=0.3,
	pars = list(boxwex = 0.18, staplewex = 0.8, outwex = 0.8),
	col="red",yaxt="n")
legend("topright",legend=c(paste("slope = ",round(summary(lm(y~x))$coef[2,1],3))),bty="n",cex=0.8)




##########################################################################################
### Evaluation of the immigration parameter for the Extra-intercept approach (proposed by Chisholm and Fung) 
### Expected slope of 1, or correlation coefficient near 1
### However, no correlation between known immigration values and those 
### estimated by the extra-intercept appraoch

rickermatrix = rickermatrix.ext.int
rickermatrix2 = rickermatrix[order(rickermatrix$knownCNDD),]
x = rickermatrix2$d
y = rickermatrix2$Int
plot(x,y,xlab="Known immigration value",ylab="Estimated immigration value",las=1, 
	ylim = c(0, 5), cex=0.5, main = "Ricker model with extra intercept", type = "n")
boxplot(y~x,las=3,at=(unique(x)),add=T,xaxt="n",pch=19,cex=0.3,
	pars = list(boxwex = 0.01, staplewex = 0.8, outwex = 0.8),
	col="red",yaxt="n")
abline(0,1,lty=2)
legend("topleft",legend=c(paste("r = ",round(cor.test(x,y)$estimate,3)),paste("slope = ",round(summary(lm(y~x))$coef[2,1],3))),
	bty="n",cex=0.8,inset=0.03)


### Full extent of estimated intercept values against known immigration values
rickermatrix = rickermatrix.ext.int
rickermatrix2 = rickermatrix[order(rickermatrix$knownCNDD),]
x = rickermatrix2$d
y = rickermatrix2$Int
plot(x,y,xlab="Known immigration value",ylab="Estimated immigration value",las=1, 
	cex=0.5, main = "Ricker model with extra intercept", type = "n")
boxplot(y~x,las=3,at=(unique(x)),add=T,xaxt="n",pch=19,cex=0.3,
	pars = list(boxwex = 0.01, staplewex = 0.8, outwex = 0.8),
	col="red",yaxt="n")
abline(0,1,lty=2)
legend("topleft",legend=c(paste("r = ",round(cor.test(x,y)$estimate,3)),paste("slope = ",round(summary(lm(y~x))$coef[2,1],3))),
	bty="n",cex=0.8,inset=0.03)
















