##- This file contains an artificial example to illustrate the functionality of the 
##- functions in file mycrr.r
##- 
##- Only limited testing of the functions has been performed and they come without any warranty 
##- 
##- If you use this code please acknowledge it and cite our paper:
##- Wolbers M, Koller MT, Witteman JCM, Steyerberg EW: 
##- Prognostic Models With Competing Risks - Methods and Application to Coronary Risk Prediction
##- Epidemiology (2009): in press
##-
##- For feedback and comments please contact Marcel Wolbers, mwolbers@oucru.org

# load functions
source("mycrr.r")

# create test data
set.seed(1)
n <- 100
x1 <- rnorm(n)
ftime <- rexp(n,rate=exp(x1))
fstatus <- sample(c(0,1,2),size=n,replace=T,prob=c(0.4,0.3,0.3)) 
dataset <- data.frame(x1=x1,x2=runif(n,min=0.1,max=10),x3=runif(n,min=0.1,max=10),x4=rnorm(n),z=rbinom(n,size=1,prob=0.5))

##--------------
##- estimate a model w/o covariate-time interactions
# with mycrr
test.mycrr <- mycrr(formula=crSurv(ftime,fstatus)~x1+x2+I(x2^2)+log(x3),data=dataset,failcode=1,cencode=0,cengroup=z,print.output=T)
par(mfrow=c(2,2)); plot.mycrr(test.mycrr) #plot Schoenfeld residuals
# same with crr
cov1 <- as.matrix(cbind(dataset[,1:2],dataset[,2]^2,log(dataset[,3])))
( test.crr <- crr(ftime=ftime,fstatus=fstatus,cov1=cov1,cengroup=dataset$z,failcode=1,cencode=0) )

##--------------
##- estimate a model with covariate-time interactions
# with mycrr
test.mycrr.td <- mycrr(formula=crSurv(ftime,fstatus)~x1+x2+I(x2^2)+log(x3),
                       covtime.formula= ~(x1)*(log(t))+(x2)*(log(t))+(log(x2))*(t) +(log(x2))*(t^2),
                       data=dataset,failcode=1,cencode=0,cengroup=z,print.output=T)
# same with crr
cov2 <- as.matrix(cbind(dataset[,1],dataset[,2],log(dataset[,2]),log(dataset[,2])))
( test.crr.td <- crr(ftime=ftime,fstatus=fstatus,cov1=cov1,cov2=cov2,tf=function(t) cbind(log(t),log(t),t,t^2),
                     cengroup=dataset$z,failcode=1,cencode=0) )

##--------------
##- predict cumulative incidence from model w/o time-dependent covariates given new data
par(mfrow=c(2,2))
newData <- data.frame(x1=rnorm(4),x2=runif(4,min=0.1,max=10),x3=runif(4,min=0.1,max=10),x4=rnorm(4),z=rbinom(4,size=1,prob=0.5))
# based on mycrr
test.predict.mycrr <- predict.mycrr(test.mycrr,newData) 
plot(test.predict.mycrr, main="mycrr\n (no cov-time interactions)",col=1:4)
# based on crr
newCov1 <- as.matrix(cbind(newData[,1:2],newData[,2]^2,log(newData[,3])))
test.predict.crr <- predict(test.crr,newCov1) 
plot(test.predict.crr, main="crr\n (no cov-time interactions)",col=1:4)
# compare the two outputs
all.equal(test.predict.mycrr,test.predict.crr)

# linear predictor for new data
predict.mycrr(test.mycrr,newData,type="linear.predictor")
# predictions at specific timepoints
predict.mycrr(test.mycrr,newData,type="subdist",times=1:5)

##--------------
##- predict cumulative incidence from model with covariate-time interactions given new data
# based on mycrr.td
test.predict.mycrr.td <- predict.mycrr(test.mycrr.td,newData) 
plot(test.predict.mycrr.td, main="mycrr\n (with cov-time interactions)",col=1:4)
# based on crr
newCov2 <- as.matrix(cbind(newData[,1],newData[,2],log(newData[,2]),log(newData[,2])))
test.predict.crr.td <- predict(test.crr.td,newCov1,newCov2) 
plot(test.predict.crr.td, main="crr\n (with cov-time interactions)",col=1:4)
# compare the two outputs
all.equal(test.predict.mycrr.td,test.predict.crr.td)

##--------------
##- C index and Royston&Sauerbrei D for test.mycrr
lp.test <- predict.mycrr(test.mycrr,dataset,type="linear.predictor") 
# C index (both lines below give the same)
cIndex.CR(predictor=lp.test,ftime=ftime,fstatus=fstatus,failcode=1,cencode=0,conf.int=T,nBoot=100)
rcorr.mycrr(linear.predictor=lp.test,ftime=ftime,fstatus=fstatus,failcode=1,cencode=0)
# Royston&Sauerbrei D
D.mycrr(linear.predictor=lp.test,ftime=ftime,fstatus=fstatus,failcode=1,cencode=0)

##--------------
##- C index for model test.mycrr.td (with time-dependent covariates)
tmp <- predict.mycrr(test.mycrr.td,covariateData=dataset) 
unique.failtimes <- tmp[,1]
cIndex.CR(predictor=t(tmp[,-1]),ftime=ftime,fstatus=fstatus,failcode=1,cencode=0,conf.int=T,nBoot=100)