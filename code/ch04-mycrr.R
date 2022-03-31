##- This is a list of functions that extend the functionality of the cmprsk library
##- See the associated example code "test.mycrr.r" for an illustration of how to use them
##- 
##- Only limited testing of the functions has been performed and they come without any warranty 
##- 
##- If you use this code, please acknowledge it and cite our paper:
##- Wolbers M, Koller MT, Witteman JCM, Steyerberg EW: 
##- Prognostic Models With Competing Risks - Methods and Application to Coronary Risk Prediction
##- Epidemiology (2009): in press
##-
##- For feedback and comments please contact Marcel Wolbers, mwolbers@oucru.org


crSurv <- function(ftime,fstatus){
  ##- simple competing risk Surv-type object
  cbind(ftime,fstatus)
}

split.covtime.formula <- function(covtime.formula){
  ##- split up covtime.formula argument of function mycrr
  ##- see comments in function mycrr below for the structure of the argument covtime.formula
  ##- example call: split.covtime.formula(~(x1)*(log(t))+(log(x2))*(t)+(log(x2))*(t^2))
  ##-
  ##- Author: Marcel Wolbers 
  ##- Last update: Jan 29, 2008 
  ##--------------------------------------------
  covtime.formula <- as.character(covtime.formula)[2]
  covtime.formula.chars <- unlist(strsplit(covtime.formula,split=""))
  # find brackets
  openingbrackets <- !is.na(match(covtime.formula.chars,"("))
  closingbrackets <- !is.na(match(covtime.formula.chars,")"))
  # count only brackets opening or closing an "outer" expression, 
  # i.e. enclosing a covariate or a function of time
  rel.openingbrackets <- openingbrackets
  rel.openingbrackets[cumsum(openingbrackets)>cumsum(closingbrackets)+1] <- 0
  rel.closingbrackets <-  closingbrackets
  rel.closingbrackets[cumsum(openingbrackets)!=cumsum(closingbrackets)] <- 0
  # covariates and functions of t
  if (sum(rel.openingbrackets)%%2!=0) stop("Error in processing argument covtime.formula!")
  numTerms <- sum(rel.openingbrackets)/2
  td.covariates <- tfunc <- inbetween.pluses <- inbetween.times <- vector("character",length=numTerms)
  for (i in 1:numTerms){
    td.covariates[i] <- substr(covtime.formula,start=which(cumsum(rel.openingbrackets)==2*i-1)[1]+1,
                                               stop=which(cumsum(rel.closingbrackets)==2*i-1)[1]-1) 
    tfunc[i] <- substr(covtime.formula,start=which(cumsum(rel.openingbrackets)==2*i)[1]+1,
                                       stop=which(cumsum(rel.closingbrackets)==2*i)[1]-1)
    # for minimal error checking;  inbetween.times should only contain "*" 
    inbetween.times[i] <- substr(covtime.formula,start=which(cumsum(rel.closingbrackets)==2*i-1)[1]+2,
                                                 stop=which(cumsum(rel.openingbrackets)==2*i)[1]-2)
    # for minimal error checking, inbetween.pluses should only contain "+" 
    inbetween.pluses[i] <- substr(covtime.formula,start=which(cumsum(rel.closingbrackets)==2*i)[1]+2,
                                                  stop=which(cumsum(rel.openingbrackets)==2*i+1)[1]-2)
  }
  inbetween.pluses[numTerms] <- "+"
  if (any(inbetween.times!="*")|any(inbetween.pluses!="+")) stop("Error in processing argument covtime.formula!")
  list(td.covariates=td.covariates,tfunc=tfunc)
} 

mycrr <- function(formula,covtime.formula=NULL,data,contrasts=NULL,subset,na.action=na.omit,
                  cengroup,failcode,cencode=0, gtol=1e-06, maxiter=10,print.output=T,...){
  ##- This is a wrapper for function crr in library(cmprsk) which has a formula interface and produces a nicer output 
  ##-
  ##- Def. of some of the arguments:
  ##- - formula: To specify time-independent covariate effects
  ##-            format: e.g. crSurv(ftime,fstatus)~x1+log(x2) 
  ##-            --> used to create argument cov1 in the call to crr
  ##- 
  ##- - covtime.formula: To specificy covariate-time interactions 
  ##             (Note: This argument does not correspond to a real R-formula but to a not very sophisticated custom-made "formula language"
  ##                    for covariate-time interactions)
  ##-            format: e.g. ~(x1)*(log(t))+(log(x2))*(t)+(log(x2))*(t^2) 
  ##-                    [interaction between x1 and log(time) as well as interaction between log(x2) and a quadratic function of time] 
  ##-            formatting conventions: 
  ##-                    - Interactions between covariate and functions of time must be of the form (cov)*(f(t)); 
  ##-                      the order (i.e. covariate first) and brackets surrounding cov and f(t) are mandatory!
  ##-                    - Different covariate-time interactions must be separated by "+"  
  ##-                    - The function of time must always be functions of "t" (i.e. "t" is the generic variable name to denote time)
  ##-            --> used to create arguments cov2 and tf in the call to crr
  ##- 
  ##-  - cengroup, failcode, cencode, maxiter: As for function crr in library(cmprsk) 
  ##
  ##- Author: Marcel Wolbers 
  ##- Last update: Dec 4, 2008 (minimal change: output gtol, maxiter, and cencode also in results) 
  ##- (limited testing only, use at your own risk!)
  ##--------------------------------------------
  require(cmprsk)
  ##- match.call etc. as in lm
  mf.base <- match.call()
  m <- match(c("formula", "data","subset","na.action","cengroup"), names(mf.base), 0)
  mf.base <- mf.base[c(1, m)]
  mf.base$drop.unused.levels <- TRUE
  mf.base[[1]] <- as.name("model.frame")
  mf <- eval(mf.base, parent.frame())
  mt <- attr(mf, "terms")
  ## -------------------
  ##- set up data for crr: response and time-independent covariate effects (argument cov1 to crr)
  ftime <- model.response(mf)[,1]
  fstatus <- model.response(mf)[,2]
  if (!missing(cengroup)) cengroup <- model.extract(mf,"cengroup")
  if (missing(cengroup)) cengroup <- 1
  cov1 <- model.matrix(mt,mf,contrasts)
  if ((dimnames(cov1)[[2]])[1]=="(Intercept)") cov1 <- cov1[,-1,drop=F] # remove intercept
  # programming note: A more 'elegant' way to get rid of the intercept would be to set "attributes(mt)$intercept <- 0" upfront.
  #                   However, this leads to problems for factors which are then coded differently, i.e. a binary factor gets 
  #                   two indicator variables in this case, not one!!
  terms.cov1 <- mt
  ## -------------------
  ##- set up arguments cov2 and tf for crr (i.e. arguments related to covariate-time interactions)
  if (is.null(covtime.formula)) {
    cov2 <- NULL
    terms.cov2 <- NULL
  }
  if (!is.null(covtime.formula)){
    #- split covtime.formula into pieces
    splitted.formula <- split.covtime.formula(covtime.formula)
    td.covariates <- splitted.formula$td.covariates
    tfunc <- splitted.formula$tfunc
    numTerms <- length(splitted.formula$td.covariates)
    #- create cov2 and tf from td.covariates and tfunc
    cov2 <- NULL
    terms.cov2 <- NULL
    tf <- "function(t) cbind("
    for (i in 1:numTerms){
       #- create cov2 part corresponding to the respective term
       mf.base[[2]] <- as.formula(paste("~",td.covariates[i]))
       mf <- eval(mf.base, parent.frame())
       mt <- attr(mf, "terms")
       terms.cov2[[i]] <- mt
       cov2.i <- model.matrix(mt,mf,contrasts)
       if ((dimnames(cov2.i)[[2]])[1]=="(Intercept)") cov2.i <- cov2.i[,-1,drop=F] #remove intercept
       colnames(cov2.i) <- paste("(",colnames(cov2.i),")*(",tfunc[i],")",sep="")
       cov2 <- cbind(cov2,cov2.i)
       #- create tf
       for (j in 1:ncol(cov2.i)){
         if ((i==numTerms)&(j==ncol(cov2.i))) tf <- paste(tf,tfunc[i],")",sep="")  
         else tf <- paste(tf,tfunc[i],",",sep="")
       }
    }
  }
  ## -------------------    
  ##- call crr 
  if (is.null(covtime.formula)){
    call.txt <- paste("crr(ftime=ftime,fstatus=fstatus,cov1=cov1",
                    ",cengroup=cengroup,failcode=failcode,cencode=cencode,gtol=gtol,maxiter=maxiter)",sep="")
  }
  if (!is.null(covtime.formula)){
    call.txt <- paste("crr(ftime=ftime,fstatus=fstatus,cov1=cov1,cov2=cov2,tf=",tf,
                    ",cengroup=cengroup,failcode=failcode,cencode=cencode,gtol=gtol,maxiter=maxiter)",sep="")
  }
  crr.result <- eval(parse(text=call.txt))
  ##- create nice output
  coef.crr <- crr.result$coef
  se.crr <-  sqrt(diag(crr.result$var))
  pVal.crr <- 2 *(1-pnorm(abs(coef.crr)/se.crr))
  result <- data.frame(coef=coef.crr,exp.coef=exp(coef.crr),se=se.crr,p=pVal.crr,
                       lower.95=exp(coef.crr-qnorm(0.975)*se.crr),
                       upper.95=exp(coef.crr+qnorm(0.975)*se.crr))
  if (is.null(covtime.formula)) coef.names <- colnames(cov1)
  if (!is.null(covtime.formula))  coef.names <- c(colnames(cov1),colnames(cov2))
  row.names(result) <- names(crr.result$coef) <- 
    rownames(crr.result$var) <- colnames(crr.result$var) <- coef.names
  #- print output
  if (print.output){
    cat("Convergence: ", crr.result$converged, "\n\n")
    cat("n=",nrow(cov1)," records used for the analysis.\n\n",sep="")
    print(round(result[,1:4],5)); cat("\n")
    print(round(result[,c(2,5,6)],5))
    #output Wald-type overall test for covariate*time interactions
    if (!is.null(covtime.formula)){
    if (ncol(cov2)>1){ 
      cat("\n")
      p1 <- ncol(cov1); p2 <- ncol(cov2)
      waldstat <-  t(as.matrix(coef.crr[(p1+1):(p1+p2)]))%*%
                   solve(crr.result$var[(p1+1):(p1+p2),(p1+1):(p1+p2)])%*%
                   as.matrix(coef.crr[(p1+1):(p1+p2)])
      cat("\nOverall Wald-type test for 'H0: all covariate*time-interactions=0': p=",1-pchisq(waldstat,df=p2),"\n")
    }
    }
  }
  invisible(list(result=result,crr.obj=crr.result,
                 formula=formula,covtime.formula=covtime.formula,
                 terms.cov1=terms.cov1,terms.cov2=terms.cov2,
                 ftime=ftime,fstatus=fstatus,cov1=cov1,cov2=cov2,
                 failcode=failcode,cencode=cencode,cengroup=cengroup,
                 gtol=gtol,maxiter=maxiter))
}

plot.mycrr <- function(mycrr.obj){
  ##- Diagnostic plot of Schoenfeld residuals against unique failure times
  res <- (mycrr.obj$crr.obj)$res
  unique.failuretimes <- sort(unique(mycrr.obj$ftime[mycrr.obj$fstatus==mycrr.obj$failcode]))
  par(ask=T)
  for (i in 1:ncol(res)){
    plot(unique.failuretimes,res[,i],main=row.names(mycrr.obj$result)[i])
    lines(lowess(unique.failuretimes,res[,i]),col=2)
    #p.val <- (summary(lm(res[,i]~unique.failuretimes))$coef)[2,4]
    #legend(x="topright",legend=paste("pVal (slope=0):",round(p.val,4)))
  }
  par(ask=F)
  NULL
}

predict.mycrr <- function(mycrr.object,covariateData,type="subdist",times=NULL){
  ##- Prediction for mycrr object (based on predict.crr)
  ##- In addition to predict.crr, it also gives predictions at specific timepoints (if times are provided)
  ##-  or the linear predictor (if type="linear.predictor")
  ##- !! Plotting (plot.crr) works only with the default arguments, i.e. type="subdist", times=NULL 
  ##-
  ##- Author: Marcel Wolbers
  ##- Last update: Jan 28, 2008 to work also with covariate-time interactions 
  ##--------------------------------------------
  # Creata cov1 based on data in covariatedata
  terms1 <- delete.response(mycrr.object$terms.cov1)
  cov1 <- model.matrix(terms1,model.frame(terms1,covariateData)) #predict.lm is more complicated also dealing properly with contrasts etc...
  if ((dimnames(cov1)[[2]])[1]=="(Intercept)") cov1 <- cov1[,-1,drop=F] #remove intercept
  # Create cov2 based on data in covariatedata 
  if (is.null(mycrr.object$covtime.formula)) cov2 <- NULL
  if (!is.null(mycrr.object$covtime.formula)){
    splitted.formula <- split.covtime.formula(mycrr.object$covtime.formula)
    td.covariates <- splitted.formula$td.covariates
    numTerms <- length(splitted.formula$td.covariates)
    cov2 <- NULL
    for (i in 1:numTerms){
      terms2.i <- mycrr.object$terms.cov2[[i]]
      cov2.i <- model.matrix(terms2.i,model.frame(terms2.i,covariateData))
      if ((dimnames(cov2.i)[[2]])[1]=="(Intercept)") cov2.i <- cov2.i[,-1,drop=F] #remove intercept
      cov2 <- cbind(cov2,cov2.i)
    }
  }
  # Predictions
  if (type=="subdist"){
    result <- predict.crr(object=mycrr.object$crr.obj,cov1=cov1,cov2=cov2)
    if (length(times)>0){
      tmp.result <- matrix(NA,ncol=length(times),nrow=nrow(cov1))
      colnames(tmp.result) <- paste("t=",times,sep="")
      for (i in 1:length(times)){
        if (times[i]<=max(result[,1])) tmp.result[,i] <- result[max(which(result[,1]<=times[i])),-1] 
      }
      result <- tmp.result
    }
  }
  if ((type=="linear.predictor")&(!is.null(cov2))) 
    stop("Linear predictor not sensible in models with covariate-time interactions!")
  if ((type=="linear.predictor")&(is.null(cov2)))  
    result <- c(cov1%*%(mycrr.object$result)$coef)
  result
} 

cIndex.CR <- function(predictor,ftime,fstatus,failcode=1,cencode=0,conf.int=F,nBoot=100){
  ##- Adapted c index for competing risk data
  ##- Similar to function rcorr.mycrr but more general as it also works for Fine&Gray models with time-dependent covariates
  ##- and cause-specific hazards competing risk models. In addition, there's an option to get bootstrap confidence intervals.
  ##-
  ##- Def. of some of the arguments:
  ##- -  predictor:
  ##-    let n be the number of observations (e.g. length(ftime) ) and
  ##-    m the number of unique failure times (i.e. length(unique(ftime[fstatus==failcode])) )
  ##-    Predictor is
  ##-    - either a numeric vector of length n containing a quantity which is monotonically related to the predicted incidence at any time
  ##-      (usually the linear predictor from a Fine&Gray model without covariate-time interactions)
  ##-    - or a matrix with n rows and m columns containing the predicted risk for each observation at each unique failure time
  ##-      (i.e. entry i,j is the predicted risk of observation i at the j-th ORDERED unique failure time)
  ##- - conf.int: if T, a bootstrap confidence interval is calculated
  ##- - nBoot: number of bootstrap samples (default 100 is rather too low, but high nBoot means longer computation time)
  ##-
  ##- Author: Marcel Wolbers, 
  ##- Last update: Feb 5, 2008 (limited testing only, use at your own risk!)
  ##--------------------------------------------
  n <- length(ftime)
  n.uniquefail <- length(unique(ftime[fstatus==failcode]))
  if (is.vector(predictor)) { if (length(predictor)!=n) stop("Argument predictor has wrong length, should be ",n) }
  if (is.matrix(predictor)) { if ((nrow(predictor)!=n)|(ncol(predictor)!=n.uniquefail))
    stop("Argument predictor has wrong dimensions, should be ",n,"*",n.uniquefail) }
  # patients with competing events censored at infinity (or, equivalently, max(ftime)+1)
  ftime2 <- ftime
  ftime2[(fstatus!=failcode)&(fstatus!=cencode)] <- max(ftime)+1
  # sort data by ftime2
  ftime2.order <- order(ftime2,partial=ifelse(fstatus==cencode,1,0)) # i.e. events before censored obs if occurring at the same time
  predictor <- predictor[ftime2.order]
  ftime <- ftime[ftime2.order]
  fstatus <- fstatus[ftime2.order]
  all.unique.failtimes <- sort(unique(ftime[fstatus==failcode]))
  cIndex <- function(predictor,ftime,fstatus,failcode,cencode){
    #- calculate cIndex (assuming that ftime is properly ordered and a variable all.unique.failtimes exists, see above)
    # which unique failure time in the full sample does each failure time belong to (only needed if predictor is a matrix to locate correct column)
    if (is.matrix(predictor)){
      failtimes <- ftime[fstatus==failcode]
      pos <- match(failtimes,all.unique.failtimes)
      j <- 1
    }
    # loop through events and calculate number of relevant and number of concordant pairs
    rel.pairs <- 0
    conc.pairs <- 0
    for (i in 1:n){
      if (fstatus[i]==failcode){
        # treat tied events and equal predictor values as in rcorr.cens of library(Hmisc)
        tied.events <- (ftime[i]==ftime[(i+1):n])&(fstatus[(i+1):n]==failcode)
        rel.pairs <- rel.pairs+(n-i)-sum(tied.events)
        if (is.vector(predictor)){
          conc.pairs <- conc.pairs+sum((predictor[i]>predictor[(i+1):n])&(!tied.events))+
                                   0.5*sum((predictor[i]==predictor[(i+1):n])&(!tied.events))
        } else {
          conc.pairs <- conc.pairs+sum((predictor[i,pos[j]]>predictor[(i+1):n,pos[j]])&(!tied.events))+
                                   0.5*sum((predictor[i,pos[j]]==predictor[(i+1):n,pos[j]])&(!tied.events))
          j <- j+1
        }
      }
    }
    # final result
    rel.pairs <- rel.pairs*2
    conc.pairs <- conc.pairs*2
    list(rel.pairs=rel.pairs,conc.pairs=conc.pairs,cIndex=conc.pairs/rel.pairs)
  }
  solution <- cIndex(predictor,ftime,fstatus,failcode,cencode)
  if (conf.int){
    boot.cIndex <- rep(NA,nBoot)
    for (k in 1: nBoot){
      boot.draw <- sort(sample(1:n,n,replace=T))
      if (is.vector(predictor)){
        boot.cIndex[k] <- cIndex(predictor[boot.draw],
                                 ftime[boot.draw],fstatus[boot.draw],failcode,cencode)$cIndex
        } else {
        boot.cIndex[k] <- cIndex(predictor[boot.draw,],
                                 ftime[boot.draw],fstatus[boot.draw],failcode,cencode)$cIndex
      }
    }
    boot.se <- sd(boot.cIndex)
    CI.lower <- 2*solution$cIndex-quantile(boot.cIndex,0.975)
    CI.upper <- 2*solution$cIndex-quantile(boot.cIndex,0.025)
  }
  if (!conf.int)  result <- c(solution$cIndex,NA,NA,NA,solution$rel.pairs,solution$conc.pairs)
  if (conf.int) result<- c(solution$cIndex,CI.lower,CI.upper,boot.se,solution$rel.pairs,solution$conc.pairs)
  names(result) <- c("cIndex","CI.lower","CI.upper","boot.se","relevant Pairs","concordant Pairs")
  result
}

rcorr.mycrr <- function(linear.predictor,ftime,fstatus,failcode=1,cencode=0,print=T){
  ##- Adapted c index for competing risk data (now largely superseded by function cIndex.CR) 
  ##- 
  ##- Author: Marcel Wolbers, 
  ##- Last update: Oct 17, 2007 (limited testing only, use at your own risk!)
  ##--------------------------------------------
  require(Design)
  # "classical" c index - censor competing events at the time they occur
  c.classical <- rcorr.cens(-linear.predictor,Surv(ftime,fstatus==failcode))[-c(2:3)]
  relPairs.classical <- c.classical["Relevant Pairs"]
  concordant.classical <- c.classical["Concordant"]
  # "new" c index for competing risks - censor competing events at infinity or, equivalently, at max(ftime)+1)
  ftime[(fstatus!=failcode)&(fstatus!=cencode)] <- max(ftime)+1
  c.compRisk <- rcorr.cens(-linear.predictor,Surv(ftime,fstatus==failcode))[-c(2:3)]
  relPairs.compRisk <- c.compRisk["Relevant Pairs"]
  concordant.compRisk <- c.compRisk["Concordant"]
  solution <- c(c(c.compRisk)["C Index"],c.compRisk["n"],sum(fstatus==failcode),sum((fstatus!=failcode)&(fstatus!=cencode)),
    relPairs.compRisk,relPairs.classical,relPairs.compRisk-relPairs.classical,c.compRisk["Uncertain"],
    concordant.compRisk,concordant.classical,concordant.compRisk-concordant.classical)
  names(solution) <- c("C index","n","# ev. of interest","# competing ev.",
                       "Relevant Pairs (total)","Relevant Pairs (Classical)","Additional relevant pairs","Uncertain pairs",
                       "Concordant (total)", "Concordant (classical)","Additional concordant pairs") 
  if (print==T) { 
    print(solution[1:4]); cat("\n")
    print(solution[5:8]); cat("\n")
    print(solution[9:11])
  } 
  invisible(solution)
}

D.mycrr <- function(linear.predictor,ftime,fstatus,failcode=1,cencode=0,cengroup=1,gtol=1e-06, maxiter=10){
  ##- Adapted Royston&Sauerbrei's D for competing risk data
  ##-
  ##- Author: Marcel Wolbers 
  ##- Last update: June 5, 2008 
  ##- (corrected scaling factor and added averaging over tied values)
  ##--------------------------------------------
  require(cmprsk)
  n <- length(linear.predictor)
  # sort data according to increasing values of linear predictor
  xord <- order(linear.predictor)
  linear.predictor <- linear.predictor[xord] 
  ftime <- ftime[xord]
  fstatus <- fstatus[xord]
  # calculate z
  kappa <- sqrt(8/pi)
  z <- kappa^(-1)*qnorm(ppoints(n))
  #avarage z over tied values of linear.predictor (if they exist)
  if (length(linear.predictor)!=length(unique(linear.predictor))){
    z.new <- z
    linear.predictor.level <- match(linear.predictor,unique(linear.predictor))
    for (i in 1:length(linear.predictor)) z.new[i] <- mean(z[linear.predictor.level==linear.predictor.level[i]])
    z <- z.new
  }
  #calculate D
  D <- crr(ftime, fstatus, cov1=z,failcode=failcode,cencode=cencode,cengroup=cengroup,gtol=gtol,maxiter=maxiter)
  c("D"=D$coef,"se"=sqrt(D$var),"CI.lower"=D$coef-qnorm(0.975)*sqrt(D$var),"CI.upper"=D$coef+qnorm(0.975)*sqrt(D$var))
}