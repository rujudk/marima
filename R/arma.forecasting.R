
##'
##' @title Forecasting of (multivariate) time series
##' of arma-type.
##'
##' @param series = matrix holding the k-variate timeseries.
##' The series is assumed to have the same format 
##' as the timeseries analysed by marima (the length, though, does
##' not need to be the same but can be shorter or longer). Results
##' from estimating the model are assumed to be saved in the
##' object 'marima' by, say, marima. 
##'
##' The series is assumed to have the total length=(nstart+nstep) (but it
##' may be longer. In any case the forecasting is starting from nstart
##' continuing to nstart+nstep. Future values already present or initialised,
##' for example, as NAs are overwritten with the forecasted values.)
##' An example of a series prepared for forcasting is in the marima library:
##' 'data(australian.killings)' 'austr' (see below, the example).  
##'
##' If future (independent) x-values for the forecasting are to be used
##' these values must be supplied in 'series' at the proper places before
##' calling 'forecast(...)'
##' 
##' @param marima = object holding the marima results to be used for the
##' forecasting, that is an output object from marima.
##' 
##' If the ar- and/or the ma-model do not include a leading unity matrix
##' this is automatically taken care of in the function (in that case the
##' dimensions of the model arrays used will be, respectively,
##' (k,k,p+1) and (k,k,q+1)) after inserting the leading unity matrix (if
##' the object 'Marima' was produced by marima, this will automatically be
##' OK.
##'
##' @param nstart = starting point for forecasting (1st forecast values
##' will be for time point t = nstart+1).
##'
##' @param nstep = length of forecast (forecasts will be for time points
##' nstart+1,...,nstart+nstep).
##'
##' @return forecasts = forecasted values following the input series.
##' The forecasted values will be (over-) written in the input series at
##' the proper future positions (if relevant)
##'
##' @return residuals = corresponding residuals for input series followed by
##' nstep future residuals (all=0) 
##'
##' @return prediction.variances = (k,k,nstep) array containing prediction
##' covariance matrices corresponding to the nstep forecasts
##'
##' @return nstart = starting point for prediction (1st prediction at point
##' nstart+1
##'
##' @return nstep = length of forecast
##'  
##' @examples
##' 
##' library(marima)
##' data(australian.killings)
##' series<-t(austr)
##' Model5 <- define.model(kvar=7,ar=1,ma=1,rem.var=1,reg.var=6:7)
##' Marima5 <- marima(series[,1:90],Model5$ar.pattern,Model5$ma.pattern,
##' penalty=1)
##' nstart  <- 90
##' nstep   <- 10
##' l       <- 91
##' Forecasts <- arma.forecast(series=series,marima=Marima5,
##'                nstart=nstart,nstep=nstep)
##' Year<-series[1,91:100];
##' Predict<-Forecasts$forecasts[2,91:100]
##' stdv<-sqrt(Forecasts$pred.var[2,2,])
##' upper.lim=Predict+stdv*1.645
##' lower.lim=Predict-stdv*1.645
##' Out<-rbind(Year,Predict,upper.lim,lower.lim)
##' print(Out)
##' # plot results:
##' plot(series[1,1:100],Forecasts$forecasts[2,],type="l",xlab="Year",
##' ylab="Rate of armed suicides",main="Prediction of suicides by firearms",
##' ylim=c(0.0,4.1))
##' lines(series[1,1:90],series[2,1:90],type="p")
##' grid(lty=2,lwd=1,col="black")
##' Years<-2005:2014
##' lines(Years,Predict,type="l")
##' lines(Years,upper.lim,type="l")
##' lines(Years,lower.lim,type="l")
##'

arma.forecast <- function(series=NULL,marima=NULL,
        nstart=NULL, nstep=1 )
{
# "["<-function(x,...).Primitive("[")(x,...,drop=FALSE)    
cat("Length of prediction series =",c(nstart+nstep)," \n")
Y  <- series[,c(1:(nstart+nstep))]
d  <- dim(Y)
if(d[1]>d[2]){ Y <- t(Y) ; d<-dim(Y) }

colnames(series)<-c(1:d[2])

means     <- marima$mean.pattern
averages  <- marima$averages
# cat("means=",means,averages,"\n")

for (i in (1:d[1])){Y[i,] <- Y[i,]-averages[i]*means[i]}
# cat("Means adjusted input, last 20 obs: \n")
# print(Y[,(d[2]-19):d[2]])

if(is.null(nstart)){nstart<-d[2]}
cat("nstart=",nstart,"\n")
Names<-rownames(marima$y.analysed)
# cat("Names=",Names,"\n")
AR <- marima$ar.estimate
MA <- marima$ma.estimate
q<-dim(AR)[3]

# print(AR)
# print(MA)

Rand <- rep(TRUE,d[1])
Regr <- !Rand
sig2<-marima$resid.cov

# print(sig2,digits=2)

for (i in 1:d[1]){
   if(sum(abs(AR[i,,]))+sum(abs(MA[i,,]))==2){
       Rand[i]<-FALSE
       sig2[i,]<-0
       sig2[,i]<-0 }
   if(sum(abs(AR[,i,]))+sum(abs(MA[,i,]))> 2){Regr[i]<-TRUE}}

# Rand = TRUE=  random variable, FALSE = deterministic variable
# Regr = TRUE= 'regressor' in model, FALSE = not 'regressor' in model.
# Variables = (FALSE,TRUE) must appear in input series (future values)

# cat("Rand",Rand,"\n")
# cat("Regr",Regr,"\n")

rownames(Y)<-Names
if(is.null(Names)){Names<-paste("Y",c(1:d[1]),sep="")}

nstart=nstart
nstep=nstep

L<-nstart+nstep
if(L>d[2]){cat("Input series length =",d[2],
                    " but (nstart+nstep)=",L,"\n")}
for ( i in 1:d[1]){
#    cat(Rand[i],Regr[i],"\n")
    if(!Rand[i]&Regr[i]){
        cat("variable no.",i,"not random and regressor.\n")
        zy<-Y[i,(nstart+1):L]
#        zy[5]<-NaN
        if(is.na(sum(zy))|is.nan(sum(zy))){
            cat("Variable no.",i,": illegal (NaN or NA) in future
            values. \n")}
            if(!is.na(sum(zy))&!is.nan(sum(zy)) ){
 cat("Variable no.",i,": future values seem OK.\n")
             } } }
# cat("Calling series : \n")
# print(Y[,(d[2]-19):d[2]])

forecast <- arma.forecastingX(series=Y,ar.poly=AR,ma.poly=MA,
   nstart=nstart,nstep=nstep,Rand=Rand,Regr=Regr,marima=marima)

pred.var<-forec.var(marima,nstep=nstep)$pred.var
 
return(list(nstart=nstart,nstep=nstep,forecasts=forecast$estimates,
        residuals=forecast$residuals,pred.var=pred.var,averages=averages,
            mean.pattern=means)) }

arma.forecastingX <- function(series=NULL,ar.poly=NULL,ma.poly=NULL,
  nstart=NULL,nstep=NULL,Rand=NULL,Regr=NULL,marima=NULL) {
"["<-function(x,...).Primitive("[")(x,...,drop=FALSE)
d <- dim(series)
if(d[1]>d[2]) { series<-t(series) }
d<-dim(series)
k<-d[1]

means<-marima$mean.pattern
averages<-marima$averages
cat("dimensions of series",d,"\n")

# for (i in (1:d[1])) {
# if  (means[i]==1)   {
#   series[i,]<-series[i,]-averages[i]
# }}
# cat("Calling series (2): \n")
# print(series[,(d[2]-19):d[2]])

 k<-dim(series)[1]
 N<-dim(series)[2]
# cat("k,N=",k,N,"\n")

if(is.null(ar.poly))ar.poly<-array(diag(d[1]),0*diag(d[1]),dim=c(k,k,2))
if(is.null(ma.poly))ma.poly<-array(diag(d[1]),0*diag(d[1]),dim=c(k,k,2))   
  ar.poly<-check.one(ar.poly)
  ma.poly<-check.one(ma.poly)
ar  <- dim(ar.poly)
ma  <- dim(ma.poly)
# cat(" AR= \n")
# print(ar)
# cat(" MA= \n")
# print(ma)
extra <- max(max(ar[3]),max(ma[3]))
# cat("min. forecast previous length =",extra,"\n")
yseries   <- matrix(c(series),nrow=d[1])
# yseries[5,]<-NaN
# cat("yseries= \n")
# print(yseries[,80:dim(yseries)[2]])

yseries[is.na(yseries)|is.nan(yseries)]<-0
residuals <- yseries*0
forecasts <- residuals
# print(residuals[,80:dim(yseries)[2]])
# s <- d
kvar<-dim(ar.poly)[1]
su0<-matrix(0,nrow=kvar,ncol=1)
# cat("su0= ",su0,"\n")

for (i in (extra:d[2])){  # Start i-loop here 
#    cat("i=",i,"\n")
suar<-su0
 if (ar[3]>1) {
     for (j in 2:ar[3]){
suar <- suar+
matrix(ar.poly[,,j],nrow=kvar)%*%matrix(yseries[,(i+1-j)],ncol=1)
if(i<0) cat("i,suar= ",i,suar,"\n")
}}
suma<-su0
 if (ma[3]>1) {
      for (j in 2:ma[3]){
suma <-
suma+matrix(ma.poly[,,j],nrow=kvar)%*%matrix(residuals[,(i+1-j)],ncol=1)
if(i<0) cat("i,suma= ",i,suma,"\n") } }

est <- suma - suar
forecasts[,i]<-est


# cat("ar- & ma-polynomier \n")
#     print(ar.poly[,,2])
#     print(ma.poly[,,2])

if (i<=nstart) { for (j in 1:k)
        { if(Rand[j]) {residuals[j,i] <- yseries[j,i]-est[j] } } }
if (i> nstart) { for (j in 1:k)
        { if(Rand[j]) {yseries[j,i]  <- est[j]
                                  residuals[j,i] <- 0 } } }

} # End i-loop here

 # cat("Estimates-88:90\n")
 # print(yseries[,88:90])
# cat("means=",means,"\n")
for (i in (1:kvar)){ if(means[i]==1)
 {yseries[i,]   <- yseries[i,]+averages[i]
  forecasts[i,] <- forecasts[i,]+averages[i] }}
# cat("y-88:90\n")
# print(yseries[,88:90])
# cat("res-88:90\n")
# print(residuals[,88:90])

return(list(estimates=forecasts,residuals=residuals)) }

##'
##' @title Function for calculation of variances of nstep forecasts
##'
##' @param marima   =  marima object (resid.cov and ar.estimates and
##'                    ma.estimates are used)
##' @param nstep    =  length of forecast
##'
##' @return = pred.var   = variance-covariances for nstep forecasts
##' @return = rand.shock = random shock representation of the marima model
##'

forec.var <- function(marima,nstep=1){
 # "["<-function(x,...).Primitive("[")(x,...,drop=FALSE)   
cat("Calculation forecasting variances.  \n")
sig2<-marima$resid.cov
ar.poly<-marima$ar.estimates
ma.poly<-marima$ma.estimates

kvar=dim(sig2)[1]

d<-dim(sig2)
if(is.null(ar.poly)){ar.poly<-diag(d[1])}
if(is.null(ma.poly)){ma.poly<-diag(d[1])} 
 ar.poly<-check.one(ar.poly)
 ma.poly<-check.one(ma.poly)
if(nstep<1){nstep<-1}
 xsi.poly<-rand.shock(ar.poly,ma.poly,nstep)
var<-xsi.poly
# cat("nsteps=",1:nstep,"\n")
for (i in 1:nstep){
#    cat("i=",i,"\n")
   var[,,i]<-
matrix(xsi.poly[,,i],nrow=kvar)%*%sig2%*%t(matrix(xsi.poly[,,i],nrow=kvar))
   if (i>1){var[,,i]<-var[,,i]+var[,,(i-1)]}}
d<-dim(var)[3]-1
# cat("d=",d,"\n")
# cat("xsi.poly=",1:nstep,"\n")
# print(round(xsi.poly,4))
# cat("var=","\n")
# print(round(var[,,1:d],4))
return(list(pred.var=var[,,1:d],rand.shock=xsi.poly[,,1:d])) }


##'
##' @title Filtering of (multivariate) time series with arma-type model.
##'
##' Calculation of forecasts and residuals for timeseries and
##' differencing of time series (used in connection with routine define.dif).
##' 
##' @param series   = matrix holding the k by n multivariate timeseries
##' (if k>n the series is transposed and a warning is given).
##'
##' @param ar.poly = (k,k,p+1) array containing autoregressive matrix
##' polynomial model part.
##'
##' @param ma.poly = (k,k,q+1) array containing moving average matrix
##' polynomial model part.
##'
##' If a leading unity matrix is not included in the ar- and/or the ma-part
##' of the model this is automatically taken care of in the function 
##' (in that case the dimensions of the model arrays used in arma.filter() 
##'  are, respectively, (k,k,p+1) and (k,k,q+1)).
##' 
##' @return estimates = estimated values for input series 
##'
##' @return residuals = corresponding residuals
##'
##' Both estimates and residuals are organised as k by n matrices.
##'
##' @examples
##'
##' library(marima)
##' data(australian.killings)
##' series<-t(austr)[,1:90]
##' # Define marima model
##' Model5 <- define.model(kvar=7,ar=1,ma=1,rem.var=1,reg.var=6:7)
##' 
##' # Estimate marima model
##' Marima5 <- marima(series,Model5$ar.pattern,Model5$ma.pattern,penalty=1)
##' 
##' # Calculate residuals by filtering
##' Resid <- arma.filter(series,Marima5$ar.estimates,
##'      Marima5$ma.estimates)
##' # Compare residuals
##'
##' plot(Marima5$residuals[2,4:89],Resid$residuals[2,5:90],
##' xlab="marima residuals", ylab="arma.filter residuals")
##'        

arma.filter <- function(series=NULL,
                  ar.poly=array(diag(kvar),dim=c(kvar,kvar,1)),
    ma.poly=array(diag(kvar),dim=c(kvar,kvar,1)), means=1 )
{
"["<-function(x,...).Primitive("[")(x,...,drop=FALSE)
cat("Start of arma.filter \n")
 if(is.null(series)){
    stop("No input series specified in input to 'arma.filter' \n") }
d <- dim(series)
if(d[1]>d[2]) {
    cat("Warning: Input series is transposed to be a k x n series. \n")
    cat("Output (estimated values and residuals) \n")
    cat("will be organised the same way (k x n). \n")
series<-t(series) }

vmeans<-means
cat("vmeans,means=",vmeans,means,"\n")
d<-dim(series)
kvar<-d[1]

averages<-rep(0,kvar)
means<-rep(1,kvar)
if(length(vmeans==1)){
if(vmeans!=1){means<-rep(0,kvar)}}
for (i in 1:kvar){
 averages[i]<-mean(series[i,])
 if (means[i]==1){series[i,]<-series[i,]-averages[i]}}
 cat("means=",means,"\n")
 cat(averages,"\n")

 ar.poly<-check.one(ar.poly)
#  cat("ar.poly=")
# print(round(short.form(ar.poly,leading=F),4))
 ma.poly<-check.one(ma.poly)
#  cat("ma.poly=")
# print(round(short.form(ma.poly,leading=F),4))

m   <- dim(ar.poly)
ma  <- dim(ma.poly)

su0 <- matrix(0,nrow=kvar,ncol=1)
extra <- 2*max(m[3],ma[3])

extray <- matrix(0,nrow=kvar,ncol=extra)
for (i in 1:kvar){if(means[i]!=1){extray[i,]<-extray[i,]+averages[i]}}

# cat("extra=",extra,"\n")
# print(extray)
 
yseries   <- cbind(extray,series)
estimates <- yseries*0
residuals <- estimates
 
# cat("Series= \n")
#  print(averages)
#  print(round(yseries[,1:20],2))
#  print(round(short.form(ar.poly,leading=F),2))

s  <- dim(yseries)
cat(" dim(yseries)",s,"\n")

for (i in (extra/2+1):s[2]){
sur<-su0
 if(m[3]>1){
 for (j in 2:m[3]){
#  cat("sur=",sur,"\n")
#  cat("m,i,j=",m,i,j,"\n")
#  cat("dimensions=",dim(ar.poly[,,j]),dim(yseries[,(i+1-j)]),"\n")
#  print(ar.poly[,,j])
#  print(yseries[,(i+1-j)])
 sur <- sur +
 matrix(ar.poly[,,j],nrow=kvar)%*%matrix(yseries[,(i+1-j)],nrow=kvar) } }

#  cat(i,j,(i+1-j),sur,"\n")
#  cat("y=",yseries,"\n")

 suma<-su0
 if(ma[3]>1){
 for (j in 2:ma[3]){
 suma<-suma +
 matrix(ma.poly[,,j],nrow=kvar)%*%matrix(residuals[,(i+1-j)],nrow=kvar) } }
    
#  cat(i,j,(i+1-j),suma,"\n")
 su <- -sur +suma
 estimates[,i] <- su
 residuals[,i] <- yseries[,i]-estimates[,i]}

estimates  <- estimates[,(extra+1):(extra+d[2])]
residuals  <- residuals[,(extra+1):(extra+d[2])]

for (i in 1:kvar){if(means[i]==1){estimates[i,]<-estimates[i,]+averages[i]}}

# cat("dimensions",dim(yseries),dim(series),dim(residuals),"\n")

return(list(estimates=estimates,residuals=residuals,averages=averages,
             means=means)) }


