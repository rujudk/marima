print.marima <-
function(obj,estimates=TRUE,pvalues=FALSE,pattern=TRUE,fvalues=TRUE){
  cat("Model dimension = ",obj$kvar, " N = ",obj$N," \n")
 # if(pattern==TRUE){
 #     cat("ar.pattern = \n")
 #     print(short.form(obj$call.ar.pattern,leading=F))
 #     cat("ar.pattern = \n")
 #     print(short.form(obj$call.ma.pattern,leading=F))}

    cat("Averages for all variales:\n",obj$averages,"\n")
    cat("Covariance(data):\n")
    print(round(obj$data.cov,4))
    cat("Covariance(residuals):\n")
    print(round(obj$resid.cov,4))
    cat("Random variables are",obj$randoms,"\n")
  
 if(pattern==TRUE){
    cat("\n")
    cat("AR definition:\n")
  { AR<-short.form(obj$call.ar.pattern,leading=FALSE)
    print(AR)}
    cat("MA definition:\n")
  { MA<-short.form(obj$call.ma.pattern,leading=FALSE)
    print(MA)}
  if(pattern!=TRUE){cat("No model identification output specified \n") }}
  
  if(estimates==TRUE) {cat("AR estimates:\n")
  print(round(short.form(obj$ar.estimates,leading=FALSE),4)) 
  if(fvalues==TRUE){
    cat("AR f-values (squared t-values):\n")
  print(round(short.form(obj$ar.fvalues,leading=FALSE),4))}
  if(pvalues==TRUE){
    cat("AR p-values (%):\n")
  print(round(100*short.form(obj$ar.pvalues,leading=FALSE),2))}
                    }
  if(estimates!=TRUE){cat("No estimated model output specified \n") }

  if(estimates==TRUE){cat("MA estimates:\n")
   print(round(short.form(obj$ma.estimates,leading=FALSE),4)) 
  if(fvalues==TRUE){
    cat("MA f-values (squared t-values):\n")
  print(round(short.form(obj$ma.fvalues,leading=FALSE),2))}
   if(pvalues==TRUE){
    cat("MA p-values (%):\n")
  print(round(100*short.form(obj$ma.pvalues,leading=FALSE),2))}}
  if(estimates!=TRUE){cat("No estimated model output specified \n") }

}

