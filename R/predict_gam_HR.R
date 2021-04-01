



#' Title
#'
#' @param fit
#' @param d
#' @param quantile_diff
#' @param Xvar
#' @param Yvar
#' @param binaryX
#'
#' @return
#' @export
#'
#' @examples

predict_gam_HR <- function(fit, d, quantile_diff=c(0.25,0.75), Xvar, Yvar, binaryX=FALSE){
  set.seed(12345)
  require(mgcv)
  require(dplyr)

  d$dummy<-0

  Wvars <- colnames(d)[!(colnames(d) %in% c("Y","X" ,"id" ,"dummy"))]
  #set covariates to the median/mode
  for(i in Wvars){
    if(class(d[,i])=="character"|class(d[,i])=="factor"){
      d[,i] <- Mode(d[,i])
    }else{
      d[,i] <- median(d[,i])
    }
  }

  d <- d[order(d$X),]

  if(binaryX==F){
    #Make sure subset has overall quantiles within it
    q1 <- unname(quantile(d$X,quantile_diff[1]))
    q3 <- unname(quantile(d$X,quantile_diff[2]))

    q1_pos <- which(abs(d$X- q1)==min(abs(d$X- q1)))[1]
    q3_pos <- which(abs(d$X- q3)==min(abs(d$X- q3)))[1]
    d$X[q1_pos] <- q1
    d$X[q3_pos] <- q3
  }

  if(binaryX==T){
    q1 <- min(d$X)
    q3 <- max(d$X)
    q1_pos <-1
    q3_pos <- nrow(d)
    d$X[q1_pos] <- min(d$X)
    d$X[q3_pos] <- max(d$X)
    #Note, can I just grab the coefficient on X from the "fit" regression object?
  }

  #get the direct prediction
  preds <- predict(fit,newdata=d,type="response")

  #get the prediction matrix
  try(Xp <- predict(fit,newdata=d,type="lpmatrix"))
  # order the prediction matrix in the order of the exposure
  Xp <- Xp[order(d$X),]



  # take difference from the 25th percentile of X
  diff <- t(apply(Xp,1,function(x) x - Xp[q1_pos,]))

  # calculate the predicted HR's
  point.diff <- diff %*% coef(fit)
  point.HR <- exp(point.diff)

  # calculate the pointwise SE - naive SE
  se.diff <- sqrt(diag( diff%*%vcov(fit)%*%t(diff) ) )

  # calculate upper and lower bounds
  lb.HR <- exp(point.diff - 1.96*se.diff)
  ub.HR <- exp(point.diff + 1.96*se.diff)
  Zval <-  abs(point.diff/se.diff)
  Pval <- exp(-0.717*Zval - 0.416*Zval^2)

  plotdf<-data.frame(Y=Yvar, X= Xvar, N=nrow(d), q1=d$X[q1_pos], q3=d$X[q3_pos],
                     point.HR=point.HR, lb.HR=lb.HR, ub.HR=ub.HR, Pval=Pval)

  if(binaryX==T){
    res <- plotdf[nrow(plotdf),]
  }else{
    #res <- plotdf[round(nrow(d)*quantile_diff[2],0),]
    res <- plotdf[q3_pos,]
  }

  return(list(res=res, plotdf=plotdf))
}
