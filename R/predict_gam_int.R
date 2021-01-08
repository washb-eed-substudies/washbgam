



predict_gam_int <- function(fit, d, quantile_diff=c(0.25,0.75), Xvar, Yvar){

  d$dummy<-0

  Wvars <- colnames(d)[!(colnames(d) %in% c("Y","X","V" ,"id" ,"dummy"))]
  #set covariates to the median/mode
  for(i in Wvars){
    if(class(d[,i])=="character"|class(d[,i])=="factor"){
      d[,i] <- Mode(d[,i])
    }else{
      d[,i] <- median(d[,i])
    }
  }

  d <- d[order(d$X),]

  q1 <- unname(quantile(d$X,quantile_diff[1]))
  q3 <- unname(quantile(d$X,quantile_diff[2]))

  reslist <- plotdf_list <- list()
  for(i in 1:length(unique(d$V))){
    diff <- NULL
    subgroup <- unique(d$V)[i]
    dsub <- d[d$V==subgroup,]

    #Make sure subset has overall quantiles within it
    q1_pos <- which(abs(dsub$X- q1)==min(abs(dsub$X- q1)))
    q3_pos <- which(abs(dsub$X- q3)==min(abs(dsub$X- q3)))
    dsub$X[q1_pos] <- q1
    dsub$X[q3_pos] <- q3
    dsub <- dsub[order(dsub$X),]

    Xp <- predict(fit,newdata=dsub,type="lpmatrix")

    # order the prediction matrix in the order of the exposure
    Xp <- Xp[order(dsub$X),]

    # take difference from the 25th percentile of X
    diff <- t(apply(Xp,1,function(x) x - Xp[q1_pos,]))

    # calculate the predicted differences
    point.diff <- diff %*% coef(fit)

    # calculate the pointwise SE - naive SE
    se.diff <- sqrt(diag( diff%*%vcov(fit)%*%t(diff) ) )

    # calculate upper and lower bounds
    lb.diff <- point.diff - 1.96*se.diff
    ub.diff <- point.diff + 1.96*se.diff

    plotdf<-data.frame(Y=Yvar, X= Xvar, subgroup=subgroup, x=dsub$X, q1=q1, q3=q3, point.diff, lb.diff=lb.diff, ub.diff=ub.diff)
    reslist[[i]] <- plotdf[q3_pos,]
    plotdf_list[[i]] <- plotdf
  }

  res <- bind_rows(reslist)

  return(list(res=res, plotdf=plotdf_list))
}
