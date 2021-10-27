
#' Title
#'
#' @param m
#' @param newdata
#' @param nreps
#' @param xlab
#' @param ylab
#' @param title
#'
#' @return
#' @export
#'
#' @examples
gam_simul_CI <- function(m,newdata,nreps=10000, xlab="", ylab="", title="", gam_diff=NULL) {
  set.seed(12345)
  require(mgcv)
  require(dplyr)

  newdata <- newdata %>% mutate(dummy=0)

  Wvars <- colnames(newdata)[!(colnames(newdata) %in% c("Y","X" ,"id" ,"dummy"))]
  #set covariates to the median/mode
  for(i in Wvars){
    if(class(newdata[,i])=="character"|class(newdata[,i])=="factor"){
      newdata[,i] <- Mode(newdata[,i])
    }else{
      newdata[,i] <- median(newdata[,i])
    }
  }

  newdata <- newdata[order(newdata$X),]

  Vb <- vcov(m,unconditional = TRUE)
  pred <- predict(m, newdata, se.fit = TRUE)
  fit <- pred$fit
  se.fit <- pred$se.fit
  BUdiff <- MASS::mvrnorm(n=nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
  Cg <- predict(m, newdata, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  pred <- data.frame(newdata,fit=pred$fit,se.fit=pred$se.fit)
  pred <- mutate(pred,
                 uprP = fit + (2 * se.fit),
                 lwrP = fit - (2 * se.fit),
                 uprS = fit + (crit * se.fit),
                 lwrS = fit - (crit * se.fit)
  )

  pred <- pred %>% arrange(X)
  p <- ggplot(pred) + geom_ribbon(aes(x=X, ymin=lwrS, ymax=uprS), alpha=0.5) +
    geom_path(aes(x=X, y=lwrS), color="blue")+
    geom_path(aes(x=X, y=uprS), color="red")+
    geom_path(aes(x=X, y=fit ), color="black") +
    xlab(xlab) + ylab(ylab) +
    ggtitle(title)

  if(!is.null(gam_diff)){
    p <- p +
      geom_vline(xintercept = gam_diff$q1) +
      geom_vline(xintercept = gam_diff$q3) +
      ggtitle(paste0(round(gam_diff$point.diff,2)," (",
                     round(gam_diff$lb.diff,2), ", ",
                     round(gam_diff$ub.diff,2), ")"))
  }

  return(list(p=p, pred=pred))
}






