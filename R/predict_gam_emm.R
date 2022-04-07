



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
predict_gam_emm <- function(fit, d, quantile_diff = c(0.25, 0.75), Xvar, Yvar, binaryX = FALSE){

  set.seed(12345)
  require(mgcv)
  require(dplyr)
  d$dummy <- 0

  Wvars <- colnames(d)[!(colnames(d) %in% c("Y", "X", "V",
                                            "id", "dummy"))]
  for (i in Wvars) {
    if (class(d[, i]) == "character" | class(d[, i]) ==
        "factor") {
      d[, i] <- Mode(d[, i])
    }
    else {
      d[, i] <- median(d[, i])
    }
  }
  d <- d[order(d$X), ]


  if (binaryX == F) {
    q1 <- unname(quantile(d$X, quantile_diff[1]))
    q3 <- unname(quantile(d$X, quantile_diff[2]))
    q1_pos <- which(abs(d$X - q1) == min(abs(d$X - q1)))[1]
    q3_pos <- which(abs(d$X - q3) == min(abs(d$X - q3)))[1]
    d$X[q1_pos] <- q1
    d$X[q3_pos] <- q3
  }
  if (binaryX == T) {
    q1 <- min(d$X)
    q3 <- max(d$X)
    q1_pos <- 1
    q3_pos <- nrow(d)
    d$X[q1_pos] <- min(d$X)
    d$X[q3_pos] <- max(d$X)
  }

  #Add specific rows for the quartile predictions
  Nrows <- nrow(d)
  d <- bind_rows(d, d[1:2,])
  d$X[c(Nrows+1,Nrows+2)] <- c(q1, q3)

  dfull <- d
  plotdf <- res <- NULL

  #Detect if modifier is
  if(class(dfull$V)!="numeric"){
    for(i in unique(dfull$V)){
      d <- dfull
      d$V <- i
      preds <- predict(fit, newdata = d, type = "response")
      Xp <- predict(fit, newdata = d, type = "lpmatrix")
      #Xp <- Xp[order(d$X), ] #caitlin debug fix... check
      diff <- t(apply(Xp, 1, function(x) x - Xp[Nrows+1, ]))
      point.diff <- diff %*% coef(fit)
      se.diff <- sqrt(diag(diff %*% vcov(fit) %*% t(diff)))
      lb.diff <- point.diff - 1.96 * se.diff
      ub.diff <- point.diff + 1.96 * se.diff
      Zval <- abs(point.diff/se.diff)
      Pval <- exp(-0.717 * Zval - 0.416 * Zval^2)
      resdf <- data.frame(Y = Yvar, X = Xvar, Vlevel=i, N = nrow(d), q1 = d$X[q1_pos],
                          q3 = d$X[q3_pos], pred.q1 = preds[q1_pos], pred.q3 = preds[q3_pos],
                          point.diff, lb.diff = lb.diff, ub.diff = ub.diff, Pval = Pval)

      # if(binaryX==T){
      #   temp_res <- resdf[1, ]
      # }else{
        temp_res <- resdf[nrow(resdf), ]
      #}


      plotdf <- bind_rows(plotdf, resdf)

      res <- bind_rows(res, temp_res)
    }

  }else{

    quartiles <- as.numeric(summary(dfull$V)[c(2,5)])

    for(i in unique(quartiles)){
      d <- dfull
      d$V <- i
      preds <- predict(fit, newdata = d, type = "response")
      Xp <- predict(fit, newdata = d, type = "lpmatrix")
      #Xp <- Xp[order(d$X), ]
      diff <- t(apply(Xp, 1, function(x) x - Xp[Nrows+1, ]))
      point.diff <- diff %*% coef(fit)
      se.diff <- sqrt(diag(diff %*% vcov(fit) %*% t(diff)))
      lb.diff <- point.diff - 1.96 * se.diff
      ub.diff <- point.diff + 1.96 * se.diff
      Zval <- abs(point.diff/se.diff)
      Pval <- exp(-0.717 * Zval - 0.416 * Zval^2)
      resdf <- data.frame(Y = Yvar, X = Xvar, Vlevel=i, N = nrow(d), q1 = d$X[q1_pos],
                          q3 = d$X[q3_pos], pred.q1 = preds[q1_pos], pred.q3 = preds[q3_pos],
                          point.diff, lb.diff = lb.diff, ub.diff = ub.diff, Pval = Pval)

      #if(binaryX == T){
      temp_res <- resdf[nrow(resdf), ]
      # }else{
      #   temp_res <- resdf[q3_pos, ]
      # }

      plotdf <- bind_rows(plotdf, resdf[1:Nrows,])

      res <- bind_rows(res, temp_res)
    }

  }

  return(list(res = res, plotdf = plotdf))
}
