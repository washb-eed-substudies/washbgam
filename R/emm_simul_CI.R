

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

emm_simul_CI <- function (m, newdata, nreps = 10000, xlab = "", ylab = "",
                          title = ""#, gam_diff = NULL
){
  set.seed(12345)
  require(mgcv)
  require(dplyr)
  newdata <- newdata %>% mutate(dummy = 0)
  Wvars <- colnames(newdata)[!(colnames(newdata) %in% c("Y", "X", "V", "id", "dummy"))]
  for (i in Wvars) {
    if (class(newdata[, i]) == "character" | class(newdata[,
                                                           i]) == "factor") {
      newdata[, i] <- Mode(newdata[, i])
    }
    else {
      newdata[, i] <- median(newdata[, i])
    }
  }
  newdata <- newdata[order(newdata$X), ]
  newdata$V <- as.character(newdata$V)
  newdata2 <- newdata1 <- newdata
  newdata1$V <- unique(newdata$V)[1]
  newdata2$V <- unique(newdata$V)[2]

  Vb <- vcov(m, unconditional = TRUE)
  pred1 <- predict(m, newdata1, se.fit = TRUE)
  pred2 <- predict(m, newdata2, se.fit = TRUE)
  fit1 <- pred1$fit
  fit2 <- pred2$fit
  se.fit1 <- pred1$se.fit
  se.fit2 <- pred2$se.fit
  BUdiff <- MASS::mvrnorm(n = nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
  Cg1 <- predict(m, newdata1, type = "lpmatrix")
  Cg2 <- predict(m, newdata2, type = "lpmatrix")
  simDev1 <- Cg1 %*% t(BUdiff)
  simDev2 <- Cg2 %*% t(BUdiff)
  absDev1 <- abs(sweep(simDev1, 1, se.fit1, FUN = "/"))
  absDev2 <- abs(sweep(simDev2, 1, se.fit2, FUN = "/"))
  masd1 <- apply(absDev1, 2L, max)
  masd2 <- apply(absDev2, 2L, max)
  crit1 <- quantile(masd1, prob = 0.95, type = 8)
  crit2 <- quantile(masd2, prob = 0.95, type = 8)
  pred1 <- data.frame(newdata1, fit = pred1$fit, se.fit = pred1$se.fit)
  pred2 <- data.frame(newdata2, fit = pred2$fit, se.fit = pred2$se.fit)
  pred1 <- mutate(pred1, uprP = fit + (2 * se.fit1), lwrP = fit -
                    (2 * se.fit1), uprS = fit + (crit1 * se.fit1), lwrS = fit -
                    (crit1 * se.fit1)) %>% arrange(X) %>% mutate(Vlevel=unique(newdata$V)[1])
  pred2 <- mutate(pred2, uprP = fit + (2 * se.fit2), lwrP = fit -
                    (2 * se.fit2), uprS = fit + (crit2 * se.fit2), lwrS = fit -
                    (crit2 * se.fit2)) %>% arrange(X) %>% mutate(Vlevel=unique(newdata$V)[2])
  pred <- bind_rows(pred1, pred2)
  p <- ggplot(pred) + geom_ribbon(aes(x = X, ymin = lwrS, ymax = uprS, group=V, fill=V, color=V),
                                  alpha = 0.5) +
    # geom_path(aes(x = X, y = lwrS, group=V, color=V)) +
    # geom_path(aes(x = X, y = uprS, group=V, color=V)) +
    geom_path(aes(x = X, y = fit, group=V, color=V)) +
    xlab(xlab) + ylab(ylab) + ggtitle(title)

  return(list(p = p, pred = pred))
}
