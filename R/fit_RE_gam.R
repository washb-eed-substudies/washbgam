


#' Title
#'
#' @param d
#' @param Y
#' @param X
#' @param W
#' @param forcedW Set to NA to skip forcing any variables, don't include to skip prescreening of "tr" and variables starting with age_, ageday_, or agedays_. Or manually set with a vector of variable names.
#' @param V
#' @param id
#' @param family
#' @param pval
#' @param print
#'
#' @return
#' @export
#'
#' @examples

fit_RE_gam <- function(d, Y, X, W=NULL,
                       forcedW=NULL,
                       V=NULL, id="clusterid", family = "gaussian", pval = 0.2, print=TRUE){

  cat("\nNon-prescreened covariates: ", paste(forcedW, sep="", collapse=", "), "\n")
  #cat("Forced covariates:", forcedW,"\n\n")

  set.seed(12345)
  require(mgcv)
  require(dplyr)
  require(faraway)
  if(!is.null(V)){
    require(lmtest)
  }

  if(!is.null(W)){
    W <- subset(d, select = W)
  }
  Y <- subset(d, select = Y)
  colnames(Y) <- "Y"
  X <- subset(d, select = X)
  colnames(X) <- "X"
  id <- subset(d, select = id)
  colnames(id) <- "id"

  if(!is.null(V)){
    Vvar <- subset(d, select = V)
    colnames(Vvar) <- "V"
  }else{
    Vvar <-data.frame(V=rep(1, nrow(d)))
  }

  if(!is.null(W)){
    gamdat <- data.frame(Y, X, id, Vvar, W)
  }else{
    gamdat <- data.frame(Y, X, id, Vvar)
  }

  if(!is.null(W)){

    if(sum(is.na(forcedW))!=0){
      colnamesW <- names(W)
    }else{
      if(is.null(forcedW)){
        Wnames <- names(W)
        forcedW <- c(Wnames[Wnames=="tr"|grepl("age_", Wnames)|grepl("agedays_", Wnames)|grepl("ageday_", Wnames)])
      }
      colnamesW <- names(W)[!(names(W) %in% forcedW)]
    }
    #cat(names(W)[!(names(W) %in% forcedW)])
    screenW <- subset(gamdat, select = colnamesW)
  }else{
    screenW <- NULL
  }

  if(!is.null(screenW)){
    if(print == TRUE){
      cat("\n-----------------------------------------\nPre-screening the adjustment covariates:\n-----------------------------------------\n")
    }
    suppressWarnings(Wscreen <- washb_prescreen(Y = gamdat$Y,
                                                Ws = screenW, family = family, pval = pval, print = print))

    if(!is.null(forcedW)){
      Wscreen <- c(as.character(Wscreen), as.character(forcedW))
      #cat("\nNon-prescreened covariates: ", paste(forcedW, sep="", collapse=", "), "\n")
    }

    #drop perfectly multicollinear variables
    W <- subset(gamdat, select = Wscreen)
    # W$constant<-rep(1,nrow(gamdat))
    # tmp<-glm(constant ~ ., data=W, family=family)
    #https://daviddalpiaz.github.io/appliedstats/collinearity.html

    collinear_vars <- NULL
    Wdf <- W
    Wdf$constant<-rep(1,nrow(gamdat))
    for(i in 1:ncol(W)){
      tmp<-glm(constant ~ ., data=Wdf, family=family)
      todrop <-  NULL
      todrop <- suppressWarnings(names(tmp$coefficients)[-1][as.vector(vif(tmp)) > 10][1])
      if(!is.null(todrop)&!is.na(todrop)){
        collinear_vars <- c(collinear_vars,todrop)
        Wdf <- Wdf[,colnames(Wdf)!=todrop]
      }
    }

    # for(i in 1:length(colnames(W))){
    #   if(sum(grepl(colnames(W)[i], todrop))>0){
    #     to_drop <- c(to_drop, colnames(W)[i])
    #   }
    # }
    #W <- subset(W, select = -c(constant))
    # to_keep<-tmp$coefficients[!is.na(tmp$coefficients)]
    # to_keep<-names(to_keep[-which(names(to_keep) == "(Intercept)")])
    to_keep <- colnames(W)[!(colnames(W) %in% collinear_vars)]
    if(length(to_keep)!=length(colnames(W))){
      cat("\nDropped for collinearity with other covariates:\n",colnames(W)[!(colnames(W) %in% to_keep)])
    }
    W_processed <- W[which(colnames(W) %in% to_keep)]

    Wscreen <- colnames(W_processed)

    cat("\n\nCovariated included in model:\n",Wscreen)

  }else{
    Wscreen = NULL
  }

  if(!is.null(Wscreen)){
    d <- subset(gamdat, select = c("Y","X","id", "V", Wscreen))
  }else{
    d <- subset(gamdat, select = c("Y","X","id", "V"))
  }



  fullrows <- nrow(d)
  d <- d %>% filter(!is.na(Y))
  Yrows <- nrow(d)
  cat("\nRows dropped due to missing outcome: ",fullrows - Yrows,"\n")
  d <- d %>% filter(!is.na(X))
  Xrows <- nrow(d)
  cat("Rows dropped due to missing exposure: ",Yrows -Xrows,"\n")
  if(!is.null(W) & length(Wscreen)>0){
    cat("Percent missingness by covariate:\n")
    print(sapply(d[,-c(1:3)], function(x) round(sum(is.na(x))/nrow(X)*100,1)))
    d <- d[complete.cases(d),]
    cat("\nRows dropped due to missing covariates: ",Xrows - nrow(d),"\n")
  }
  cat("Final sample size: ", nrow(d),"\n")


  d$dummy<-1

  if(!is.null(W) & length(Wscreen)>0){

    #Make formula for adjusted model
    Ws <- subset(gamdat, select = c(Wscreen))

    #seperate factors and numeric
    W_factors <- colnames(Ws)[(grepl("factor", sapply(Ws, class))|grepl("character", sapply(Ws, class)))]
    W_numeric <- colnames(Ws)[(grepl("integer", sapply(Ws, class))|grepl("numeric", sapply(Ws, class)))]

    #seperate numeric indicators/few levels from continious
    indicator_vec <- rep(TRUE, length(W_numeric))
    for(i in 1:length(W_numeric)){
      N_unique <- length(unique(Ws[,W_numeric[i]]))
      if(N_unique>20){
        indicator_vec[i] <- FALSE
      }
    }

    W_indicator <- W_numeric[indicator_vec]
    W_continious <- W_numeric[!indicator_vec]

    #Create GAM equation
    if(length(W_continious)>0){
      eq_num <- paste0("s(", W_continious, ", bs=\"cr\")", collapse=" + ")
    }else{
      eq_num=NULL
    }
    if(length(W_factors)+length(W_indicator)>0){
      eq_fact <- paste0(" + ",paste0(c(W_factors,W_indicator), collapse=" + "))
    }else{
      eq_fact=NULL
    }

    #---------------------------------
    #fit model
    #---------------------------------

    #Check if X is binary or continious
    if(length(unique(d$X))>2){
      if(!is.null(V)){
        form <- paste0("Y~s(X, bs=\"cr\")+ V + X*V+",eq_fact," +",eq_num,"+ s(id,bs=\"re\",by=dummy)")
        form <- gsub("+ +","+",form, fixed=TRUE)
        equation <- as.formula(form)

        form_null <- paste0("Y~s(X, bs=\"cr\")+ V + ",eq_fact," +",eq_num,"+ s(id,bs=\"re\",by=dummy)")
        form_null <- gsub("+ +","+",form_null, fixed=TRUE)
        equation_null <- as.formula(form_null)
      }else{
        form <- paste0("Y~s(X, bs=\"cr\")+",eq_fact," +",eq_num,"+ s(id,bs=\"re\",by=dummy)")
        form <- gsub("+ +","+",form, fixed=TRUE)
        equation <- as.formula(form)
      }
    }else{
      if(!is.null(V)){
        form <- paste0("Y~X+ V + X*V +",eq_fact," +",eq_num,"+ s(id,bs=\"re\",by=dummy)")
        form <- gsub("+ +","+",form, fixed=TRUE)
        equation <- as.formula(form)

        form_null <- paste0("Y~Xc+ V + ",eq_fact," +",eq_num,"+ s(id,bs=\"re\",by=dummy)")
        form_null <- gsub("+ +","+",form_null, fixed=TRUE)
        equation_null <- as.formula(form_null)
      }else{
        form <- paste0("Y~X+",eq_fact," +",eq_num,"+ s(id,bs=\"re\",by=dummy)")
        form <- gsub("+ +","+",form, fixed=TRUE)
        equation <- as.formula(form)
      }
    }
    fit <- mgcv::gam(formula = equation,data=d)
    if(!is.null(V)){
      fit_null <- mgcv::gam(formula = equation_null,data=d)
      LRp <- lrtest(fit, fit_null)[2, 5]
    }
  }else{

    if(length(unique(d$X))>2){
      if(!is.null(V)){
        equation <- as.formula(paste0("Y~s(X, bs=\"cr\")+ V + X*V + s(id,bs=\"re\",by=dummy)"))
        fit <- mgcv::gam(formula = equation,data=d)

        equation_null <- as.formula(paste0("Y~s(X, bs=\"cr\")+ V + s(id,bs=\"re\",by=dummy)"))
        fit_null <- mgcv::gam(formula = equation_null,data=d)
        LRp <- lrtest(fit, fit_null)[2, 5]

      }else{
        fit <- mgcv::gam(Y~s(X, bs="cr")+s(id,bs="re",by=dummy),data=d)
      }
    }else{
      if(!is.null(V)){
        equation <- as.formula(paste0("Y~X + V + X*V + s(id,bs=\"re\",by=dummy)"))
        fit <- mgcv::gam(formula = equation,data=d)

        equation_null <- as.formula(paste0("Y~X + V + s(id,bs=\"re\",by=dummy)"))
        fit_null <- mgcv::gam(formula = equation_null,data=d)

        LRp <- lrtest(fit, fit_null)[2, 5]


      }else{
        fit <- mgcv::gam(Y~X +s(id,bs="re",by=dummy),data=d)
      }
    }
  }

  if(!is.null(V)){
    cat("\nInteraction p-value: ",LRp,"\n")
    return(list(fit=fit, dat=d, int.p=LRp))
  }else{
    return(list(fit=fit, dat=d))
  }
}

