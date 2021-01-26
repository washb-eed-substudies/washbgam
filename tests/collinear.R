
library(corrplot)
x <- seq(0, 100, 1)
# colinear with x
y <- x + 2.3
# almost colinear with x / some small gaussian noise
z <- x + rnorm(mean = 0, sd = 5, n = 101)
# uncorrrelated gaussian
w <- rnorm(mean = 0, sd = 1, n = 101)

# this frame is made to exemplify the procedure
df <- data.frame(x = x, y = y, z = z, w = w, clusterid=1:length(x))

library(washb)
library(washbgam)
res <- fit_RE_gam(d=df, Y="y", X="x", W=c("y","x","z","w"), forcedW = NULL, family="gaussian")


d=df
Y="y"
X="x"
W=c("y","x","z","w")
forcedW = NULL
family="gaussian"
V=NULL
id="clusterid"
pval = 0.2
print=TRUE

corr.matrix <- cor(df)
corrplot.mixed(corr.matrix)

library(corrplot)
corrplot(cor(df))
library(caret)
indexesToDrop <- findCorrelation(cor(df), cutoff = 0.8)
corrplot(cor(df[,-indexesToDrop]))



# #drop perfectly multicollinear variables
# constant<-rep(1,nrow(df))
# tmp<-lm(constant ~ ., data=df)
# to_keep<-tmp$coefficients[!is.na(tmp$coefficients)]
# to_keep<-names(to_keep[-which(names(to_keep) == "(Intercept)")])
# df_result<-df[to_keep]

#Use vif
library(faraway)
tmp<-lm(constant ~ ., data=df)
#https://daviddalpiaz.github.io/appliedstats/collinearity.html
vif(tmp)
todrop <- suppressWarnings(names(tmp$coefficients)[as.vector(vif(tmp)) > 10])



#
library(FSelector)
linear.correlation(constant ~ ., data=df)



#
x1 = runif(1000)
x2 = runif(1000)
x3 = x1 + x2
x4 = runif(1000)
x5 = runif(1000)*0.00000001 +x4
x6 = x5 + x3
x = data.frame(x1, x2, x3, x4, x5, x6)
X <- as.matrix(x)
  qr.X <- qr(X, tol=1e-9, LAPACK = FALSE)
(rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
(keep <- qr.X$pivot[seq_len(rnkX)])
## 1 2 4 5
X2 <- X[,keep]



        df
        X <- as.matrix(df)
        qr.X <- qr(X, tol=1e-3, LAPACK = FALSE)
        (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
        (keep <- qr.X$pivot[seq_len(rnkX)])
        ## 1 2 4 5
        X2 <- X[,keep]
