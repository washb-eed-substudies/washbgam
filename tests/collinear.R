
library(corrplot)
x <- seq(0, 100, 1)
# colinear with x
y <- x + 2.3
# almost colinear with x / some small gaussian noise
z <- x + rnorm(mean = 0, sd = 5, n = 101)
# uncorrrelated gaussian
w <- rnorm(mean = 0, sd = 1, n = 101)

# this frame is made to exemplify the procedure
df <- data.frame(x = x, y = y, z = z, w = w)

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

tmp<-glm(constant ~ ., data=df, family=family)
#https://daviddalpiaz.github.io/appliedstats/collinearity.html
todrop <- suppressWarnings(names(tmp$coefficients)[vif(tmp) > 5])

