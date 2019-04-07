## R functions a lot like a calculator! ##

## You can write simple mathematical expressions for analysis:
1+1
abs(-100)
sqrt(25)

## You can save these results to specific variables using <- :
a <- 1+1
b <- abs(-100)
c <- sqrt(25)

## You can create a vector of values in the same way by using c():
a <- c(1,2,3,4,5,6,7,8,9,10)
b <- c(10,9,8,7,6,5,4,3,2,1)

## You can combine these vectors in to a matrix:
dat <- cbind(a, b)
dat

## You can plot these values using plot():
plot(dat[,1], dat[,2], xlab="Toxicity", ylab="Abundance")

## You can assess a simple linear model using lm():
mod <- lm(dat[,1] ~ dat[,2])
summary(mod)

## And you can add this line to the plot you've already made:
plot(dat[,1], dat[,2], xlab="Toxicity", ylab="Abundance")
abline(mod)

## You can install a package with functions you might need to use later:
install.packages("compute.es")

## You load these packages into each new R-session using library:
library(compute.es)