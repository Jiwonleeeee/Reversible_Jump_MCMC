rm(list=ls())
load("/Users/wonny/Downloads/result0107.rda")
dataset <- read.csv("/Users/wonny/Downloads/actigraph_csv.csv")

subject1 <- subset(dataset, Index==1)
data <- subject1$axis1

# data <- (subject1$axis1+subject1$axis2+subject1$axis3)/3


len <- length(data)/10
newdata <- numeric(len)

for(i in 1:len){
  newdata[i] <- mean(data[(10*i-9):(10*i)])
}

y <- round(newdata)
t <- 1:length(y)
library(mvtnorm)

sigma2 <- mean(sigma2_store[2000:20000])
l <- mean(l_store[2000:20000])
#l <- 2
mu <- mean(mu_store[2000:20000])
g <- 10^-4

kernel_ftn <- function(t_input, l_input, g_input){
  K <- exp(-(outer(t_input,t_input,FUN = "-")^2)/l_input)+g_input*diag(length(t_input))
  return(K)
}


L <- length(t)

X <- t[t_post==1]
y <- y[t_post==1]

XX <- matrix(seq(0,L,length=1000), ncol=1)
SXX <- sigma2 * kernel_ftn(as.numeric(XX), l, g)
Sigma <- sigma2 * kernel_ftn(as.numeric(X), l, g)
SX <- sigma2 * (exp(-outer(as.numeric(XX), X, FUN = "-")^2/l))

Si <- solve(Sigma)
mup <- mu + SX %*% Si %*% (y-mu)
Sigmap <- SXX - SX %*% Si %*% t(SX)

S <- svd(Sigmap)
Sigmap <- S$u %*% diag(S$d) %*% t(S$u) # enforce the symmetry

YY <- rmvnorm(100, mup, Sigmap)

q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))


matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y")
points(X, y, pch=20, cex=0.5, col="blue")
lines(XX, mup, lwd=1)
lines(XX, q1, lwd=1, lty=2, col=2)
lines(XX, q2, lwd=1, lty=2, col=2)


matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y", xlim=c(0,100))
points(XX,rep(10,length(XX)))
points(X, y, pch=20, cex=1, col="blue")
lines(XX, mup, lwd=2)
lines(XX, q1, lwd=2, lty=2, col="red")
lines(XX, q2, lwd=2, lty=2, col="red")


lines(XX,sqrt(diag(SXX)), col=5, lwd=2)
lines(XX,sqrt(diag(Si)), col=6, lwd=1)
lines(XX,sqrt(diag(SX)), col=7, lwd=1)

lines(XX,sqrt(diag(Sigmap)), col="green", lwd=2)




matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y", xlim=c(400,600))
points(X, y, pch=20, cex=0.5, col="blue")
lines(XX,sqrt(diag(Sigmap)), col="green", lwd=2)
abline(h=sqrt(sigma2))
lines(XX, mup, lwd=1)
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)

