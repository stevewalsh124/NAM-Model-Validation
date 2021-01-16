# By M.A.R. Ferreira, 2014. Updated by M.A.R. Ferreira, 2016.

xaxis <- seq(0,1,0.025)
yaxis <- seq(0,1,0.025)

s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
# plot(s)

t <- as.matrix(dist(s))
# image(t)

###################################
# Exponential covariance function #
###################################

covfunc.exponential <- function(t,phi,sigma2) {sigma2 * exp(-phi*t)}

# With nugget effect

tau2 <- 0
sigma2 <- 1
phi <- 1
plot(seq(0,1,0.001),c(tau2,rep(0,length(seq(0,1,0.001))-1))+covfunc.exponential(seq(0,1,0.001),phi,sigma2))

covmatrix <- diag(tau2,nrow(t)) + covfunc.exponential(t,phi,sigma2)
#image(covmatrix)
x <- t(chol(covmatrix)) %*% rnorm(nrow(s))
z <- matrix(x,nrow=length(xaxis),ncol=length(yaxis),byrow=FALSE)
image(z)
# persp(x=xaxis, y=yaxis, z=z)

########################
# Anisotropic Modeling #
########################

# C(s) = sigma2 * rho((s-s')^T %*% B %*% (s-s'))
# Banerjee et al pg 31, B must be positive definite

# Higdon Temperatures in the North Atlantic
# d transformed to r: r = diag(c(sigma1, sigma2)) %*% matrix(cos(theta), -sin(theta), sin(theta), cos(theta)) %*% t(d)

sigvec <- c(1,3)
theta <- pi/2

tt <- matrix(NA, nrow = nrow(s), ncol = nrow(s))

for (i in 1:nrow(s)) {
  for (j in i:nrow(s)) {
    d <- s[i,] - s[j,]
    r = diag(sigvec) %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)),2,2) %*% (d)
    tt[i,j] <- tt[j,i] <- sqrt(r[1]^2+r[2]^2)#tau2 + covfunc.exponential(abs(r[1]-r[2]),phi,sigma2)#abs(r[1]-r[2])
  }
}

covmatrix_anisop <- diag(tau2,nrow(tt)) + covfunc.exponential(tt,phi,sigma2)
#image(covmatrix)
x <- t(chol(covmatrix_anisop)) %*% rnorm(nrow(s))
z <- matrix(x,nrow=length(xaxis),ncol=length(yaxis),byrow=FALSE)
image(z)


## Ergodicity assumed in next comment:
## C(h) = (lim u -> Inf gamma(u)) - gamma(h), gamma is the semivariogram


# Nonstationary example

yield <- function(xi1, xi2)  
{
  xi1 <- 7*xi1+1     ## these two
  xi2 <- 900*xi2+100 ## allow seq(0,1) below
  xi1 <- 3*xi1 - 15
  xi2 <- xi2/50 - 13
  xi1 <- cos(0.5)*xi1 - sin(0.5)*xi2
  xi2 <- sin(0.5)*xi1 + cos(0.5)*xi2
  y <- 1.15^(-xi1^2/80 - 0.5*(xi2 + 0.03*xi1^2 - 40*0.03)^2)
  return(100*y)
}

xi1 <- xi2 <- seq(0,1,length.out = 100)
g <- expand.grid(xi1, xi2)
y <- yield(g[,1], g[,2]) + rnorm(10000,0,8)
# persp(xi1, xi2, matrix(y, ncol=length(xi2)), theta=45, phi=45, 
#       lwd=0.5, xlab="xi1 : time", ylab="xi2 : temperature", 
#       zlab="yield", expand=0.4)


cols <- heat.colors(128)
image(xi1, xi2, matrix(y, ncol=length(xi2)), col=cols, 
      xlab="xi1 : time", ylab="xi2 : temperature")
# contour(xi1, xi2, matrix(y, ncol=length(xi2)), nlevels=4, add=TRUE)
