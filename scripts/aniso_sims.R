# Steve Walsh, June 2021, anisotropic random fields with MLEs
# computed from geoR::likfit and convoSPAT::Aniso_fit

# blas_set_num_threads(1)

Nsims <- 30
seed <- 1
set <- 1
lo <- 80

write.pdf <- F

## truestart:    the aniso params only are true starting values
## alltruestart: aniso params as well as sigma2 and phi are at true starting values
## neither:      sigma2, phi, maj.min = 1; theta = 0
truestart <- T
alltruestart <- F
if(truestart & alltruestart){stop("truestart or alltruestart, not both")}

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

set.seed(seed)

box <- 20
if(set==12){box <- 200}
if(set==13){box <- 1}
xaxis <- seq(0,box,length.out = lo)
yaxis <- seq(0,box,length.out = lo)

s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
# plot(s)

t <- as.matrix(dist(s))
# image(t)

###################################
# Exponential covariance function #
###################################

covfunc.exponential <- function(t,phi,sigma2) {sigma2 * exp(-t/phi)}

# With nugget effect
if(set == 1){truths <- c(beta=0, tau2=0.01, sigma2=0.8, phi=0.4, theta=pi/4, maj.min=3)}
if(set == 2){truths <- c(beta=0, tau2=0.1, sigma2=6, phi=1, theta=pi/8, maj.min=2)}
if(set == 3){truths <- c(beta=0, tau2=0.2, sigma2=3, phi=2, theta=pi/3, maj.min=1.3)}
if(set == 4){truths <- c(beta=0, tau2=0.2, sigma2=3, phi=3, theta=pi/2, maj.min=2.5)}
if(set == 5){truths <- c(beta=0, tau2=0.05, sigma2=1, phi=4, theta=pi/4, maj.min=2.5)}
if(set == 6){truths <- c(beta=0, tau2=0.02, sigma2=4, phi=1, theta=pi/4, maj.min=2)}
if(set == 7){truths <- c(beta=0, tau2=0.04, sigma2=5, phi=2, theta=pi/8, maj.min=2)}
if(set == 8){truths <- c(beta=0, tau2=0.04, sigma2=5, phi=1.5, theta=pi/4, maj.min=1)}
if(set == 9){truths <- c(beta=0, tau2=0.03, sigma2=4.5, phi=2, theta=pi/3, maj.min=3)}
if(set == 10){truths <- c(beta=0, tau2=0.03, sigma2=4.5, phi=5, theta=pi/3, maj.min=3)}
if(set == 11){truths <- c(beta=0, tau2=0.03, sigma2=4.5, phi=0.9, theta=pi/3, maj.min=3)}
if(set == 12){truths <- c(beta=0, tau2=0.03, sigma2=4.5, phi=30, theta=pi/3, maj.min=3)}
if(set == 13){truths <- c(beta=0, tau2=0.03, sigma2=4.5, phi=0.2, theta=pi/3, maj.min=3)}

if(!exists("truths")){stop("True values have not been set")}

if(!dir.exists(paste0("~/NAM-Model-Validation/csv/aniso/set",set))){
  dir.create(paste0("~/NAM-Model-Validation/csv/aniso/set",set))}

if(write.pdf) write.csv(truths, file = paste0("~/NAM-Model-Validation/csv/aniso/set",set,"/set",set,"truths.csv"))

tau2 <- truths["tau2"]
sigma2 <- truths["sigma2"]
phi <- truths["phi"]
maj.min <- c(1, truths["maj.min"])
theta <- truths["theta"]

######################
# Isotropic Modeling #
######################

# plot(seq(0,1,0.001),c(tau2,rep(0,length(seq(0,1,0.001))-1))+covfunc.exponential(seq(0,1,0.001),phi,sigma2))
# 
# covmatrix <- diag(tau2,nrow(t)) + covfunc.exponential(t,phi,sigma2)
# #image(covmatrix)
# x <- t(chol(covmatrix)) %*% rnorm(nrow(s))
# z <- matrix(x,nrow=length(xaxis),ncol=length(yaxis),byrow=FALSE)
# image(z)
# # persp(x=xaxis, y=yaxis, z=z)

########################
# Anisotropic Modeling #
########################

# C(s) = sigma2 * rho((s-s')^T %*% B %*% (s-s'))
# Banerjee et al pg 31, B must be positive definite

# Higdon Temperatures in the North Atlantic (1998)
# d transformed to r: 
# r = diag(c(sigma1, sigma2)) %*% matrix(cos(theta), -sin(theta), sin(theta), cos(theta)) %*% t(d)

tt <- matrix(NA, nrow = nrow(s), ncol = nrow(s))

for (i in 1:nrow(s)) {
  for (j in i:nrow(s)) {
    d <- s[i,] - s[j,]
    r = diag(maj.min) %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)),2,2) %*% (d)
    tt[i,j] <- tt[j,i] <- sqrt(r[1]^2+r[2]^2)#tau2 + covfunc.exponential(abs(r[1]-r[2]),phi,sigma2)#abs(r[1]-r[2])
  }
}

covmatrix_anisop <- diag(tau2,nrow(tt)) + covfunc.exponential(tt,phi,sigma2)
#image(covmatrix)

## Ergodicity assumed in next comment:
## C(h) = (lim u -> Inf gamma(u)) - gamma(h), gamma is the semivariogram

# #########################
# # Nonstationary example #
# #########################
#
# yield <- function(xi1, xi2)  
# {
#   xi1 <- 7*xi1+1     ## these two
#   xi2 <- 900*xi2+100 ## allow seq(0,1) below
#   xi1 <- 3*xi1 - 15
#   xi2 <- xi2/50 - 13
#   xi1 <- cos(0.5)*xi1 - sin(0.5)*xi2
#   xi2 <- sin(0.5)*xi1 + cos(0.5)*xi2
#   y <- 1.15^(-xi1^2/80 - 0.5*(xi2 + 0.03*xi1^2 - 40*0.03)^2)
#   return(100*y)
# }
# 
# xi1 <- xi2 <- seq(0,1,length.out = 100)
# g <- expand.grid(xi1, xi2)
# y <- yield(g[,1], g[,2]) + rnorm(10000,0,8)
# # persp(xi1, xi2, matrix(y, ncol=length(xi2)), theta=45, phi=45, 
# #       lwd=0.5, xlab="xi1 : time", ylab="xi2 : temperature", 
# #       zlab="yield", expand=0.4)
# 
# 
# cols <- heat.colors(128)
# image(xi1, xi2, matrix(y, ncol=length(xi2)), col=cols, 
#       xlab="xi1 : time", ylab="xi2 : temperature")
# # contour(xi1, xi2, matrix(y, ncol=length(xi2)), nlevels=4, add=TRUE)


##################################
# MLEs for Anisotropy/Stationary #
##################################

myMLEs <- mySpats <- matrix(NA, Nsims, 6)
colnames(myMLEs) <- names(truths)
colnames(mySpats) <- c("beta", "tausq", "sigma2", "lam1", "lam2", "eta")

times <- times2 <- c()

library(geoR)
library(fields)
library(convoSPAT)

if(write.pdf) {pdf(paste0("~/NAM-Model-Validation/pdf/aniso/set",set,
                          "_",Nsims,"_",nrow(s),"_box",box,"_seed",seed,".pdf"))}

for (i in 1:Nsims) {
  print(paste0("sim #", i))
  
  x <- t(chol(covmatrix_anisop)) %*% rnorm(nrow(s))
  z <- matrix(x,nrow=length(xaxis),ncol=length(yaxis),byrow=FALSE)
  
  if(i < 10){
    image.plot(z, main=paste0("sigmasq",sigma2," phi",phi, 
                              " ratio",maj.min[2], " angle",round(180/pi*theta,3)))
  }

  tic <- proc.time()[3]
  if(truestart){
    myMLE <- likfit(as.geodata(cbind(s,x)), ini.cov.pars = c(1,1), fix.psiA = F, 
                    psiA = theta, fix.psiR = F, psiR = maj.min[2])} 
  if(alltruestart){
    myMLE <- likfit(as.geodata(cbind(s,x)), ini.cov.pars = c(sigma2,phi), fix.psiA = F, 
                    psiA = theta, fix.psiR = F, psiR = maj.min[2])}
  if(!truestart & !alltruestart){
    myMLE <- likfit(as.geodata(cbind(s,x)), ini.cov.pars = c(1,1), fix.psiA = F, 
                    psiA = 0, fix.psiR = F, psiR = 1)}
  
  myMLEs[i,] <-  c(myMLE$beta, myMLE$tausq, myMLE$sigmasq, myMLE$phi, myMLE$aniso.pars[1], myMLE$aniso.pars[2])
  toc <- proc.time()[3]
  print(toc - tic)
  times[i] <- toc - tic
  
  
  tic <- proc.time()[3]
  my_spat <- Aniso_fit(geodata = as.geodata(cbind(s,x)))
  mySpats[i,] <-  c(my_spat$beta.GLS, my_spat$tausq.est, my_spat$sigmasq.est, my_spat$aniso.pars)
  toc <- proc.time()[3]
  print(toc - tic)
  times2[i] <- toc - tic
}

if(write.pdf){
  write.csv(cbind(myMLEs, times),
            file = paste0("~/NAM-Model-Validation/csv/aniso/set",set,
                          "/aniso_sim_results_",Nsims,"_",nrow(x),"_box",box,"_seed",seed,
                          if(truestart){"_truestart"}, if(alltruestart){"_alltruestart"},".csv"))
  
  write.csv(cbind(mySpats, times2),
            file = paste0("~/NAM-Model-Validation/csv/aniso/set",set,
                          "/aniso_sim_results_",Nsims,"_",nrow(x),"_box",box,"_seed",seed,
                          "_spat.csv"))
}


par(mfrow=c(2,3))
for (i in 1:6) {
  hist(myMLEs[,i], main = paste(names(truths)[i], round(truths[i],3)))
  abline(v=truths[i], col="blue", lwd=2)
  if(i==6)   abline(v=1/truths[i], col="red", lwd=2, lty=3)
  if(i==5)   abline(v=pi/2 - truths[i], col="red", lwd=2, lty=3)
  
}

truths_spat <- c(truths[1:4], truths[6], truths[5])
spats <- cbind(mySpats[,1:3], sqrt(mySpats[,4]), sqrt(mySpats[,5]/mySpats[,4]), mySpats[,6])
par(mfrow=c(2,3))
for (i in 1:6) {
  hist(spats[,i], main = paste(names(truths_spat)[i], round(truths_spat[i],3)))
  if(i==6)   abline(v=pi/2 - truths_spat[i], col="red", lwd=2, lty=3)
  if(i==5)   abline(v=1/truths_spat[i], col="red", lwd=2, lty=3)
  abline(v=truths_spat[i], col="blue", lwd=2)
}

if(write.pdf) dev.off()
