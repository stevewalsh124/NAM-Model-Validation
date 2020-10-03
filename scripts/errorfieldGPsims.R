# Simulating 2004-2017 precip error fields based on actual landfall locations
# Steve Walsh August 2020; GP simulation based on Marco's code from Spatial stats course

library(raster)
library(fields)
library(LaplacesDemon)

seed <- 1
set.seed(seed)

real.dom <- F
lo <- 80 #length.out
lb <- -15
ub <-  15

# Make directories for simulated fields and output, if they do not exist
seed.dir <- paste0("~/NAM-Model-Validation/csv/sim_df/seed",seed)
seed.dir.out <- paste0("~/NAM-Model-Validation/csv/simsMLEout/seed",seed)
seed.RDa.out <- paste0("~/NAM-Model-Validation/RData/RDatasim0/seed",seed)
if(!dir.exists(seed.dir)){dir.create(seed.dir)}
if(!dir.exists(seed.dir.out)){dir.create(seed.dir.out)}
if(!dir.exists(seed.RDa.out)){dir.create(seed.RDa.out)}

# Functions to make a matrix symmetric (pos def) and convert a vector to square matrix
symmtrz <- function(X){as.matrix((X+t(X))/2)}
as.square <- function(mat){matrix(mat, nrow=sqrt(length(mat)),ncol=sqrt(length(mat)))}

# Load in the true generating B and Sigma_theta matrices (estimated from true data)
load("~/NAM-Model-Validation/RData/truthB")
load("~/NAM-Model-Validation/RData/truthSig")

# Load in the landfall locations and intensities to use for Bx_i mean structure
loc_int <- read.csv("~/NAM-Model-Validation/csv/storm_levels_and_locs.csv", row.names = 1)

# If I wanted to have 47*n simulated storms total instead of 47
n <- 1
if(n != 1){loc_int <-  do.call("rbind", replicate(n, loc_int, simplify = FALSE))}

# Number of simulated storms (N) and parameters per storm to be estimated (P)
N <- 47 * n; P <- 3

# Design matrix based on landfall locations
xi <- with(loc_int, model.matrix(~ storm_locs))

# Generate true theta values based on loaded B and Sigma_theta
theta_i_sim <- matrix(nrow=N, ncol=P)
for (i in 1:N) {theta_i_sim[i,] <-  rmvn(n=1, mu=t(truthB%*%xi[i,]), Sigma = symmtrz(truthSig))}
# save(theta_i_sim, file=paste0("~/NAM-Model-Validation/RData/RDatasim0/theta_i_sim",seed))

# Look at histogram with overlays attempting to correspond to the mixture over locations (Bx_i)
xseq <- seq(-2,2,by=0.01)
par(mfrow=c(1,3))
for(i in 1:3){
  hist(theta_i_sim[,i], freq = F)
  seq1 <- dnorm(xseq, mean=truthB[i,]%*%c(1,0,0), sd=sqrt(truthSig[1,1]))*sum((loc_int$storm_locs=="ATL"))/nrow(loc_int)
  seq2 <- dnorm(xseq, mean=truthB[i,]%*%c(1,1,0), sd=sqrt(truthSig[2,2]))*sum((loc_int$storm_locs=="FL"))/nrow(loc_int)
  seq3 <- dnorm(xseq, mean=truthB[i,]%*%c(1,0,1), sd=sqrt(truthSig[3,3]))*sum((loc_int$storm_locs=="GULF"))/nrow(loc_int)
  lines(xseq, seq1+seq2+seq3)
}

# s <- 11
#
# args <- commandArgs(TRUE)
# if(length(args) > 0)
#   for(i in 1:length(args))
#     eval(parse(text=args[[i]]))

par(mfrow=c(2,2))
for (st in 1:N) {
  print(paste("this is storm", st))
  storm <- read.csv(list.files("~/NAM-Model-Validation/csv/error_df/subtractPWmeanF/", full.names = T)[st], row.names = 1)
  
  if(real.dom){
    # Simulated locations based on actual locations
    s <- storm[,c(2,3)]
    # plot(s)
  } else {
    # Simulated locations from unit square
    xaxis <- seq(lb, ub, length.out = lo)#storm$x#seq(0,1,0.01)
    yaxis <- seq(lb, ub, length.out = lo)#storm$y#seq(0,1,0.01)
    s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
    for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
  }

  t <- as.matrix(dist(s))
  # image(t)

  ##############################
  # Matern covariance function #
  ##############################

  library(geoR)

  # With nugget effect
  tau2 <- 0
  sigma2 <- exp(theta_i_sim[st,1])
  phi <- exp(theta_i_sim[st,2])
  nu <- exp(theta_i_sim[st,3])
  # plot(seq(0,1,0.001),c(tau2,rep(0,length(seq(0,1,0.001))-1))+sigma2*matern(seq(0,1,0.001),phi=(1/phi),kappa=nu))
  covmatrix <- diag(tau2,nrow(t)) + sigma2*matern(t,phi=phi,kappa=nu)
  # image(covmatrix)
  z1 <- t(chol(covmatrix)) %*% rnorm(nrow(s))
  # z <- matrix(NA,nrow=length(xaxis),ncol=length(yaxis))
  # for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) z[i,j] <- x[i+(j-1)*length(xaxis)]
  # persp(x=storm$x, y=storm$y, z=z)
  # image(z)
  sim.out <- cbind(x=s[,1],y=s[,2], z=z1)
  colnames(sim.out) <- c("x","y","z")
  write.csv(sim.out,
            file = paste0(seed.dir,"/","sim",if(st<10){"0"},st,".csv"), row.names = F)

  # plot(rasterFromXYZ(cbind(s, z1)), main = st); US(add=T)

  # plot(rasterFromXYZ(storm[,c(2,3,1)]), main = loc_int$storm_locs[st]); US(add=T)
}
