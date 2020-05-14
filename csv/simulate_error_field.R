#Simualate error fields for each hurricane
library(SpatialEpi) #latlong2grid

setwd("/nas/home/grad/walsh124/")

#check memory
as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE))

params_deg_csv <- read.csv("NAM-Model-Validation/svg.param.ests.error_deg.csv")
params_csv <- read.csv("NAM-Model-Validation/svg.param.ests.error4.1.19.csv")
nam_df_csv <- read.csv("NAM-Model-Validation/namdf_nate2017_250km.csv")
storms1 <- read.csv("NAM-Model-Validation/storms1.csv")[,2]
storms2 <- read.csv("NAM-Model-Validation/storms2.csv")[,2]
storms3 <- read.csv("NAM-Model-Validation/storms3.csv")[,2]
storms4 <- read.csv("NAM-Model-Validation/storms4.csv")[,2]
storms5 <- read.csv("NAM-Model-Validation/storms5.csv")[,2]
stormsTS <- read.csv("NAM-Model-Validation/stormsTS.csv")[,2]
stormsATL <- read.csv("NAM-Model-Validation/stormsATL.csv")[,2]
stormsGULF <- read.csv("NAM-Model-Validation/stormsGULF.csv")[,2]
stormsFL <- read.csv("NAM-Model-Validation/stormsFL.csv")[,2]
head(params_csv)
phivec  <- params_deg_csv[,2]
prRangevec <- params_deg_csv[,3]
tau2vec <- params_deg_csv[,4]
sig2vec <- params_deg_csv[,5]
summary(nam_df_csv)
# Getting error with persp() final command, this sorting didn't help 
# (both xaxis and yaxis have to strictly increase?)
nam_df_sort <- nam_df_csv[order(nam_df_csv$x, nam_df_csv$y),]
xaxis <- nam_df_sort[,3]
yaxis <- nam_df_sort[,4]

# s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
# for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
s <- nam_df_csv[,3:4]
plot(s)
t <- as.matrix(dist(s))
# image(t)

################################
# Gaussian covariance function #
################################

covfunc.Gaussian <- function(t,phi,sigma2) {sigma2 * exp(- phi^2 * t^2)}

matrix.sqrt <- function(H)
{
  # Computes square root of nonnegative definite symmetric matrix using spectral decomposition
  
  if(nrow(H)==1) {H.sqrt = matrix(sqrt(H),nrow=1,ncol=1)} else
  {
    H.eigen = eigen(H)
    H.eigen.values = H.eigen$values    
    H.eigen.values[abs(H.eigen$values) < 10^(-10)] = 0
    H.sqrt = H.eigen$vectors %*% diag(sqrt(H.eigen.values)) %*% t(H.eigen$vectors)
  }  
  
  H.sqrt
}

# With nugget effect

tau2 <-mean(tau2vec[stormsGULF]) #previous 0.5
sigma2 <- mean(sig2vec[stormsGULF]) #1
phi <- mean(phivec[stormsGULF]) #5
plot(seq(0,1,0.001),c(tau2,rep(0,length(seq(0,1,0.001))-1))+covfunc.Gaussian(seq(0,1,0.001),phi,sigma2))

covmatrix <- diag(tau2,nrow(t)) + covfunc.Gaussian(t,phi,sigma2)
# image(covmatrix)
x <- matrix.sqrt(covmatrix) %*% rnorm(nrow(s))
z <- matrix(NA,nrow=length(xaxis),ncol=length(yaxis))
for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) z[i,j] <- x[i+(j-1)*length(xaxis)]
image(z)
persp(x=xaxis, y=yaxis, z=z)
# persp(x=unique(xaxis), y=unique(yaxis), z=z)
