#Simualate error fields for each hurricane
library(SpatialEpi) #latlong2grid
library(dplyr)
getwd()

#check memory
# as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE))

params_deg_csv <- read.csv("csv/svg.param.ests.error_deg.csv")
params_csv <- read.csv("csv/svg.param.ests.error4.1.19.csv")
nam_df_csv <- read.csv("csv/namdf_nate2017_500km.csv")
nam_df_full <- read.csv("csv/namdf_nate2017.csv")
storms1 <- read.csv("csv/storms1.csv")[,2]
storms2 <- read.csv("csv/storms2.csv")[,2]
storms3 <- read.csv("csv/storms3.csv")[,2]
storms4 <- read.csv("csv/storms4.csv")[,2]
storms5 <- read.csv("csv/storms5.csv")[,2]
stormsTS <- read.csv("csv/stormsTS.csv")[,2]
stormsATL <- read.csv("csv/stormsATL.csv")[,2]
stormsGULF <- read.csv("csv/stormsGULF.csv")[,2]
stormsFL <- read.csv("csv/stormsFL.csv")[,2]
head(params_csv)
phivec  <- params_deg_csv[,2]
prRangevec <- params_deg_csv[,3]
tau2vec <- params_deg_csv[,4]
sig2vec <- params_deg_csv[,5]
summary(nam_df_csv)
xaxis <- nam_df_csv[,3]
yaxis <- nam_df_csv[,4]
xaxis_full <- nam_df_full[,3]
yaxis_full <- nam_df_full[,4]
nam_df_sort <- nam_df_csv[order(nam_df_csv$x, nam_df_csv$y),]
# s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
# for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
s <- nam_df_csv[,3:4]
s_full <- nam_df_full[,3:4]
plot(s)
plot(s_full)

##############################################
# Subsetting so that the area is a rectangle #
##############################################
xmin.loc <- -91
xmax.loc <- -81
ymin.loc <- 31
ymax.loc <- 40
s <- filter(s, x > xmin.loc & x < xmax.loc & y < ymax.loc & y > ymin.loc)
xaxis <- sort(unique(xaxis[xaxis > xmin.loc & xaxis < xmax.loc]))
yaxis <- sort(unique(yaxis[yaxis >  ymin.loc & yaxis <  ymax.loc]))

s_full <- filter(s_full, x > xmin.loc & x < xmax.loc & y < ymax.loc & y > ymin.loc)
xaxis_full <- sort(unique(xaxis_full[xaxis_full > xmin.loc & xaxis_full < xmax.loc]))
yaxis_full <- sort(unique(yaxis_full[yaxis_full >  ymin.loc & yaxis_full <  ymax.loc]))
plot(s_full)
s <- s_full
yaxis <- yaxis_full
xaxis <- xaxis_full
# xaxis <- sort(unique(xaxis))
# yaxis <- sort(unique(yaxis))
# plot(s)

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

#NAM_plotter loaded from PlotAllLogStormsHTML_Mac_copy.Rmd#############
plot(NAM_plotter, xlim= c(xmin.loc,xmax.loc), ylim=c(ymin.loc,ymax.loc),
     main="Original Forecast")

covmatrix <- diag(tau2,nrow(t)) + covfunc.Gaussian(t,phi,sigma2)
# image(covmatrix)
iters <- 100
forecast_plus_error_rasters <- list()
for (k in 1:iters) {
  print(k)
  x <- matrix.sqrt(covmatrix) %*% rnorm(nrow(s))
  z <- matrix(NA,nrow=length(xaxis),ncol=length(yaxis))
  for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) z[i,j] <- x[i+(j-1)*length(xaxis)]
  # par(mfrow=c(1,2))
  # image(z, col=terrain.colors(100))
  # image(-z, col=terrain.colors(100), main="image(-z, col=terrain.colors(100))")
  
  # ; dim(z)
  # persp(x=xaxis, y=yaxis, z=z, phi = 40)
  # plot(rasterFromXYZ(cbind(s,as.vector(z))), main="rasterFromXYZ(cbind(s,as.vector(z)))")
  forecast_plus_error_rasters[[k]] <- NAM_plotter + rasterFromXYZ(cbind(s,as.vector(z)))
  plot(NAM_plotter + rasterFromXYZ(cbind(s,as.vector(z))),
       main="Forecast + random error field")
  if(k==1){ total_them_up <-forecast_plus_error_rasters[[k]]}
  if(k>=2){ total_them_up <-total_them_up + forecast_plus_error_rasters[[k]]}
}

par(mfrow=c(2,3))
plot(NAM_plotter, main="NAM", xlim=c(xmin.loc,xmax.loc), ylim=c(ymin.loc,ymax.loc))
plot(forecast_plus_error_rasters[[1]], main="NAM + a random error field")
plot(forecast_plus_error_rasters[[1]]-NAM_plotter, main="A random error field")
plot(total_them_up/k,main="Avg of 100 Random EFs + NAM")
plot(NAM_plotter-total_them_up/k,main="Thus, avg of 100 Random EFs")
plot(k*NAM_plotter-total_them_up,main="100 Random EFs")

# total_them_up <- raster()
# total_them_up <-forecast_plus_error_rasters[[1]]

# surface2 <- interp(x=s[,1], y=s[,2], z=z, duplicate = "mean")
# persp(surface2)


#somethings wrong with z2 (row/column mismatch?)
z2 <- matrix(NA,nrow=length(unique(nam_df_csv$x)),ncol=length(unique(nam_df_csv$y)))
for (i in 1:length(unique(nam_df_csv$x))) for (j in 1:length(unique(nam_df_csv$y))){
  z2[i,j] <- x[i+(j-1)*length(unique(nam_df_csv$x))]
}
# image(sort(unique(nam_df_csv$x)),sort(unique(nam_df_csv$y)),z2)
