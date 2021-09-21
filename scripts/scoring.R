# scoringRules

library(scoringRules)
library(raster)

pdf("~/NAM-Model-Validation/pdf/scoringCRPS.pdf")

# Plot some CRPSs for standard Gaussian distribution with truth at 0
# Then vary the scale from 1 to c(0,2,10)
x <- seq(-5,5,.1)
plot(x,crps_norm(x), type="l", main = "Negative CRPS for Gaussian, truth = 0")
lines(x, crps_norm(x, scale=0), col="blue")
lines(x, crps_norm(x, scale=2), col="orange")
lines(x, crps_norm(x, scale=10), col="red")
legend("topright", legend=paste("sd =",c(0,1,2,10)), lty=1, 
       col=c("blue","black","orange","red"))

# check to make sure my CRPS is the same
my_crps <- function(x, mu=0, sig=1){
  sig * (1/sqrt(pi) - 2*dnorm((x-mu)/sig)-(x-mu)/sig*(2*pnorm((x-mu)/sig)-1))
}

# it is typically used in negative orientation
all.equal(-my_crps(x), crps_norm(x))

# pkg is more robust than mine, of course:
all.equal(-my_crps(x, sig=0), crps_norm(x, scale=0))

# nonparametric version of CRPS
# based on equation 2 of Jordan etal 2018; scoringRules paper
my_np_crps <- function(y, scen){
  m <- length(scen)
  summin <- matrix(NA, m, m)
  for (i in 1:m) {
    summin[i,] <- abs(scen[i] - scen)
  }
  mean(abs(scen - y)) - 1/(2*m^2)*sum(summin)
}

for (s in 1:6) {
  print(s)
  load(paste0("~/NAM-Model-Validation/RData/prediction",s))
  dim(simvals)
  
  # scenarios is simvals combined with NAM
  scenarios <- simvals+NAM_pred$value
  scenarios <-ifelse(scenarios > 0, scenarios, 0)
  means <- rowMeans(scenarios)
  sds <- apply(scenarios, 1, sd)
  
  
  orig_NAM_crps <- crps_norm(y = ST4_pred$value, mean = NAM_pred$value, scale = 0)
  UQ_NAM_crps <- crps_norm(y = ST4_pred$value, mean = means, sd = sds)
  # plot(orig_NAM_crps,type="l")
  # plot(UQ_NAM_crps, type="l")
  
  # CRPS is exactly the absolute error of the forecast when deterministic prediction
  zlim <- c(0, max(orig_NAM_crps, UQ_NAM_crps))
  par(mfrow=c(2,2))
  plot(ST4_r - NAM_r, main=paste0("storm #", s, " error field"))
  plot(abs(ST4_r - NAM_r), main="abs error field\n CRPS for determ NAM", zlim = zlim)
  plot(rasterFromXYZ(cbind(coords, UQ_NAM_crps)), main="CRPS for our UQ", zlim = zlim)
  plot(rasterFromXYZ(cbind(coords, orig_NAM_crps - UQ_NAM_crps)), 
       main="diff of CRPSs\n positive is we win")
  plot(rasterFromXYZ(cbind(coords, ifelse(orig_NAM_crps - UQ_NAM_crps > 0, 1, -1))), 
       main="diff of CRPSs\n det vs Gaussian")
  
  
  plot(orig_NAM_crps, UQ_NAM_crps, main = "compare CRPSs: det vs UQ")
  abline(0,1, col="blue")
  
  np_CRPS <- c()
  for (i in 1:length(ST4_pred$value)) {
    # if(i %% 1000 == 0) print(i)
    np_CRPS[i] <- my_np_crps(ST4_pred$value[i], scenarios[i,])
  } 
  
  plot(UQ_NAM_crps, np_CRPS, main = "UQ: normal vs nonparametric", pch=".")
  abline(0,1,col="red")
  plot(rasterFromXYZ(cbind(coords, ifelse(orig_NAM_crps - np_CRPS > 0, 1, -1))), 
       main="diff of CRPSs\n det vs nonparam")
  
  # assumption of a normal predictive density
  mean(UQ_NAM_crps < orig_NAM_crps)
  mean(UQ_NAM_crps[ST4_pred$value > 1] < orig_NAM_crps[ST4_pred$value > 1])
  mean(UQ_NAM_crps[ST4_pred$value > 5] < orig_NAM_crps[ST4_pred$value > 5])
  mean(UQ_NAM_crps[ST4_pred$value > 10] < orig_NAM_crps[ST4_pred$value > 10])
  
  # nonparametric, using an empirical CDF
  mean(np_CRPS < orig_NAM_crps)
  mean(np_CRPS[ST4_pred$value > 1] < orig_NAM_crps[ST4_pred$value > 1])
  mean(np_CRPS[ST4_pred$value > 5] < orig_NAM_crps[ST4_pred$value > 5])
  mean(np_CRPS[ST4_pred$value > 10] < orig_NAM_crps[ST4_pred$value > 10])
  
  # do stuff on the original scale (square it)
  
  means2 <- rowMeans(scenarios^2)
  sds2 <- apply(scenarios^2, 1, sd)
  
  orig_NAM2_crps <- crps_norm(y = ST4_pred$value^2, mean = NAM_pred$value^2, scale = 0)
  UQ_NAM2_crps <- crps_norm(y = ST4_pred$value^2, mean = means2, sd = sds2)
  zlim2 = c(0, max(orig_NAM2_crps, UQ_NAM2_crps))
  
  par(mfrow=c(2,2))
  plot(ST4_r^2 - NAM_r^2, main = paste0("storm #", s, " orig scale EF"))
  plot(rasterFromXYZ(cbind(coords, orig_NAM2_crps)), main = "orig scale abs error")
  plot(rasterFromXYZ(cbind(coords, UQ_NAM2_crps)), main="CRPS for our UQ", zlim = zlim2)
  plot(rasterFromXYZ(cbind(coords, orig_NAM2_crps - UQ_NAM2_crps)), 
       main="diff of CRPSs\n positive is we win")
  plot(rasterFromXYZ(cbind(coords, ifelse(orig_NAM2_crps - UQ_NAM2_crps > 0, 1, -1))), 
       main="diff of CRPSs\n det vs Gaussian")
  
  
  np_CRPS2 <- c()
  for (i in 1:length(ST4_pred$value)) {
    # if(i %% 1000 == 0) print(i)
    np_CRPS2[i] <- my_np_crps(ST4_pred$value[i]^2, scenarios[i,]^2)
  } 
  
  plot(UQ_NAM2_crps, np_CRPS2, main = "UQ: normal vs nonparametric", pch=".")
  abline(0,1,col="red")
  plot(rasterFromXYZ(cbind(coords, orig_NAM2_crps - np_CRPS2)), 
       main="diff of CRPSs\n positive is we win")
  plot(rasterFromXYZ(cbind(coords, ifelse(orig_NAM2_crps - np_CRPS2 > 0, 1, -1))), 
       main="diff of CRPSs\n det vs nonparam")
  
  # assumption of a normal predictive density
  print("% of Gaussian UQ better than orig NAM: 0, 1, 25 mm")
  print(mean(UQ_NAM2_crps < orig_NAM2_crps))
  print(mean(UQ_NAM2_crps[ST4_pred$value > 1] < orig_NAM2_crps[ST4_pred$value > 1]))
  print(mean(UQ_NAM2_crps[ST4_pred$value > 5] < orig_NAM2_crps[ST4_pred$value > 5]))

  # nonparametric, using an empirical CDF
  print("% of nonparametric UQ better than orig NAM: 0, 1, 25 mm")
  print(mean(np_CRPS2 < orig_NAM2_crps))
  print(mean(np_CRPS2[ST4_pred$value > 1] < orig_NAM2_crps[ST4_pred$value > 1]))
  print(mean(np_CRPS2[ST4_pred$value > 5] < orig_NAM2_crps[ST4_pred$value > 5]))

  # note the change in the 0+ bar
  hist(orig_NAM2_crps - UQ_NAM2_crps, main = "CRPS: determ - normal UQ")
  hist(orig_NAM2_crps - np_CRPS2, main = "CRPS: determ - nonparam UQ")
  
}

dev.off()
