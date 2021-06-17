###  Loading a list "barry_EFs" with 5 elements: coordinates, 1000 error fields, 
###  NAM forecast, and ST4 observation, pointwise mean
###  All precip is in log mm

# install.packages(c("raster","fields"))
library(raster)
library(fields)

## prediction storm to analyze (1-6)
s <- 5

# May need to change this file path
load(paste0("RData/prediction",s))

# coords  <- barry_EFs[[1]]
# simvals <- barry_EFs[[2]]
# NAM <- barry_EFs[[3]]
# ST4 <- barry_EFs[[4]]
# PWmean <- barry_EFs[[5]]

coords <- NAM_pred[,1:2]
NAM <- NAM_pred
ST4 <- ST4_pred[,c(3,4,2)]
PWmean <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_mean.grd")
  
par(mfrow=c(2,2))
plot(rasterFromXYZ(NAM), main = "NAM, bias adjusted"); US(add=T, col="gray80")
plot(PWmean, main = "Pointwise Mean map"); US(add=T, col="gray80")
plot(rasterFromXYZ(NAM) + PWmean, main = "NAM, no bias adjusted"); US(add=T, col="gray80")
plot(rasterFromXYZ(ST4), main = "ST4 observed precip"); US(add=T, col="gray80")

for (i in 1:4) {
  plot(rasterFromXYZ(cbind(coords, simvals[,i])), main = paste("Error field", i))
  US(add=T, col="gray80")
}

sqNAM <- cbind(NAM[1:2],(NAM[,3])^2)
plot(rasterFromXYZ(sqNAM), main = "NAM, orig scale (mm)")

sqST4 <- cbind(ST4[1:2],(ST4[,3])^2)
plot(rasterFromXYZ(sqST4), main = "ST4, orig scale (mm)")


# One way to visualize error fields on mm scale
# exp_EFs <- ifelse(simvals < 0, -exp(-simvals), exp(simvals))
exp_EFs <- ifelse(simvals < 0, -simvals^2, simvals^2)
for (i in 1:2) {plot(rasterFromXYZ(cbind(coords, exp_EFs[,i])), main= paste0("Sqrd Error field #", i))}

par(mfrow=c(2,2))
# # Bad
# for (i in 1:4) {
#   plot(rasterFromXYZ(sqNAM) + rasterFromXYZ(cbind(coords, exp_EFs[,i])),
#        main= paste0("exp NAM + exp Error field #", i))
# }

# Good
for (i in sample(1:ncol(simvals),4)) {
  plot(rasterFromXYZ(cbind(coords, (NAM[,3] + simvals[,i])^2)), 
       main= paste0("(NAM + Error field)^2 #", i))
}

## NAM w/ bias adj, NAM w/o bias adj, ST4 max precip values
max((NAM[,3])^2)
max(values((rasterFromXYZ(NAM)+PWmean)), na.rm = T)^2
max((ST4[,3])^2)
maxes <- apply((NAM[,3] + simvals)^2, 2, max)
par(mfrow=c(1,1))
hist(maxes, breaks=40, main = paste("Maxes for storm", s))
abline(v = max(ST4[,3])^2, col="blue", lwd=2)
abline(v = max(NAM[,3])^2, col="green", lwd=2)
legend(legend = c("NAM","ST4"), "topright", col=c("green","blue"), lty=1, lwd=2)
