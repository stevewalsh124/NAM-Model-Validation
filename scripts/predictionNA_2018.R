##  Prediction for the future storms (validate with 2018 and 2019 method)
##  Establish regressors for an incoming storm (eg: $\boldsymbol{x}_p^T = [1,1,0]$ for a Florida storm)
##  Use x_p to generate theta_i^*=Bx_p + omega_i for i in 1:1000
##    From these theta_i^*, generate 1000 corresponding error fields
##  Add these generated error fields to the available NAM forecast 
##  Check to see if pointwise 95% percentile is greater than $\approx$ 95\% of pixels for the corresponding ST4. 
##  Compute probability maps that rain totals will surpass a particular threshold, etc.

before <- Sys.time()

#Change NAM_pred, ST4_pred, and x_pred
#Change name of PDF

NAM_pred <- read.csv("~/NAM-Model-Validation/csvtest/prediction/2018michael/NAM_GIS.csv")
ST4_pred <- read.csv("~/NAM-Model-Validation/csvtest/prediction/2018michael/ST4_GIS.csv")
x_pred <- c(1,0,1) # Alberto 2018: GULF storm

Ngen <- 5000
subtractPWmean <- T
PWstamp <- ifelse(subtractPWmean, "subtractpw", "nopw")

# pdf(paste0("~/NAM-Model-Validation/pdf/prediction/prediction_Michael2018_GIS_GHiG_NA_",PWstamp, Ngen,".pdf"))
library(LaplacesDemon)
library(geoR)
library(raster)
library(RandomFields)

mask <- raster("~/NAM-Model-Validation/lsmask.nc")
mask[mask==-1]  <- NA
extent(mask)[1] <- extent(mask)[1]-360
extent(mask)[2] <- extent(mask)[2]-360
mask.regrid <- resample(mask, projectRaster(raster(
  "~/NAM-Model-Validation/nam_218_20050829_1200_f012.grib"),
  crs = "+proj=longlat +datum=WGS84"), method='ngb') 

PW_mean <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_mean.grd")*mask.regrid

# run GibbsSamplerHurrRegr first for B and as.square
B_pred <- as.square(apply(B[burn:iters,], 2, median))

Sigma_theta_pred <- as.square(apply(Sigma_theta[burn:iters,], 2, median))

theta_pred <- rmvn(Ngen, B_pred %*% x_pred, Sigma_theta_pred)

# (neg_rows <- which(theta_pred < 0))
# for (i in 1:length(neg_rows)) {
#   col2chg <- 1
#   count <- Ngen
#   # print(c(col2chg,count))
#   while(neg_rows[i] > count){
#     count <- count + Ngen
#     col2chg <- col2chg + 1
#     # print(c(count, col2chg))
#   }
#   theta_pred[neg_rows[i]] <- median(theta_pred[,col2chg])
# }
# (neg_rows <- which(theta_pred < 0))

NAM_r <- rasterFromXYZ(NAM_pred[,c(3,4,2)])
ST4_r <- rasterFromXYZ(ST4_pred[,c(3,4,2)]) #row order messed up from logprecipvariogram write.csv
par(mfrow=c(1,2))
plot(NAM_r, main = "NAM forecast")
plot(ST4_r, main = "ST4 forecast")

if(subtractPWmean){
  NAM_r <- NAM_r - PW_mean
  NAM_pred <- as.data.frame(rasterToPoints(NAM_r))
  colnames(NAM_pred) <- c("x","y","value")
}

coords <- cbind(NAM_pred$x, NAM_pred$y)
coordsST4 <- cbind(ST4_pred$x, ST4_pred$y)

simvals <- matrix(NA, nrow = nrow(NAM_pred), ncol = Ngen) 
for (g in 1:Ngen) {
  if(g %% 50 == 0) cat(g,"\n")
  simvals[,g] <- RFsimulate(model = RMwhittle(nu = exp(theta_pred[g,3]), var = exp(theta_pred[g,1]), 
                                              scale = exp(theta_pred[g,2])), 
                            err.model = RMnugget(var = unique(all_storm_res[,"MLEnugget"])),
                            x = cbind(NAM_pred$x, NAM_pred$y))$variable1
}

# example error field
sim <- cbind(coords, simvals[,1])
sim_r <-rasterFromXYZ(sim)

plot(NAM_r-ST4_r, main="Actual Error Field")
plot(sim_r, main="Simulated Error Field")

# calculate 95%ile for each point
ub_rain <- apply(simvals, 1, function(x) quantile(x, .95))
ub_rain99 <- apply(simvals, 1, function(x) quantile(x, .99))
ub_rain100 <- apply(simvals, 1, max)
plot(rasterFromXYZ(cbind(coords, ub_rain+NAM_pred$value)), main = "NAM + 95% UB on EFs")

# compare NAM with ST4, and NAM+UB with ST4
par(mfrow=c(1,3))
off_base <- ifelse(NAM_pred$value < ST4_pred$value, 1, 0)
mean(off_base)
plot(rasterFromXYZ(cbind(coords,off_base)), main= bquote("NAM < ST4"~.(round(mean(off_base),4))))

off_est <- ifelse(ub_rain+NAM_pred$value < ST4_pred$value, 1, 0)
mean(off_est)
plot(rasterFromXYZ(cbind(coords,off_est)), main=bquote("NAM + 95% UB < ST4"~.(round(mean(off_est),4))))

off_est99 <- ifelse(ub_rain99+NAM_pred$value < ST4_pred$value, 1, 0)
mean(off_est99)
plot(rasterFromXYZ(cbind(coords,off_est99)), main=bquote("NAM + 99% UB < ST4"~.(round(mean(off_est99),4))))

off_est100 <- ifelse(ub_rain100+NAM_pred$value < ST4_pred$value, 1, 0)
mean(off_est100)
plot(rasterFromXYZ(cbind(coords,off_est100)), main=bquote("NAM + 100% UB < ST4"~.(round(mean(off_est100),4))))


#2in = 50.8mm
sim_prob <- matrix(NA, nrow = nrow(simvals), ncol = Ngen)
for (i in 1:Ngen) {
  sim_prob[,i] <- ifelse(simvals[,i]+NAM_pred$value > log(50.8), 1, 0)
}

two_in_rain_f <- ifelse(NAM_pred$value > log(50.8), 1, 0)
two_in_rain_o <- ifelse(ST4_pred$value > log(50.8), 1, 0)
two_in_rain_oNA <- ifelse(ST4_pred$value > log(50.8), 1, NA)
two_in_rain_less2 <- ifelse(ST4_pred$value < log(50.8), 1, NA)
prob_map_vals <- apply(sim_prob, 1, mean)

par(mfrow=c(1,3))
plot(rasterFromXYZ(cbind(coords, two_in_rain_f)), main="rain > 2\" in forecast", cex.main=1.5, legend=F, cex.axis=1.5)
plot(rasterFromXYZ(cbind(coords,prob_map_vals)), main = "prob map of rain > 2\"",
     cex.main=1.5, zlim=c(0,1), cex.axis=1.5)
# addRasterLegend(rasterFromXYZ(cbind(coords, prob_map_vals)),zlim=c(0,1),cex.axis = 1.4)
plot(rasterFromXYZ(cbind(coords, two_in_rain_o)), main="rain > 2\" in observed", cex.main=1.5, legend=F, cex.axis=1.5)

par(mfrow=c(1,2))
plot(rasterFromXYZ(cbind(coords, prob_map_vals)), main="all rain probs")
hist(prob_map_vals, main="all rain probs")

par(mfrow=c(2,2))
plot(rasterFromXYZ(cbind(coords, prob_map_vals*two_in_rain_oNA)))#, main="rain > 2\" in observed\n plotted by prob of > 2\"")
hist(prob_map_vals*two_in_rain_oNA, main=NULL, xlab=NULL)#, main="prob of > 2\" associated with pts \nthat actually were over 2\"")

# par(mfrow=c(1,2))
plot(rasterFromXYZ(cbind(coords, prob_map_vals*two_in_rain_less2)))#, main="rain < 2\" in observed\n plotted by prob of < 2\"")
hist(prob_map_vals*two_in_rain_less2, main=NULL,xlab=NULL)#, main="prob of < 2\" associated with pts \nthat actually were under 2\"")

# dev.off()

after <- Sys.time()
after - before

## Dave plots for water folks, from 4/28/20 notes
# pdf("~/NAM-Model-Validation/pdf/prediction/NAM_ST4_with_error_sims.pdf")

min_p <- min(values(NAM_r+PW_mean), values(ST4_r), values(NAM_r), na.rm = T)
max_p <- max(values(NAM_r+PW_mean), values(ST4_r), values(NAM_r), na.rm = T)
for (i in (1:3)+1) {
  sim <- cbind(coords, simvals[,i])
  sim_r <-rasterFromXYZ(sim)
  min_p <- min(min_p, values(NAM_r+sim_r), na.rm = T)
  max_p <- max(max_p, values(NAM_r+sim_r), na.rm = T)
}
par(mfrow=c(2,3))
plot(NAM_r + PW_mean, zlim=c(min_p,max_p), main="NAM forecast"); US(add=T, col="lightgray")
plot(ST4_r, zlim=c(min_p,max_p), main="ST4 observed"); US(add=T, col="lightgray")
plot(NAM_r, zlim=c(min_p,max_p), main="NAM bias adj"); US(add=T, col="lightgray")
for (i in (1:3)+1) {
  sim <- cbind(coords, simvals[,i])
  sim_r <-rasterFromXYZ(sim)
  plot(NAM_r+sim_r, main=paste("NAM bias adj + error field",i),zlim=c(min_p,max_p))
  US(add=T, col="lightgray")
}

#plot with c(0,7 bounds)
par(mfrow=c(2,3))
plot(NAM_r + PW_mean, zlim=c(0,7), main="NAM forecast"); US(add=T, col="lightgray")
plot(ST4_r, zlim=c(0,7), main="ST4 observed"); US(add=T, col="lightgray")
plot(NAM_r, zlim=c(0,7), main="NAM bias adj"); US(add=T, col="lightgray")
for (i in (1:3)+1) {
  sim <- cbind(coords, simvals[,i])
  sim_r <-rasterFromXYZ(sim)
  temp <- NAM_r + sim_r
  values(temp)[values(temp)>=7] = 7
  plot(temp, main=paste("NAM bias adj + error field",i),zlim=c(0,7))
  US(add=T, col="lightgray")
}

# dev.off()
