# Posterior for the pointwise (PW) mean field

library(geoR)
library(fields)
library(raster)

### Prior: mu ~ N(0, \Sigma_0)) 
# Sigma_0 is from the PWmean where pixel counts are geq20 MLE params

load("~/NAM-Model-Validation/RData/PWmean_geq20_paramsMLE.RData") 

source("~/NAM-Model-Validation/scripts/mosaicList.R")

# Read in error counts
count_files <- list.files("~/NAM-Model-Validation/error_rasters_counts/", pattern = ".grd", full.names = T)
error_counts <- mosaicList(count_files)
error_files <- list.files("~/NAM-Model-Validation/error_rasters/", pattern = ".grd", full.names = T)
error_sum <- mosaicList(error_files)

# only look at >= 1
geq1 <- error_counts
geq1[geq1 < 1] <- NA
# geq1[geq1 >= 1] <- 1
# plot(error_sum/geq1)

# Convert raster to data frame for ggplot
PWmean1 <- as(error_sum/geq1, "SpatialPixelsDataFrame")
PWM1_df <- as.data.frame(PWmean1)
# write.csv(PWM1_df, file = "~/NAM-Model-Validation/csv/PWM1_df.csv")

## Prior mean: 0
x0 <- rep(0,nrow(PWM1_df))

# For the prior covariance structure using geq20 params
s <- PWM1_df[,c(2,3)]
D <- as.matrix(dist(s))
sigma2 <- MLtemp$sigmasq
phi <- MLtemp$phi
nu <- MLtemp$kappa
# plot(seq(0,10,0.01),sigma2*matern(seq(0,10,0.01),phi=phi,kappa=nu),type="l") #(1/phi)
priorcovmatrix <- sigma2*matern(D,phi=phi,kappa=nu) # + diag(tau2,nrow(t))
prior.prec.mtx <- solve(priorcovmatrix)

# loads post_cov_mtx
load("~/NAM-Model-Validation/RData/post_cov_mtx.RData")

# loads sum_cov_mtx
load("~/NAM-Model-Validation/RData/sum_cov_matx.RData")
data.prec.mtx <- solve(sum_cov_mtx)


## Obtain the posterior mean, m
m = 0 + post_cov_mtx %*% data.prec.mtx %*% PWM1_df[,1]

PW_post <- rasterFromXYZ(cbind(PWM1_df[,2:3],m))
PW_data <- rasterFromXYZ(PWM1_df[,c(2,3,1)])
# writeRaster(PW_post, "~/NAM-Model-Validation/error_rasters_summary/PW_post")

brks <- c(-5, -2, seq(-1,1,by=0.25), 2, 5)
pdf("~/NAM-Model-Validation/pdf/data_vs_post_PWmean_t.pdf", pointsize = 6)
par(mfrow=c(1,2), mar=c(5,4,4,6) + 0.1)
plot(PW_data, main="Data PW mean", breaks = brks, 
     col = c("blue", "green", cm.colors(length(brks)-5), "orange", "red"))  # col = cm.colors(length(brks)-1)
US(add=T, col="gray80")
plot(PW_post, main="Posterior PW mean", breaks = brks, 
     col = c("blue", "green", cm.colors(length(brks)-5), "orange", "red"))
US(add=T, col="gray80")
dev.off()

## See which pairs of points have the highest covariance
hi_cov <- which(post_cov_mtx > 5e-2)
head(hi_cov, 30)

# Could map these to the lat lons from PWM1_df
pairs_hicov <- cbind((hi_cov%/%nrow(PWM1_df) + 1), hi_cov%%nrow(PWM1_df))



##### Update posterior mean map
# Make it so that m = 0 + post_cov_mtx + sum_i=1^47{prec_mtx_i %*% EF_i}

PWPmeanr <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_mean.grd")
PWPmeanr[!is.na(PWPmeanr)] <- 0
# writeRaster(PWPmeanr, "~/NAM-Model-Validation/error_rasters_summary/PW_mean_0")
PWPmeanr[PWPmeanr == 0] <- 1

# ~20 min with 6 cores
sum_prec_EF <- matrix(0, nrow(PWmean1), ncol(PWmean1))
for (i in 1:47) {
  print(i)
  EF <- mosaicList(c(error_files[i], "~/NAM-Model-Validation/error_rasters_summary/PW_mean_0.grd"))*PWPmeanr
  plot(EF, main=i)
  EFdf <- as.data.frame(as(EF, "SpatialPixelsDataFrame"))
  load(paste0("~/NAM-Model-Validation/RData/all_prc_mtx/",ifelse(i < 10, paste0("0",i),i),".RData"))
  sum_prec_EF <- sum_prec_EF + (all_prc_mtx %*% EFdf[,1])
}
load("~/NAM-Model-Validation/RData/post_cov_mtx.RData")
load("~/NAM-Model-Validation/RData/sum_cov_mtx.RData")
m_new <- post_cov_mtx %*% sum_prec_EF
m_flat <- sum_cov_mtx %*% sum_prec_EF
# writeRaster(rasterFromXYZ(cbind(EFdf[,2:3], m_new)), "~/NAM-Model-Validation/error_rasters_summary/PW_post_new")
# writeRaster(rasterFromXYZ(cbind(EFdf[,2:3], m_flat)), "~/NAM-Model-Validation/error_rasters_summary/PW_post_flat")

#tomorrow: run the prediction code with the pwpostnew and pwpostflat; use the post_cov_mtx for pw_uq
#have to rerun all of the LPVgrams, etc again?


new_postPW <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_post_new.grd")
old_postPW <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_post.grd")
dataPW <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_mean.grd")
flatPW <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_post_flat.grd")

min_v <- min(c(values(dataPW),values(flatPW),values(old_postPW),values(new_postPW)), na.rm = T)
max_v <- max(c(values(dataPW),values(flatPW),values(old_postPW),values(new_postPW)), na.rm = T)

pdf("~/NAM-Model-Validation/pdf/compare_PWmeans.pdf")
par(mfrow=c(2,2))
plot(new_postPW, main="new postPW",zlim=c(min_v,max_v))
plot(old_postPW, main = "old postPW",zlim=c(min_v,max_v))
plot(dataPW, main="dataPW", xlim=extent(flatPW)[1:2], ylim=extent(flatPW)[3:4],zlim=c(min_v,max_v))
plot(flatPW, main="flatPW",zlim=c(min_v,max_v))

plot(abs(new_postPW - old_postPW), main="abs(newpostPW - oldpostPW)")
plot(abs(new_postPW - dataPW), main="abs(newpostPW - dataPW)")
plot(abs(new_postPW - flatPW), main="abs(newpostPW - flatPW)")
plot(abs(dataPW - flatPW), main="abs(dataPW - flatPW)")

plot((new_postPW - old_postPW), main="(newpostPW - oldpostPW)")
plot((new_postPW - dataPW), main="(newpostPW - dataPW)")
plot((new_postPW - flatPW), main="(newpostPW - flatPW)")
plot((dataPW - flatPW), main="(dataPW - flatPW)")

dev.off()

dataPW01 <- dataPW
dataPW01[dataPW01 <0] <- -1
dataPW01[dataPW01 >0] <- 1
plot(dataPW01)