# Posterior for the pointwise (PW) mean field

library(geoR)
library(fields)
library(raster)

source("~/NAM-Model-Validation/scripts/mosaicList.R")

# CHANGE THIS FOR LOADING INFORMED VS FLAT PRIOR!
# loads post_cov_mtx
# "/work/dragonstooth/walsh124/post_cov_mtx_expntl_informedPrior.RData" # geq20
# "/work/dragonstooth/walsh124/post_cov_mtx_expntl_informedPrior15.RData" # geq15
load("~/NAM-Model-Validation/RData/post_cov_mtx_expntl.RData")

# loads sum_prec_mtx
load("~/NAM-Model-Validation/RData/sum_prec_mtx_expntl.RData")
data.prec.mtx <- sum_prec_mtx#solve(sum_cov_mtx)

# Convert PWmean raster to data frame for ggplot
PWmean <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_mean.grd")
PWmean1 <- as(PWmean, "SpatialPixelsDataFrame")
PWM1_df <- as.data.frame(PWmean1)
# write.csv(PWM1_df, file = "~/NAM-Model-Validation/csv/PWM1_df.csv")

## Obtain the posterior mean, m

all_prc_files <- list.files("/work/dragonstooth/walsh124/myMLE_bigprecs", full.names = T)
EF_files <- list.files("~/NAM-Model-Validation/csv/error_df_sqrt/subtractPWmeanF", full.names = T)
sum_all_prcXy <- rep(0, nrow(PWM1_df))

for(i in 1:length(all_prc_files)){
  
  print(i)
  
  ## these load all_prc_mtx for i in 1:47
  load(all_prc_files[i])
  
  EF_csv <- read.csv(EF_files[i], row.names = 1)
  new_EF <- rep(0, nrow(PWM1_df))
  
  for (j in 1:nrow(EF_csv)) {
    k <- which(abs(PWM1_df$x - EF_csv$x[j]) < 1e-3 & abs(PWM1_df$y - EF_csv$y[j]) < 1e-3)
    new_EF[k] <- PWM1_df$value[k]
  }
  sum_all_prcXy <- sum_all_prcXy + all_prc_mtx %*% new_EF
}

# plot(rasterFromXYZ(cbind(PWM1_df$x, PWM1_df$y, sum_all_prcXy)))
new_m <- post_cov_mtx %*% sum_all_prcXy
PW_post_new <- rasterFromXYZ(cbind(PWM1_df$x, PWM1_df$y, new_m))
writeRaster(PW_post_new, "~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_post_informPrior",
            overwrite=T)

par(mfrow=c(1,2))
z = range(values(PWmean), na.rm = T)
plot(PW_post_new, main = "new m", zlim=c(-3,3))
plot(PWmean, main = "old m", zlim=z)

#equivalent, alternative way to get the posterior PWmean
m = 0 + post_cov_mtx %*% data.prec.mtx %*% PWM1_df[,1]

PW_post <- rasterFromXYZ(cbind(PWM1_df[,2:3],m))
all.equal(values(PW_post), values(PW_post_new))

PW_data <- rasterFromXYZ(PWM1_df[,c(2,3,1)])
writeRaster(PW_post, "~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_post_15")

brks <- c(-5, -2, seq(-1,1,by=0.25), 2, 5)
brks2 <- c(-10, -4, seq(-2,2,by=0.5), 4, 10)
pdf("~/NAM-Model-Validation/pdf/data_vs_post_PWmean_expntl_newer.pdf", pointsize = 8)

par(mfrow=c(1,2), mar=c(5,4,4,6) + 0.1)

plot(PW_data, main="Data PW mean", breaks = brks, 
     col = c("blue", "green", cm.colors(length(brks)-5), "orange", "red"))  # col = cm.colors(length(brks)-1)
US(add=T, col="gray80")

plot(PW_post, main="Posterior PW mean (old m)", breaks = brks, 
     col = c("blue", "green", cm.colors(length(brks)-5), "orange", "red"))
US(add=T, col="gray80")

plot(PW_data, main="Data PW mean", breaks = brks, 
     col = c("blue", "green", cm.colors(length(brks)-5), "orange", "red"))  # col = cm.colors(length(brks)-1)
US(add=T, col="gray80")

plot(PW_post_new, main="Posterior PW mean (new m)", breaks = brks, 
     col = c("blue", "green", cm.colors(length(brks)-5), "orange", "red"))
US(add=T, col="gray80")

par(mfrow=c(1,1))
plot(PW_post - PW_data, main="Subtract PW means (old m)")
US(add=T, col="gray80")

plot(PW_post_new - PW_data, main="Subtract PW means (new m)")
US(add=T, col="gray80")

PW_post_df <- rasterToPoints(PW_post)
sds <- sqrt(diag(post_cov_mtx))
PW_post_sds <- rasterFromXYZ(cbind(PW_post_df[,1:2], sds))
# writeRaster(PW_post_sds, "~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_post_sds")
plot(PW_post_sds, main = "posterior SDs")
US(add=T, col="gray80")

plot(PW_post_sds^2, main = "posterior Vars")
US(add=T, col="gray80")

# writeRaster(PW_post/PW_post_sds, "~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_post_stdz")
plot(PW_post/PW_post_sds, main = "standardized mean (old m)",breaks = brks2, 
     col = c("blue", "green", cm.colors(length(brks2)-5), "orange", "red"))
US(add=T, col="gray80")

# writeRaster(PW_post_new/PW_post_sds, "~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_post_stdz_newm")
plot(PW_post_new/PW_post_sds, main = "standardized mean (new m)",breaks = brks2, 
     col = c("blue", "green", cm.colors(length(brks2)-5), "orange", "red"))
US(add=T, col="gray80")

dev.off()





# ## See which pairs of points have the highest covariance
# hi_cov <- which(post_cov_mtx > 5e-2)
# head(hi_cov, 30)
# 
# # Could map these to the lat lons from PWM1_df
# pairs_hicov <- cbind((hi_cov%/%nrow(PWM1_df) + 1), hi_cov%%nrow(PWM1_df))
# 
# 
# 
# ##### Update posterior mean map
# # Make it so that m = 0 + post_cov_mtx + sum_i=1^47{prec_mtx_i %*% EF_i}
# 
# PWPmeanr <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_mean.grd")
# PWPmeanr[!is.na(PWPmeanr)] <- 0
# # writeRaster(PWPmeanr, "~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_mean_0")
# PWPmeanr[PWPmeanr == 0] <- 1
# 
# # ~20 min with 6 cores
# sum_prec_EF <- matrix(0, nrow(PWmean1), ncol(PWmean1))
# for (i in 1:47) {
#   print(i)
#   EF <- mosaicList(c(error_files[i], "~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_mean_0.grd"))*PWPmeanr
#   plot(EF, main=i)
#   EFdf <- as.data.frame(as(EF, "SpatialPixelsDataFrame"))
#   load(paste0("~/NAM-Model-Validation/RData/all_prc_mtx/",ifelse(i < 10, paste0("0",i),i),".RData"))
#   sum_prec_EF <- sum_prec_EF + (all_prc_mtx %*% EFdf[,1])
# }
# load("~/NAM-Model-Validation/RData/post_cov_mtx.RData")
# load("~/NAM-Model-Validation/RData/sum_cov_mtx.RData")
# m_new <- post_cov_mtx %*% sum_prec_EF
# m_flat <- sum_cov_mtx %*% sum_prec_EF
# # writeRaster(rasterFromXYZ(cbind(EFdf[,2:3], m_new)), "~/NAM-Model-Validation/error_rasters_summary/PW_post_new")
# # writeRaster(rasterFromXYZ(cbind(EFdf[,2:3], m_flat)), "~/NAM-Model-Validation/error_rasters_summary/PW_post_flat")
# 
# #tomorrow: run the prediction code with the pwpostnew and pwpostflat; use the post_cov_mtx for pw_uq
# #have to rerun all of the LPVgrams, etc again?
# 
# 
# new_postPW <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_post_new.grd")
# old_postPW <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_post.grd")
# dataPW <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_mean.grd")
# flatPW <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_post_flat.grd")
# 
# min_v <- min(c(values(dataPW),values(flatPW),values(old_postPW),values(new_postPW)), na.rm = T)
# max_v <- max(c(values(dataPW),values(flatPW),values(old_postPW),values(new_postPW)), na.rm = T)
# 
# pdf("~/NAM-Model-Validation/pdf/compare_PWmeans2.pdf")
# par(mfrow=c(2,2))
# plot(new_postPW, main="new postPW",zlim=c(min_v,max_v))
# plot(old_postPW, main = "old postPW",zlim=c(min_v,max_v))
# plot(dataPW, main="dataPW", xlim=extent(flatPW)[1:2], ylim=extent(flatPW)[3:4],zlim=c(min_v,max_v))
# plot(flatPW, main="flatPW",zlim=c(min_v,max_v))
# 
# plot(abs(new_postPW - old_postPW), main="abs(newpostPW - oldpostPW)")
# plot(abs(new_postPW - dataPW), main="abs(newpostPW - dataPW)")
# plot(abs(new_postPW - flatPW), main="abs(newpostPW - flatPW)")
# plot(abs(dataPW - flatPW), main="abs(dataPW - flatPW)")
# 
# plot((new_postPW - old_postPW), main="(newpostPW - oldpostPW)")
# plot((new_postPW - dataPW), main="(newpostPW - dataPW)")
# plot((new_postPW - flatPW), main="(newpostPW - flatPW)")
# plot((dataPW - flatPW), main="(dataPW - flatPW)")
# 
# dev.off()
# 
# dataPW01 <- dataPW
# dataPW01[dataPW01 <0] <- -1
# dataPW01[dataPW01 >0] <- 1
# plot(dataPW01)