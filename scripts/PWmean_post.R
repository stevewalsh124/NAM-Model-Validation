# Posterior for the pointwise (PW) mean field

library(geoR)
library(fields)

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
writeRaster(PW_post, "~/NAM-Model-Validation/error_rasters_summary/PW_post")

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
