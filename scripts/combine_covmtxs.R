## combine covmtxs: 
## run after aggregate_covmtxs_expntl.R
## run before PWmean_post_expntl.R

library(raster)

PWmean <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_mean_1to15.grd")
PWmean1 <- as(PWmean, "SpatialPixelsDataFrame")
PWM1_df <- as.data.frame(PWmean1)
colnames(PWM1_df) <- c("value", "x", "y")
a <- nrow(PWM1_df)

# IF YOU WANT A FLAT PRIOR, USE sum_prec_mtx <- matrix(0,a,a)
# the informed prior precision based on geq20 for sqrt expntl
# load("/work/dragonstooth/walsh124/bigPriorPrec.RData")
# sum_prec_mtx <- prior.prec.mtx
# sum_prec_mtx <- matrix(0, a, a) #is used instead for flat prior
sum_prec_mtx <- diag(.Machine$double.eps^0.5, a, a)

# this loads 47 big precision matrices for each storm 
precs <- list.files("~/NAM-Model-Validation/RData/myMLE_bigprecs_1to15", full.names = T)
for(i in 1:length(precs)){	
  print(i)	
  load(precs[i])	
  sum_prec_mtx <- sum_prec_mtx + all_prc_mtx	
}	

system.time(post_cov_mtx <- solve(sum_prec_mtx))
save(post_cov_mtx, file="~/NAM-post_cov_mtx_expntl_1to15_diagjitter.RData")

# times <- c(1140, 1155, 1156, 1208, 1209, 1210, 1224, 1254, 1301, 1313, 1314, 1410, 1420, 1422, 1443, 
#            1450, 1502, 1512, 1518, 1522, 1528, 1538, 1541, 1548, 1551, 1601, 1626, 1634, 1645, 1659, 
#            1714, 1721, 1735, 1736, 1807, 1808, 1853, 1901, 1915, 1915, 2042, 2116, 2147, 2242, 2502,
#            2508, 2724) - 1130
# plot(sort(n_pixels), times, main="Size vs Compute Time")
# plot(sort(n_pixels), sqrt(times), main="Size vs Sqrt Time")
# abline(lm(sqrt(times)~sort(n_pixels)))
