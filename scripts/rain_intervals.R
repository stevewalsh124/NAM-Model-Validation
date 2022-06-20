# Compare margins of error from UQ with variability among test storms

library(raster)

par(mfrow=c(2,2))
UBs_int95s <- LBs_int95s <- avg_int95s <- max_int95s <- 
  diffs_in <- diffs_in_gr0 <- diffs_in_gr1 <- diffs_in_gp1i <- diffs_in_gp3i <- 
  act_rain <- act_rain_gr0 <- act_rain_gr1 <- act_rain_gp1i <- 
  act_rain_UB <- act_rain_UB_gp1i <- c()
npixs <- c()
for (ste in 1:6) {
  print(ste)
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",ste,"nopw.RData"))
  
  # 1000 rain values for each grid point
  rain_scen <- ifelse(simvals+NAM_pred$value>0, simvals+NAM_pred$value, 0)^2
  hist(rain_scen)
  
  # look at entire range for each grid point
  rain_ints <- apply(rain_scen, 1, function(X){diff(range(X))})
  hist(rain_ints, breaks=40)
  
  # look at 95% CI at each grid point
  rain_ints95 <- apply(rain_scen, 1, function(X){quantile(X,0.975)-quantile(X,0.025)})
  hist(rain_ints95, breaks=40)
  
  # Find the 95% upper bound for these 95% CIs
  UBs_int95s[ste] <- quantile(rain_ints95,0.95)
  # Find the average length of these 95% CIs
  avg_int95s[ste] <- mean(rain_ints95)
  # Find the 5% lower bound for these 95% CIs
  LBs_int95s[ste] <- quantile(rain_ints95,0.05)
  # Find the maximum CI
  max_int95s[ste] <- max(rain_ints95)
  # Look how forecasted value affects CIs
  plot(NAM_pred$value^2, rain_ints95)
  
  # Gather the mean amount of rain for each storm in buffer region
  act_rain[ste] <- mean(ST4_pred$value^2)/25.4
  act_rain_gr0[ste] <- mean(ST4_pred$value[ST4_pred$value>0]^2)/25.4
  act_rain_gr1[ste] <- mean(ST4_pred$value[ST4_pred$value>1]^2)/25.4
  act_rain_gp1i[ste] <- mean(ST4_pred$value[ST4_pred$value>2.54]^2)/25.4
  act_rain_UB[ste] <-  quantile(ST4_pred$value^2, 0.95)/25.4
  act_rain_UB_gp1i[ste] <-  quantile(ST4_pred$value[ST4_pred$value>2.54]^2, 0.95)/25.4
  
  # Gather the differences between obs and forecast
  
  # Gather number of grid points to weight the average
  npixs[ste] <- length(NAM_pred$value)
  
  print(paste("max NAM (in):",max(NAM_pred$value)^2/25.4))
  print(paste("max ST4 (in):",max(ST4_pred$value)^2/25.4))
  print(paste("max abs diff (in):",max(abs(ST4_pred$value^2-NAM_pred$value^2)/25.4)))
  diffs_in <- c(diffs_in, abs(ST4_pred$value^2-NAM_pred$value^2)/25.4)
  diffs_in_gr0 <- c(diffs_in_gr0, abs(ST4_pred$value[ST4_pred$value>0]^2-NAM_pred$value[ST4_pred$value>0]^2)/25.4)
  diffs_in_gr1 <- c(diffs_in_gr1, abs(ST4_pred$value[ST4_pred$value>1]^2-NAM_pred$value[ST4_pred$value>1]^2)/25.4)
  diffs_in_gp1i <- c(diffs_in_gp1i, abs(ST4_pred$value[ST4_pred$value>2.54]^2-NAM_pred$value[ST4_pred$value>2.54]^2)/25.4)
  diffs_in_gp3i <- c(diffs_in_gp3i, abs(ST4_pred$value[ST4_pred$value>2.54*3]^2-NAM_pred$value[ST4_pred$value>2.54*3]^2)/25.4)
}

# Weighted average of 95% upper bound for length of 95% confidence intervals
# Also, Weighted average of 95% confidence intervals
# 25.4mm for every inch
sum(npixs*UBs_int95s)/sum(npixs)/25.4 # 4.322 inches
sum(npixs*avg_int95s)/sum(npixs)/25.4 # 1.605 inches
max(max_int95s)/25.4 # 1.605 inches

# Margin of error is half the length of CI
sum(npixs*UBs_int95s)/sum(npixs)/25.4/2 # 4.322 inches
sum(npixs*avg_int95s)/sum(npixs)/25.4/2 # 1.605 inches
max(max_int95s)/25.4/2 # 1.605 inches

# Look at 95%, 99% UBs and max absolute difference of NAM & ST4
quantile(diffs_in, 0.95)
quantile(diffs_in, 0.99)
max(diffs_in)

# Look at actual rain amount under different thresholds
act_rain
act_rain_gr0
act_rain_gr1
act_rain_gp1i

act_rain_UB
act_rain_UB_gp1i

# Average absolute difference b/w observed and forecasted precip,
# with different thresholds for observed precip (>0, >1mm, >.1in, >.3in)
mean(diffs_in)
mean(diffs_in_gr0)
mean(diffs_in_gr1)
mean(diffs_in_gp1i)
mean(diffs_in_gp3i)
