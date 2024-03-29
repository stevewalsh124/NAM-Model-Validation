# Compare margins of error from UQ with variability among test storms

library(raster)

pdf("pdf/rain_intervals.pdf")
# png("png/rain_intervals.png", 2100, 2100, res = 350)
# par(mfrow=c(3,2), mar=c(5, 4, 4, 2) + 0.1)

pred_names <- c("Alberto","Florence","Gordon","Michael","Barry","Dorian")

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
  # hist(rain_scen)
  
  # look at entire range for each grid point
  rain_ints <- apply(rain_scen, 1, function(X){diff(range(X))})
  # hist(rain_ints, breaks=40)
  
  # look at 95% CI at each grid point
  rain_ints95 <- apply(rain_scen, 1, function(X){quantile(X,0.975)-quantile(X,0.025)})
  rain_UBs <- apply(rain_scen, 1, function(X){quantile(X,0.975)})
  rain_LBs <- apply(rain_scen, 1, function(X){quantile(X,0.025)})  
  rain_UB99s <- apply(rain_scen, 1, function(X){quantile(X,0.995)})
  rain_LB99s <- apply(rain_scen, 1, function(X){quantile(X,0.005)})
  rain_SDs <-  apply(rain_scen, 1, sd)
  diffins <- rain_UBs - rain_LBs
  # hist(rain_ints95, breaks=40)
  
  # Find the 95% upper bound for these 95% CIs
  UBs_int95s[ste] <- quantile(rain_ints95,0.95)
  # Find the average length of these 95% CIs
  avg_int95s[ste] <- mean(rain_ints95)
  # Find the 5% lower bound for these 95% CIs
  LBs_int95s[ste] <- quantile(rain_ints95,0.05)
  # Find the maximum CI
  max_int95s[ste] <- max(rain_ints95)
  # Look how forecasted value affects CIs
  # plot(rowMeans(rain_scen)/25.4, rain_ints95/25.4/2, main = paste0("Storm ",ste,": ",pred_names[ste]),
  #      xlab = "predictive mean (in)", ylab = "margin of error (in)")
  
  plot(rowMeans(rain_scen)/25.4, rain_LBs/25.4, main = paste0("Storm ",ste,": ",pred_names[ste]),
       xlab = "predictive mean (in)", ylab = "bounds and ST4 (in)", ylim=c(0,max(rain_UBs/25.4)))
  points(rowMeans(rain_scen)/25.4, rain_UBs/25.4, col="red")
  points(rowMeans(rain_scen)/25.4, NAM_pred$value^2/25.4, col="blue")
  points(rowMeans(rain_scen)/25.4, ST4_pred$value^2/25.4, col="GREEN", pch="*")
  points(rowMeans(rain_scen)/25.4, rain_LBs/25.4, pch=".")
  points(rowMeans(rain_scen)/25.4, rain_UBs/25.4, col="red", pch=".")
  points(rowMeans(rain_scen)/25.4, rain_LB99s/25.4, pch=".")
  points(rowMeans(rain_scen)/25.4, rain_UB99s/25.4, col="red", pch=".")
  points(rowMeans(rain_scen)/25.4, NAM_pred$value^2/25.4, col="blue", pch=".")
  legend("topleft", legend = c("ub","lb","NAM","ST4"), col=c("black","red","blue","green"),pch=1)
  # lines((0:200)/10, sqrt((0:200)/10)*1.5, col="green")
  
  # Gather the mean amount of rain for each storm in buffer region (25.4mm per inch)
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
  
  # par(mfrow=c(1,2), mar=c(5, 4, 4, 2) + 0.1)
  # plot(rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, rowMeans(rain_scen)/25.4)), 
  #      main = paste0("storm ",ste,": pred mean (in)"))
  # plot(rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, rain_ints95/25.4/2)), 
  #      main = paste0("storm ",ste,": margin of error (in)"))
}

# Weighted average of 95% upper bound for length of 95% confidence intervals
# Also, Weighted average of 95% confidence intervals
# 25.4mm for every inch
sum(npixs*UBs_int95s)/sum(npixs)/25.4 # 4.322 inches
sum(npixs*avg_int95s)/sum(npixs)/25.4 # 1.605 inches
max(max_int95s)/25.4 # 13.215 inches

# Margin of error is half the length of CI
sum(npixs*UBs_int95s)/sum(npixs)/25.4/2 # 2.161 inches
sum(npixs*avg_int95s)/sum(npixs)/25.4/2 # 0.803 inches
max(max_int95s)/25.4/2 # 6.607 inches

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

dev.off()

plot(rowMeans(rain_scen)/25.4, rain_ints95/25.4/2, main = paste0("Storm ",ste,": ",pred_names[ste]),
     xlab = "predictive mean (in)", ylab = "margin of error (in)")
plot(rowMeans(rain_scen)/25.4, rain_UBs/25.4, main = paste0("Storm ",ste,": ",pred_names[ste]),
     xlab = "predictive mean (in)", ylab = "ub 97.5% (in)")
plot(rowMeans(rain_scen)/25.4, rain_LBs/25.4, main = paste0("Storm ",ste,": ",pred_names[ste]),
     xlab = "predictive mean (in)", ylab = "lb 2.5% (in)")

plot(rowMeans(rain_scen)/25.4, rain_LBs/25.4, main = paste0("Storm ",ste,": ",pred_names[ste]),
     xlab = "predictive mean (in)", ylab = "bounds and ST4 (in)", ylim=c(0,max(rain_UBs/25.4)))
points(rowMeans(rain_scen)/25.4, rain_UBs/25.4, col="red")
points(rowMeans(rain_scen)/25.4, NAM_pred$value^2/25.4, col="blue")
points(rowMeans(rain_scen)/25.4, ST4_pred$value^2/25.4, col="GREEN", pch="*")
points(rowMeans(rain_scen)/25.4, rain_LBs/25.4, pch=".")
points(rowMeans(rain_scen)/25.4, rain_UBs/25.4, col="red", pch=".")
points(rowMeans(rain_scen)/25.4, NAM_pred$value^2/25.4, col="blue", pch=".")
legend("topleft", legend = c("ub","lb","NAM","ST4"), col=c("black","red","blue","green"),pch=1)

plot(rowMeans(rain_scen)/25.4, rain_SDs/25.4, main = "predictive SD vs predictive mean")
