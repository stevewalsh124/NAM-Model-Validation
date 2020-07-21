library(raster) #raster, extent
library(fields) #US()

################
# Plot PWmeans #
################

# these generated from LogPrecipVariograms and PWmean_post, respectively
dataPWmean <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_mean.grd")
flatPWmean <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_post_flat.grd")

png("~/NAM-Model-Validation/png/dataPWmean.png", width=1000, height=700)
plot(dataPWmean, cex.axis=2, 
     xlim=extent(flatPWmean)[1:2], ylim=extent(flatPWmean)[3:4])
US(add=T, col="darkgray")
dev.off()

png("~/NAM-Model-Validation/png/flatPWmean.png", width=1000, height=700)
plot(flatPWmean, cex.axis=2, 
     xlim=extent(flatPWmean)[1:2], ylim=extent(flatPWmean)[3:4])
US(add=T, col="darkgray")
dev.off()


#####################
# post vs MLE plots #
#####################

# from GibbsSamplerHurrRegrNA.R
load("~/NAM-Model-Validation/RData/Gibbs_flatPW.RData")

# Lower and upper bound extra space to see end of distbns
up_b <- c(0.5,0.7,0.1); lw_b <- c(0.3,0.3,0.1)

## Plot MLE vs posterior estimates for the 47 storms
png("~/NAM-Model-Validation/png/NMV_compareMLEpost.png", width = 1400, height = 700)
par(mar=c(4.5,1.5,1.5,1.5), mfrow=c(1,3))
for (i in 1:P) {
  plot(0, 0, col = "white", ylab = "", 
       xlim=c(min(theta_hat[,i], emp_thetaMED[,i])-lw_b[i],
              max(theta_hat[,i], emp_thetaMED[,i])+up_b[i]),
       ylim=c(0,2), yaxt='n', cex = 3, cex.axis =2, cex.lab = 2,
       xlab = bquote(theta[.(i)]~"bottom and"~hat(theta)[.(i)]~"top"))
  # mult_seg <- data.frame(x0 = c(0.7, 0.2, - 0.9, 0.5, - 0.2),    # Create data frame with line-values
  #                        y0 = c(0.5, 0.5, 0.6, - 0.3, 0.4),
  #                        x1 = c(0, 0.2, 0.4, - 0.3, - 0.6),
  #                        y1 = c(- 0.1, 0.3, - 0.6, - 0.8, 0.9))
  
  segments(x0 = emp_thetaMED[,i],                            # Draw multiple lines
           y0 = 0,
           x1 = theta_hat[,i],
           y1 = 1, col = loc_int$loc, lwd=3)
  points(x=emp_thetaMED[,i], y=rep(0, length(emp_thetaMED[,i])), col= loc_int$loc)
  points(x=theta_hat[,i],    y=rep(1, length(emp_thetaMED[,i])), col= loc_int$loc)
  for (j in 1:N) {
    #top row, theta hats
    xseq <- seq(theta_hat[j,i]-3*sqrt(solve(hessians[[j]])[i,i]),
                theta_hat[j,i]+3*sqrt(solve(hessians[[j]])[i,i]),
                length.out = 1000)
    lines(xseq,smush[i]*dnorm(xseq,
                              theta_hat[j,i],
                              sqrt(solve(hessians[[j]])[i,i]))+1, col=loc_int$loc[j],
          lty = as.integer(loc_int$loc)[j], lwd=3)
    
    #bottom row, thetas from MCMC
    xseq <- seq(min(theta_burn[[j]][,i]),
                max(theta_burn[[j]][,i]), 
                length.out = 1000)
    lines(xseq, smush[i]*dnorm(xseq, 
                               mean(theta_burn[[j]][,i]), 
                               sd(theta_burn[[j]][,i])), col=loc_int$loc[j],
          lty = as.integer(loc_int$loc)[j], lwd=3)
  }
  legend("topright", 
         2, lwd=3,
         legend=c(levels(loc_int$loc)[unique(as.integer(loc_int$loc))]),
         col=unique(as.integer(loc_int$loc)), lty=unique(as.integer(loc_int$loc)), cex=2)
}
dev.off()

png("~/NAM-Model-Validation/png/NMV_compareMLEpost2.png", width = 1400, height = 700)
par(mar=c(4.5,1.5,1.5,1.5), mfrow=c(1,3))
for (i in 1:P) {
  plot(0, 0, col = "white", ylab = "", 
       xlim=c(min(theta_hat[,i], emp_thetaMED[,i])-lw_b[i],
              max(theta_hat[,i], emp_thetaMED[,i])+up_b[i]),
       ylim=c(0,2), yaxt='n', cex = 3, cex.axis =2, cex.lab = 2,
       xlab = bquote(theta[.(i)]~"bottom and"~hat(theta)[.(i)]~"top"))
  # mult_seg <- data.frame(x0 = c(0.7, 0.2, - 0.9, 0.5, - 0.2),    # Create data frame with line-values
  #                        y0 = c(0.5, 0.5, 0.6, - 0.3, 0.4),
  #                        x1 = c(0, 0.2, 0.4, - 0.3, - 0.6),
  #                        y1 = c(- 0.1, 0.3, - 0.6, - 0.8, 0.9))
  
  segments(x0 = emp_thetaMED[,i],                            # Draw multiple lines
           y0 = 0,
           x1 = theta_hat[,i],
           y1 = 1, col = loc_int$loc)
  points(x=emp_thetaMED[,i], y=rep(0, length(emp_thetaMED[,i])), col= loc_int$loc)
  points(x=theta_hat[,i],    y=rep(1, length(emp_thetaMED[,i])), col= loc_int$loc)
  for (j in 1:N) {
    #top row, theta hats
    xseq <- seq(theta_hat[j,i]-3*sqrt(solve(hessians[[j]])[i,i]),
                theta_hat[j,i]+3*sqrt(solve(hessians[[j]])[i,i]),
                length.out = 1000)
    lines(xseq,smush[i]*dnorm(xseq,
                              theta_hat[j,i],
                              sqrt(solve(hessians[[j]])[i,i]))+1, col=loc_int$loc[j],
          lty = as.integer(loc_int$loc)[j], lwd=3)
    
    #bottom row, thetas from MCMC
    xseq <- seq(min(theta_burn[[j]][,i]),
                max(theta_burn[[j]][,i]), 
                length.out = 1000)
    lines(xseq, smush[i]*dnorm(xseq, 
                               mean(theta_burn[[j]][,i]), 
                               sd(theta_burn[[j]][,i])), col=loc_int$loc[j],
          lty = as.integer(loc_int$loc)[j], lwd=3)
  }
  legend("topright", 
         2, lwd=3,
         legend=c(levels(loc_int$loc)[unique(as.integer(loc_int$loc))]),
         col=unique(as.integer(loc_int$loc)), lty=unique(as.integer(loc_int$loc)), cex=2)
}
dev.off()

png("~/NAM-Model-Validation/png/NMV_compareMLEpost3.png", width = 1400, height = 700)
par(mar=c(4.5,1.5,1.5,1.5), mfrow=c(1,3))
for (i in 1:P) {
  plot(0, 0, col = "white", ylab = "", 
       xlim=c(min(theta_hat[,i], emp_thetaMED[,i])-lw_b[i],
              max(theta_hat[,i], emp_thetaMED[,i])+up_b[i]),
       ylim=c(0,2), yaxt='n', cex = 3, cex.axis =2, cex.lab = 2,
       xlab = bquote(theta[.(i)]~"bottom and"~hat(theta)[.(i)]~"top"))
  
  segments(x0 = emp_thetaMED[,i],                            # Draw multiple lines
           y0 = 0,
           x1 = theta_hat[,i],
           y1 = 1, col = loc_int$loc, lwd=3)
  points(x=emp_thetaMED[,i], y=rep(0, length(emp_thetaMED[,i])), col= loc_int$loc)
  points(x=theta_hat[,i],    y=rep(1, length(emp_thetaMED[,i])), col= loc_int$loc)
  for (j in 1:N) {
    #top row, theta hats
    xseq <- seq(theta_hat[j,i]-3*sqrt(solve(hessians[[j]])[i,i]),
                theta_hat[j,i]+3*sqrt(solve(hessians[[j]])[i,i]),
                length.out = 1000)
    lines(xseq,smush[i]*dnorm(xseq,
                              theta_hat[j,i],
                              sqrt(solve(hessians[[j]])[i,i]))+1, col=loc_int$loc[j],
          lwd=3) #lty = as.integer(loc_int$loc)[j], 
    
    #bottom row, thetas from MCMC
    xseq <- seq(min(theta_burn[[j]][,i]),
                max(theta_burn[[j]][,i]), 
                length.out = 1000)
    lines(xseq, smush[i]*dnorm(xseq, 
                               mean(theta_burn[[j]][,i]), 
                               sd(theta_burn[[j]][,i])), col=loc_int$loc[j],
          lwd=3) #lty = as.integer(loc_int$loc)[j], 
  }
  legend("topright", 
         2, lwd=3,
         legend=c(levels(loc_int$loc)[unique(as.integer(loc_int$loc))]),
         col=unique(as.integer(loc_int$loc)), lty=unique(as.integer(loc_int$loc)), cex=2)
}
dev.off()