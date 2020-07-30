library(raster) #raster, extent
library(fields) #US()

########################
# Plot NAM, ST4, error #
########################

# run LogPrecipVariograms with desired storm number for storm.dirs (45 for harvey, etc)
# pick subtractPWmean T or F
png(paste0("~/NAM-Model-Validation/png/NAM_ST4_error_", storm_year, storm_name,".png"),
    width=760, height=260)
multiplot(g1, g2, g3, cols=3)
dev.off()

error_counts_plot <- error_counts*mask.regrid
error_counts_plot[error_counts_plot <= 0] <- NA
plot(error_counts_plot, 
     xlim = extent(PW_mean)[1:2], 
     ylim = extent(PW_mean)[3:4],
     cex.axis = 3)
US(add=T, col="gray")

png("~/NAM-Model-Validation/png/error_counts.png", height = 500, width = 700)
par(mfrow=c(1,1))
r <- raster(nrows=5, ncols=5, vals=1:25)
plot(error_counts_plot, 
     xlim = extent(PW_mean)[1:2], 
     ylim = extent(PW_mean)[3:4], legend=F, cex.axis=2)
r.range <- c(minValue(error_counts_plot), maxValue(error_counts_plot))
plot(error_counts_plot, legend.only=TRUE,
     legend.width=2, legend.shrink=1,
     axis.args=list(#at=seq(r.range[1]-1, r.range[2], 5),
                    #labels=seq(r.range[1]-1, r.range[2], 5), 
                    cex.axis=1.45))#,
     # legend.args=list(text='Error counts', side=4, font=2, line=2.5, cex=1))
US(add=T, col="darkgray")
dev.off()


################
# Plot PWmeans #
################

# these generated from LogPrecipVariograms and PWmean_post, respectively
dataPWmean <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_mean.grd")
flatPWmean <- raster("~/NAM-Model-Validation/error_rasters_summary/PW_post_flat.grd")

val <- max(abs(c(floor(4*minValue(pwms))/4, ceiling(4*maxValue(pwms))/4)))
plot(PW_mean, main = paste0("PW mean when n >= ",n),
     col=c("blue",cm.colors(length(seq(-1.5,1.5,.25)[abs(seq(-1.5,1.5,.25)) < val ])-1),"red"),
     breaks=c(-val,seq(-1.5,1.5,.25)[abs(seq(-1.5,1.5,.25)) < val ],val))
US(add=T, col="gray75")

plot(PW_mean, main = paste0("PW mean when n >= ",n),
     col=cm.colors(18),
     breaks=seq(-3,3,.25))
US(add=T, col="gray75")

plot(PW_mean, main="Pointwise Mean")
PW_mean_spdf <- as((PW_mean), "SpatialPixelsDataFrame")
PW_mean_df <- as.data.frame(PW_mean_spdf)
colnames(PW_mean_df) <- c("value", "x", "y")

g4= ggplot(aes(x=x,y=y,fill=value),data=PW_mean_df) +
  geom_tile() + theme_classic() +
  geom_polygon(data=subset(map_data("state"), region %in% regions),
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white") +
  # labs(title =paste("Pointwise Mean"))+
  coord_fixed(xlim=extent(PW_mean)[1:2],ylim=extent(PW_mean)[3:4], ratio = 1)+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"), 
        legend.title = element_text(size=18), legend.text = element_text(size=16))
plot(g4)

png("~/NAM-Model-Validation/png/PWmeanggplot.png", width=800, height=500)
plot(g4)
dev.off()


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















#PWmean, PWvar, PWstandardizederror
par(mfrow=c(1,3),mar=c(5,4,4,2)+.1)
pwms <- vars <- stds <- list()
for (n in 5) { 
  # par(mfrow=c(1,2))
  Nmap <- error_counts
  Nmap[Nmap < n] <- NA
  # Nmap[Nmap >= n] <- 1
  
  pwms <- PW_mean_geq5 <- error_sum/Nmap
  
  val <- max(abs(range(values(pwms), na.rm = T)))
  plot(pwms, main = paste0("PW mean when n >= ",n),
       col=cm.colors(10),
       zlim=c(-val,val),
       xlim=extent(PW_mean)[1:2], ylim=extent(PW_mean)[3:4], cex.axis=2); US(add=T, col="gray75")
  
  vars <- S2_geq5 <- ((error_sum_sq - (error_sum*error_sum)/Nmap)/(Nmap-1))*mask.regrid
  plot(S2_geq5, main=paste0("Pointwise Variance when n >= ", n),
       xlim=extent(PW_mean)[1:2], ylim=extent(PW_mean)[3:4], cex.axis=2);US(add=T, col="gray49")
  
  stds <- std_error_geq5 <- PW_mean_geq5/(sqrt(S2_geq5/Nmap))
  # png("stdgeq5.png",width = 1000, height = 1000) #cex.main = 2.5, cex.axis = 2
  val <- max(abs(range(values(stds), na.rm = T)))
  plot(stds, main = paste0("Std err when n >= ",n), legend=F,
       col=cm.colors(10),
       xlim=extent(PW_mean)[1:2], ylim=extent(PW_mean)[3:4],zlim=c(-val,val),
        cex.axis=2); US(add=T, col="gray75")
  plot(stds, legend.only=TRUE,
       legend.width=2, legend.shrink=1, col=cm.colors(10),
       axis.args=list(at=c(-val,val),
         #labels=seq(r.range[1]-1, r.range[2], 5), 
         cex.axis=1.45))#,
  # legend.args=list(text='Error counts', side=4, font=2, line=2.5, cex=1))
}





















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

png("~/NAM-Model-Validation/png/NMV_compareMLEpost2.png", width = 1500, height = 500)
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