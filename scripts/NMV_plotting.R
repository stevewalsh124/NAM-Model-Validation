library(raster) #raster, extent
library(fields) #US()

black <- rgb(c(0,0,0)/255,1)
orange <- rgb(230/255,159/255,0,1)
# skyblue <- rgb(86/255,180/255,233/255,1)
# bluegreen <- rgb(0,158/255,115/255,1)
yellow <- rgb(240/255,228/255,66/255,1)
blue <- rgb(0,114/255,178/255,1)
verm <- rgb(213/255, 94/255, 0,1)
redpurp <- rgb(204/255, 121/255, 167/255, 1)
plot(1:10, col=c(black,yellow,blue, redpurp,verm))
###################################
# Figure 1: landfalls by location #
###################################

all_content = readLines("csv/US_20042017_landfalls_mod.csv")
skip_second = all_content[-6] #6th storm had no significant landfall
landfalls = read.csv(textConnection(skip_second), header = T, stringsAsFactors = FALSE)
colnames(landfalls) <- c("Year","Storm","Landfall Time (UTC)","BT.Lat","BTLon","SSHWS")
landfalls$inten <- factor(ifelse(landfalls$SSHWS%in%c("TS","1","2"),"TS,1,2","3,4,5"))
landfalls$inten_alt <- factor(ifelse(landfalls$SSHWS%in%c("TS"),"TS",
                                     ifelse(landfalls$SSHWS%in%c("1","2"),"1,2","3,4,5")))
landfalls$loc   <- factor(ifelse(landfalls$BT.Lat>32,"ATL",
                                 ifelse(landfalls$BTLon>-88,"FL","GULF")))
landfalls$SSHWS=factor(landfalls$SSHWS, levels = levels(factor(landfalls$SSHWS))[c(5,1:4)])
landfalls$inten=factor(landfalls$inten , levels=levels(landfalls$inten)[c(2,1)])
landfalls$inten_alt=factor(landfalls$inten_alt , levels=levels(landfalls$inten_alt)[c(3,1,2)])


levels(landfalls$SSHWS)
stormsTS<- which(landfalls$SSHWS=="TS")
storms1 <- which(landfalls$SSHWS==1)
storms2 <- which(landfalls$SSHWS==2)
storms3 <- which(landfalls$SSHWS==3)
storms4 <- which(landfalls$SSHWS==4)
storms5 <- which(landfalls$SSHWS==5)
storms_TS_1_2 <- c(stormsTS,storms1,storms2)
storms_3_4_5  <- c(storms3,storms4,storms5)

stormsATL <- which(landfalls$BT.Lat>32)
stormsFL  <- which(landfalls$BT.Lat<32&landfalls$BTLon> -88)
stormsGULF<- which(landfalls$BT.Lat<32&landfalls$BTLon< -88)
dim(landfalls[stormsATL,]);dim(landfalls[stormsFL,]);dim(landfalls[stormsGULF,])

class(landfalls$BT.Lat) <- "numeric"
class(landfalls$BTLon) <- "numeric"

g_byint <- ggplot() + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  geom_point(data = landfalls[c(stormsTS),], aes(x=BTLon, y = BT.Lat, color=" TS"))+
  geom_point(data = landfalls[c(storms1,storms2),], aes(x=BTLon, y = BT.Lat, color="1/2"))+
  geom_point(data = landfalls[c(storms3,storms4,storms5),], aes(x=BTLon, y = BT.Lat, color="3/4/5"))+
  labs(title ="Landfalls by Intensity", x = "Longitude", y="Latitude") + 
  coord_fixed(xlim=c(-100, -68), ylim = c(24,45), ratio = 1)+
  scale_colour_manual(name="Intensity", values = c("green","blue","red"))

g_byint <- ggplot() + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  geom_point(data = landfalls[c(stormsTS),], 
             aes(x=BTLon, y = BT.Lat, color=" TS"),size=4,shape=15,fill=skyblue)+
  geom_point(data = landfalls[c(storms1,storms2),], 
             aes(x=BTLon, y = BT.Lat, color="1/2"), size=4,shape=16,fill=yellow)+
  geom_point(data = landfalls[c(storms3,storms4,storms5),], 
             aes(x=BTLon, y = BT.Lat, color="3/4/5"),size=4,shape=18,fill=verm)+
  labs(title ="Landfalls by Intensity", x = "Longitude", y="Latitude") + 
  coord_fixed(xlim=c(-98, -70), ylim = c(24.5,42), ratio = 1)+
  scale_colour_manual(name="Intensity", values = c(skyblue,yellow,verm),
                      guide = guide_legend(override.aes = list(shape = c(15,16,18), 
                                                               color = c(skyblue,yellow,verm))))
  
  # scale_shape_manual(values=c(16,17,18),color=c("green","blue","red"))+

g_byint

png(paste0("~/NAM-Model-Validation/png/g_byint.png"),width=4000, height=2400, res=700)
g_byint
dev.off()


##################################
# Figure 3: Plot NAM, ST4, error #
##################################

source("scripts/multiplot.R")

# run SqrtPrecipVariograms.R with:
# * write.pdf <- F
# * line 40, chg to: storm.dirs <- list.dirs("~/NAMandST4", recursive = F)[47] eg for Nate
# * pick subtractPWmean T or F
g1= ggplot(aes(x=x,y=y,fill=value),data=ST4_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max)) + 
  labs(title =paste(paste0(str_to_title(storm_name),": Stage IV Data")),x = "Longitude", y="Latitude") + 
  coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                     max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
              ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                     max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)

g2= ggplot(aes(x=x,y=y,fill=value),data=NAM_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max)) +
  labs(title =paste0(str_to_title(storm_name),": NAM Forecast"), x = "Longitude", y="Latitude") + 
  coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                     max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
              ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                     max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)

g3= ggplot(aes(x=x,y=y,fill=value),data=error_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white") +
  labs(title =paste(paste0(str_to_title(storm_name),": Stage IV - NAM")),
       x = "Longitude", y="Latitude")+
  coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                     max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
              ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                     max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)

ggsave(
  "png/NAM_ST4_error_gg.png",
  multiplot(g1,g2,g3,cols=3),
  width = 6.5,
  height = 2.5,
  dpi = 300
)

png(paste0("~/NAM-Model-Validation/png/NAM_ST4_error_", storm_year, storm_name,".png"),
    width=1800, height=525, res=175)
multiplot(g1, g2, g3, cols=3)
dev.off()

png(paste0("~/NAM-Model-Validation/png/NAM_ST4_error_", storm_year, storm_name,"in.png"),
    width=6.5, height=2.5, units = "in", res=100)
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
dataPWmean <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_mean.grd")
flatPWmean <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_post_newm.grd")

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







# # older PWmean, PWvar, PWstandardizederror
# par(mfrow=c(1,3),mar=c(5,4,4,2)+.1)
# pwms <- vars <- stds <- list()
# for (n in 5) { 
#   # par(mfrow=c(1,2))
#   Nmap <- error_counts
#   Nmap[Nmap < n] <- NA
#   # Nmap[Nmap >= n] <- 1
#   
#   pwms <- PW_mean_geq5 <- error_sum/Nmap
#   
#   val <- max(abs(range(values(pwms), na.rm = T)))
#   plot(pwms, main = paste0("PW mean when n >= ",n),
#        col=cm.colors(10),
#        zlim=c(-val,val),
#        xlim=extent(PW_mean)[1:2], ylim=extent(PW_mean)[3:4], cex.axis=2); US(add=T, col="gray75")
#   
#   vars <- S2_geq5 <- ((error_sum_sq - (error_sum*error_sum)/Nmap)/(Nmap-1))*mask.regrid
#   plot(S2_geq5, main=paste0("Pointwise Variance when n >= ", n),
#        xlim=extent(PW_mean)[1:2], ylim=extent(PW_mean)[3:4], cex.axis=2);US(add=T, col="gray49")
#   
#   stds <- std_error_geq5 <- PW_mean_geq5/(sqrt(S2_geq5/Nmap))
#   # png("stdgeq5.png",width = 1000, height = 1000) #cex.main = 2.5, cex.axis = 2
#   val <- max(abs(range(values(stds), na.rm = T)))
#   plot(stds, main = paste0("Std err when n >= ",n), legend=F,
#        col=cm.colors(10),
#        xlim=extent(PW_mean)[1:2], ylim=extent(PW_mean)[3:4],zlim=c(-val,val),
#         cex.axis=2); US(add=T, col="gray75")
#   plot(stds, legend.only=TRUE,
#        legend.width=2, legend.shrink=1, col=cm.colors(10),
#        axis.args=list(at=c(-val,val),
#          #labels=seq(r.range[1]-1, r.range[2], 5), 
#          cex.axis=1.45))#,
#   # legend.args=list(text='Error counts', side=4, font=2, line=2.5, cex=1))
# }



# newer PWmean, PWvar, PWstandardizederror
PW_post_new <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_post_newm.grd")
PW_post_sds <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_post_sds.grd")

png("~/NAM-Model-Validation/png/NMV_post_mean_sds_stdz.png", width = 1200*1.5, height = 400*1.5)
par(mfrow=c(1,3),mar=c(6,6,4,11)+.1)
lmag <- 2.5 #letter magnification
val <- ceiling(max(abs(range(values(PW_post_new), na.rm = T))/0.1))*0.1
plot(PW_post_new, 
     col=c("blue",cm.colors(length(seq(-3,3,.5))-1),"red"),
     breaks=c(-val,seq(-3,3,.5),val), xlab="longitude", 
     ylab="latitude", cex.axis=lmag, cex.lab=lmag, legend=F)
plot(PW_post_new, 
     col=c("blue",cm.colors(length(seq(-3,3,.5))-1),"red"),
     breaks=c(-val,seq(-3,3,.5),val),
     legend.width=2, legend.shrink=1,
     legend.only=T, axis.args=list(cex.axis=lmag))
US(add=T, col="gray75")

plot(PW_post_sds, xlab="longitude", ylab="latitude", cex.axis=lmag, cex.lab=lmag, legend=F)
plot(PW_post_sds, legend.only=T, legend.width=2, legend.shrink=1, 
     axis.args=list(cex.axis=lmag))
US(add=T, col="gray75")

val_stdz <- ceiling(max(abs(range(values(PW_post_new/PW_post_sds), na.rm = T))/0.1))*0.1
plot(PW_post_new/PW_post_sds, 
     col=c("blue",cm.colors(length(seq(-3,3,.5))-1),"red"),
     breaks=c(-val_stdz,seq(-3,3,.5),val_stdz), xlab="longitude", 
     ylab="latitude", cex.axis=lmag, cex.lab=lmag, legend=F)
plot(PW_post_new/PW_post_sds, legend.only=TRUE,
     legend.width=2, legend.shrink=1,
     axis.args=list(cex.axis=lmag),
     col=c("blue",cm.colors(length(seq(-3,3,.5))-1),"red"),
     breaks=c(-val_stdz,seq(-3,3,.5),val_stdz))
US(add=T, col="gray75")
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