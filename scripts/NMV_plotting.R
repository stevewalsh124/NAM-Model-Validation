library(raster) #raster, extent
library(fields) #US()

black <- rgb(0,0,0,1)
orange <- rgb(230/255,159/255,0,1)
skyblue <- rgb(86/255,180/255,233/255,1)
bluegreen <- rgb(0,158/255,115/255,1)
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
  geom_point(data = landfalls[c(stormsTS),], 
             aes(x=BTLon, y = BT.Lat, color=" TS"),size=4,shape=15,fill=skyblue)+
  geom_point(data = landfalls[c(storms1,storms2),], 
             aes(x=BTLon, y = BT.Lat, color="1/2"), size=4,shape=16,fill=yellow)+
  geom_point(data = landfalls[c(storms3,storms4,storms5),], 
             aes(x=BTLon, y = BT.Lat, color="3/4/5"),size=4,shape=18,fill=verm)+
  geom_segment(aes(x = -88, y = 26, xend = -88, yend = 32), color="red")+
  geom_segment(aes(x = -78, y = 32, xend = -88, yend = 32), color="red")+
  labs(x = "Longitude", y="Latitude") + 
  theme(text = element_text(size=11)) +
  coord_fixed(xlim=c(-98, -70), ylim = c(24.5,42), ratio = 1)+
  scale_colour_manual(name="Intensity", values = c(skyblue,yellow,verm),
                      guide = guide_legend(override.aes = list(shape = c(15,16,18), 
                                                               color = c(skyblue,yellow,verm))))
  
g_byint

png(paste0("~/NAM-Model-Validation/png/g_byint_part.png"),width=1775, height=1000, res=350)
g_byint
dev.off()

ggsave(
  "png/g_byint_gg_part.png",
  g_byint,
  width = 7.75,
  height = 4.5,
  dpi = 350
)

#################################
# Figure 2: preprocessing steps #
#################################



##################################
# Figure 3: Plot NAM, ST4, error #
##################################

source("scripts/multiplot.R")

# run SqrtPrecipVariograms.R with:
# * write.pdf <- F
# * line 40, chg to: storm.dirs <- list.dirs("~/NAMandST4", recursive = F)[47] eg for Nate
# * pick subtractPWmean T or F

# Removed some of the greener ones to get more yellow/red
precipcolors <- c("#FFFFFF", "#E0F8E0", "#81F781",
                  "#2EFE2E", "#C8FE2E", "#FFFF00",
                  "#FACC2E", "#FFBF00", "#FF8000", "#FF4000","#FF0000")

g1= ggplot(aes(x=x,y=y,fill=value),data=ST4_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max),
                       name = expression(sqrt(mm))) + 
  labs(title =paste(paste0(str_to_title(storm_name),": Stage IV Data")),x = "Longitude", y="Latitude") + 
  coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                     max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
              ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                     max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)


g2= ggplot(aes(x=x,y=y,fill=value),data=NAM_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max),
                       name = expression(sqrt(mm))) +
  labs(title =paste0(str_to_title(storm_name),": NAM Forecast"), x = "Longitude", y="Latitude") + 
  coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                     max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
              ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                     max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)


g3= ggplot(aes(x=x,y=y,fill=value),data=error_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white",
                       name = expression(sqrt(mm))) +
  labs(title =paste(paste0(str_to_title(storm_name),": Stage IV - NAM")),
       x = "Longitude", y="Latitude")+
  coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                     max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
              ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                     max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)


png(paste0("~/NAM-Model-Validation/png/NAM_ST4_error_", storm_year, storm_name,".png"),
    width=3600, height=1050, res=350)
multiplot(g1, g2, g3, cols=3)
dev.off()

#######################
# Figure 3: No titles #
####################### 

g1= ggplot(aes(x=x,y=y,fill=value),data=ST4_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradientn(colors = precipcolors ,na.value = "white",
                       limits=c(0,precip.max), name = expression(sqrt(mm))) + 
  labs(x = "Longitude\n(a)", y="Latitude") + 
  theme(text = element_text(size=11)) +
  coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                     max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
              ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                     max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)

g2= ggplot(aes(x=x,y=y,fill=value),data=NAM_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradientn(colors = precipcolors ,na.value = "white",
                       limits=c(0,precip.max), name = expression(sqrt(mm))) +
  labs(x = "Longitude\n(b)", y="Latitude") + 
  theme(text = element_text(size=11)) +
  coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                     max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
              ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                     max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)

g3= ggplot(aes(x=x,y=y,fill=value),data=error_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",
                       na.value = "white", name = expression(sqrt(mm))) +
  labs(x = "Longitude\n(c)", y="Latitude")+
  theme(text = element_text(size=11)) +
  coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                     max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
              ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                     max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)

# grid.arrange(g3, bottom = textGrob("(a)", vjust = 0))

ggsave(
  "png/NAM_ST4_error_gg.png",
  multiplot(g1,g2,g3,cols=3),
  width = 12,
  height = 3.5,
  dpi = 350
)

png(paste0("~/NAM-Model-Validation/png/NAM_ST4_error_", storm_year, storm_name,"_notitle.png"),
    width=3600, height=1000, res=350)
multiplot(g1, g2, g3, cols=3)
dev.off()



########################
# Figure: error counts #
########################

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


##########################
# Figure 4: Plot PWmeans #
##########################

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

plot(PW_mean, main="Pointwise Mean")
PW_post_spdf <- as((PW_post_new), "SpatialPixelsDataFrame")
PW_post_df <- as.data.frame(PW_post_spdf)
colnames(PW_post_df) <- c("value", "x", "y")

plot(PW_post_sds, main="Pointwise SDs")
PW_sds_spdf <- as((PW_post_sds), "SpatialPixelsDataFrame")
PW_sds_df <- as.data.frame(PW_sds_spdf)
colnames(PW_sds_df) <- c("value", "x", "y")

plot(PW_post_new/PW_post_sds, main="Standardized Pointwise Mean")
PW_stdz_spdf <- as((PW_post_new/PW_post_sds), "SpatialPixelsDataFrame")
PW_stdz_df <- as.data.frame(PW_stdz_spdf)
colnames(PW_stdz_df) <- c("value", "x", "y")

g5= ggplot(aes(x=x,y=y,fill=value),data=PW_post_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white",
                       limit = range(c(PW_post_df$value, PW_sds_df$value, PW_stdz_df$value)),
                       name = expression(sqrt(mm))) +
    labs(x = "Longitude\n(a)", y="Latitude") + 
  theme(text = element_text(size=11))+  coord_fixed(xlim=extent(PW_post_new)[1:2]+c(1,-1), 
              ylim = extent(PW_post_new)[3:4]+c(1,-1), ratio = 1.33)

g5

g6= ggplot(aes(x=x,y=y,fill=value),data=PW_sds_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white", 
                       limit = range(c(PW_post_df$value, PW_sds_df$value, PW_stdz_df$value)),
                       name = expression(sqrt(mm))) +
    labs(x = "Longitude\n(b)", y="Latitude") + 
  theme(text = element_text(size=11))+  
  coord_fixed(xlim=extent(PW_post_new)[1:2]+c(1,-1), ylim = extent(PW_post_new)[3:4]+c(1,-1), ratio = 1.33)

g6

g7= ggplot(aes(x=x,y=y,fill=value),data=PW_stdz_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white", 
                       limit = range(c(PW_post_df$value, PW_sds_df$value, PW_stdz_df$value)),
                       name = expression(sqrt(mm))) +
    labs(x = "Longitude\n(c)", y="Latitude") + 
  theme(text = element_text(size=11))+
  coord_fixed(xlim=extent(PW_post_new)[1:2]+c(1,-1), ylim = extent(PW_post_new)[3:4]+c(1,-1), ratio = 1.33)

g7

multiplot(g5,g6,g7,cols = 3)

library(lemon)
grid_arrange_shared_legend(g5,g6,g7,position = "right")

ggsave(
  "png/PWmean_sds_stdz_gg.png",
  grid_arrange_shared_legend(g5,g6,g7,position = "right"),
  width = 10*1.5,
  height = 2.5*1.5,
  dpi = 350
)

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







###############################
# Figure 5: post vs MLE plots #
###############################

# from GibbsSamplerHurrRegrNA.R
load("~/NAM-Model-Validation/RData/Gibbs_sqrt.RData")

my_cols <- c(skyblue,yellow,verm)

# Lower and upper bound extra space to see end of distbns
up_b <- c(0,0.7); lw_b <- c(0,0.3)
smush <- c(0.02, 0.15)
## Plot MLE vs posterior estimates for the 47 storms
png("~/NAM-Model-Validation/png/NMV_compareMLEpost.png", width = 1400, height = 700)
par(mar=c(4.5,1.5,1.5,1.5), mfrow=c(1,2))
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
           y1 = 1, col = my_cols[factor(loc_int$loc)], lwd=1)
  points(x=emp_thetaMED[,i], y=rep(0, length(emp_thetaMED[,i])), col= my_cols[factor(loc_int$loc)])
  points(x=theta_hat[,i],    y=rep(1, length(emp_thetaMED[,i])), col= my_cols[factor(loc_int$loc)])
  for (j in 1:N) {
    #top row, theta hats
    xseq <- seq(theta_hat[j,i]-3*sqrt(solve(hessians[[j]])[i,i]),
                theta_hat[j,i]+3*sqrt(solve(hessians[[j]])[i,i]),
                length.out = 1000)
    lines(xseq,smush[i]*dnorm(xseq,
                              theta_hat[j,i],
                              sqrt(solve(hessians[[j]])[i,i]))+1, col=my_cols[factor(loc_int$loc)[j]],
          lty = 1, lwd=3)
    
    #bottom row, thetas from MCMC
    xseq <- seq(min(theta_burn[[j]][,i]),
                max(theta_burn[[j]][,i]), 
                length.out = 1000)
    lines(xseq, smush[i]*dnorm(xseq, 
                               mean(theta_burn[[j]][,i]), 
                               sd(theta_burn[[j]][,i])), col=my_cols[factor(loc_int$loc)[j]],
          lty = 1, lwd=3)
  }
  legend("topright",
         2, lwd=3,
         legend=levels(factor(loc_int$loc)),
         col=my_cols, lty=1, cex=2)
}
dev.off()



my_cols <- c(skyblue,yellow,verm)

# Lower and upper bound extra space to see end of distbns
up_b <- c(0,0.7); lw_b <- c(0,0.3)
smush <- c(0.02, 0.15)
## Plot MLE vs posterior estimates for the 47 storms
png("~/NAM-Model-Validation/png/NMV_compareMLEpost_newcap.png", width = 1400, height = 700)
par(mar=c(5.5,1.5,1.5,1.5), mfrow=c(1,2))
for (i in 1:P) {
  plot(0, 0, col = "white", ylab = "", 
       xlim=c(min(theta_hat[,i], emp_thetaMED[,i])-lw_b[i],
              max(theta_hat[,i], emp_thetaMED[,i])+up_b[i]),
       ylim=c(0,2), yaxt='n', cex = 3, cex.axis =2, cex.lab = 2, xlab="")
  if(i==1){  mtext("(a)", side=1, line=3, cex=2)} 
  if(i==2){  mtext("(b)", side=1, line=3, cex=2)} 
  # if(i==1){  mtext(expression(paste("(a): ",hat(theta)[1], " top and ", theta[1], " bottom")), side=1, line=4.5, cex=2)} 
  # if(i==2){  mtext(expression(paste("(b): ",hat(theta)[2], " top and ", theta[2], " bottom")), side=1, line=4.5, cex=2)}
       # xlab = bquote(theta[.(i)]~"bottom and"~hat(theta)[.(i)]~"top"))
  # mult_seg <- data.frame(x0 = c(0.7, 0.2, - 0.9, 0.5, - 0.2),    # Create data frame with line-values
  #                        y0 = c(0.5, 0.5, 0.6, - 0.3, 0.4),
  #                        x1 = c(0, 0.2, 0.4, - 0.3, - 0.6),
  #                        y1 = c(- 0.1, 0.3, - 0.6, - 0.8, 0.9))
  
  segments(x0 = emp_thetaMED[,i],                            # Draw multiple lines
           y0 = 0,
           x1 = theta_hat[,i],
           y1 = 1, col = my_cols[factor(loc_int$loc)], lwd=1)
  points(x=emp_thetaMED[,i], y=rep(0, length(emp_thetaMED[,i])), col= my_cols[factor(loc_int$loc)])
  points(x=theta_hat[,i],    y=rep(1, length(emp_thetaMED[,i])), col= my_cols[factor(loc_int$loc)])
  for (j in 1:N) {
    #top row, theta hats
    xseq <- seq(theta_hat[j,i]-3*sqrt(solve(hessians[[j]])[i,i]),
                theta_hat[j,i]+3*sqrt(solve(hessians[[j]])[i,i]),
                length.out = 1000)
    lines(xseq,smush[i]*dnorm(xseq,
                              theta_hat[j,i],
                              sqrt(solve(hessians[[j]])[i,i]))+1, col=my_cols[factor(loc_int$loc)[j]],
          lty = 1, lwd=3)
    
    #bottom row, thetas from MCMC
    xseq <- seq(min(theta_burn[[j]][,i]),
                max(theta_burn[[j]][,i]), 
                length.out = 1000)
    lines(xseq, smush[i]*dnorm(xseq, 
                               mean(theta_burn[[j]][,i]), 
                               sd(theta_burn[[j]][,i])), col=my_cols[factor(loc_int$loc)[j]],
          lty = 1, lwd=3)
  }
  legend("topright",
         2, lwd=3,
         legend=levels(factor(loc_int$loc)),
         col=my_cols, lty=1, cex=2)
}
dev.off()




# png("~/NAM-Model-Validation/png/NMV_compareMLEpost2.png", width = 1500, height = 500)
# par(mar=c(4.5,1.5,1.5,1.5), mfrow=c(1,3))
# for (i in 1:P) {
#   plot(0, 0, col = "white", ylab = "", 
#        xlim=c(min(theta_hat[,i], emp_thetaMED[,i])-lw_b[i],
#               max(theta_hat[,i], emp_thetaMED[,i])+up_b[i]),
#        ylim=c(0,2), yaxt='n', cex = 3, cex.axis =2, cex.lab = 2,
#        xlab = bquote(theta[.(i)]~"bottom and"~hat(theta)[.(i)]~"top"))
#   # mult_seg <- data.frame(x0 = c(0.7, 0.2, - 0.9, 0.5, - 0.2),    # Create data frame with line-values
#   #                        y0 = c(0.5, 0.5, 0.6, - 0.3, 0.4),
#   #                        x1 = c(0, 0.2, 0.4, - 0.3, - 0.6),
#   #                        y1 = c(- 0.1, 0.3, - 0.6, - 0.8, 0.9))
#   
#   segments(x0 = emp_thetaMED[,i],                            # Draw multiple lines
#            y0 = 0,
#            x1 = theta_hat[,i],
#            y1 = 1, col = my_cols[loc_int$loc])
#   points(x=emp_thetaMED[,i], y=rep(0, length(emp_thetaMED[,i])), col= loc_int$loc)
#   points(x=theta_hat[,i],    y=rep(1, length(emp_thetaMED[,i])), col= loc_int$loc)
#   for (j in 1:N) {
#     #top row, theta hats
#     xseq <- seq(theta_hat[j,i]-3*sqrt(solve(hessians[[j]])[i,i]),
#                 theta_hat[j,i]+3*sqrt(solve(hessians[[j]])[i,i]),
#                 length.out = 1000)
#     lines(xseq,smush[i]*dnorm(xseq,
#                               theta_hat[j,i],
#                               sqrt(solve(hessians[[j]])[i,i]))+1, col=loc_int$loc[j],
#           lty = as.integer(loc_int$loc)[j], lwd=3)
#     
#     #bottom row, thetas from MCMC
#     xseq <- seq(min(theta_burn[[j]][,i]),
#                 max(theta_burn[[j]][,i]), 
#                 length.out = 1000)
#     lines(xseq, smush[i]*dnorm(xseq, 
#                                mean(theta_burn[[j]][,i]), 
#                                sd(theta_burn[[j]][,i])), col=loc_int$loc[j],
#           lty = as.integer(loc_int$loc)[j], lwd=3)
#   }
#   legend("topright", 
#          2, lwd=3,
#          legend=c(levels(loc_int$loc)[unique(as.integer(loc_int$loc))]),
#          col=unique(as.integer(loc_int$loc)), lty=unique(as.integer(loc_int$loc)), cex=2)
# }
# dev.off()
# 
# png("~/NAM-Model-Validation/png/NMV_compareMLEpost3.png", width = 1400, height = 700)
# par(mar=c(4.5,1.5,1.5,1.5), mfrow=c(1,3))
# for (i in 1:P) {
#   plot(0, 0, col = "white", ylab = "", 
#        xlim=c(min(theta_hat[,i], emp_thetaMED[,i])-lw_b[i],
#               max(theta_hat[,i], emp_thetaMED[,i])+up_b[i]),
#        ylim=c(0,2), yaxt='n', cex = 3, cex.axis =2, cex.lab = 2,
#        xlab = bquote(theta[.(i)]~"bottom and"~hat(theta)[.(i)]~"top"))
#   
#   segments(x0 = emp_thetaMED[,i],                            # Draw multiple lines
#            y0 = 0,
#            x1 = theta_hat[,i],
#            y1 = 1, col = loc_int$loc, lwd=3)
#   points(x=emp_thetaMED[,i], y=rep(0, length(emp_thetaMED[,i])), col= loc_int$loc)
#   points(x=theta_hat[,i],    y=rep(1, length(emp_thetaMED[,i])), col= loc_int$loc)
#   for (j in 1:N) {
#     #top row, theta hats
#     xseq <- seq(theta_hat[j,i]-3*sqrt(solve(hessians[[j]])[i,i]),
#                 theta_hat[j,i]+3*sqrt(solve(hessians[[j]])[i,i]),
#                 length.out = 1000)
#     lines(xseq,smush[i]*dnorm(xseq,
#                               theta_hat[j,i],
#                               sqrt(solve(hessians[[j]])[i,i]))+1, col=loc_int$loc[j],
#           lwd=3) #lty = as.integer(loc_int$loc)[j], 
#     
#     #bottom row, thetas from MCMC
#     xseq <- seq(min(theta_burn[[j]][,i]),
#                 max(theta_burn[[j]][,i]), 
#                 length.out = 1000)
#     lines(xseq, smush[i]*dnorm(xseq, 
#                                mean(theta_burn[[j]][,i]), 
#                                sd(theta_burn[[j]][,i])), col=loc_int$loc[j],
#           lwd=3) #lty = as.integer(loc_int$loc)[j], 
#   }
#   legend("topright", 
#          2, lwd=3,
#          legend=c(levels(loc_int$loc)[unique(as.integer(loc_int$loc))]),
#          col=unique(as.integer(loc_int$loc)), lty=unique(as.integer(loc_int$loc)), cex=2)
# }
# dev.off()

#####################################
# Figure 6: Real vs Sim Error Field #
#####################################

# Run predictionNA_2018_SE.R with s <- 4 for Michael
# produces the NAM_r and ST4_r below, also needed for Figures 7 & 8

plot(NAM_r-ST4_r, main="Actual Error Field"); US(add=T, col="gray")
plot(sim_r, main="Simulated Error Field"); US(add=T, col="gray")

error_r_spdf <- as((NAM_r - ST4_r), "SpatialPixelsDataFrame")
error_r_df <- as.data.frame(error_r_spdf)
colnames(error_r_df) <- c("value", "x", "y")

sim_r_spdf <- as((sim_r), "SpatialPixelsDataFrame")
sim_r_df <- as.data.frame(sim_r_spdf)
colnames(sim_r_df) <- c("value", "x", "y")


g8= ggplot(aes(x=x,y=y,fill=value),data=error_r_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white", 
                       limit = range(c(sim_r_df$value, error_r_df$value)), 
                       name = expression(sqrt(mm))) +
    labs(x = "Longitude\n(a)", y="Latitude") + 
  theme(text = element_text(size=11))+
  coord_fixed(xlim=extent(sim_r)[1:2], ylim = extent(sim_r)[3:4], ratio = 1)

g9= ggplot(aes(x=x,y=y,fill=value),data=sim_r_df) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white", 
                       limit = range(c(sim_r_df$value, error_r_df$value)),
                       name = expression(sqrt(mm))) +
    labs(x = "Longitude\n(b)", y="Latitude") + 
  theme(text = element_text(size=11))+
  coord_fixed(xlim=extent(sim_r)[1:2], ylim = extent(sim_r)[3:4], ratio = 1)


library(lemon)
grid_arrange_shared_legend(g8,g9,position = "right")

ggsave(
  "png/real_vs_sim_errors.png",
  grid_arrange_shared_legend(g8,g9,position = "right"),
  width = 6*1.5,
  height = 2.5*1.5,
  dpi = 350
)

#############################
# Figure 7: 2" and prob map #
#############################

par(mfrow=c(1,3))
plot(rasterFromXYZ(cbind(coords, two_in_rain_f)), main="rain > 2\" in forecast", cex.main=1.5, legend=F, cex.axis=1.5)
plot(rasterFromXYZ(cbind(coords,prob_map_vals)), main = "prob map of rain > 2\"",
     cex.main=1.5, zlim=c(0,1), cex.axis=1.5)
# addRasterLegend(rasterFromXYZ(cbind(coords, prob_map_vals)),zlim=c(0,1),cex.axis = 1.4)
plot(rasterFromXYZ(cbind(coords, two_in_rain_o)), main="rain > 2\" in observed", cex.main=1.5, legend=F, cex.axis=1.5)

g10= ggplot(aes(x=x,y=y,fill=value),data=data.frame(cbind(x=coords[,1],y=coords[,2], value=two_in_rain_f))) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "white", high = "blue",na.value = "white",
                       limit = 0:1, name = "P(rain > 50.8mm)") +
  scale_color_hue(direction = -1) +
  theme(text = element_text(size=11))+
    labs(x = "Longitude\n(a)", y="Latitude") + 
  coord_fixed(xlim=extent(sim_r)[1:2], ylim = extent(sim_r)[3:4], ratio = 1)

g11= ggplot(aes(x=x,y=y,fill=value),data=data.frame(cbind(x=coords[,1],y=coords[,2], value=two_in_rain_o))) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "white", high = "blue",na.value = "white",
                       limit = 0:1, name = "P(rain > 50.8mm)") +
  theme(text = element_text(size=11))+
    labs(x = "Longitude\n(b)", y="Latitude") + 
  coord_fixed(xlim=extent(sim_r)[1:2], ylim = extent(sim_r)[3:4], ratio = 1)

g12= ggplot(aes(x=x,y=y,fill=value),data=data.frame(cbind(x=coords[,1],y=coords[,2], value=prob_map_vals))) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradient2(low = "white", high = "blue",na.value = "white",
                       limit = 0:1, name = "P(rain > 50.8mm)") +
  scale_color_hue(direction = -1, h.start=90) +
  theme(text = element_text(size=11))+
    labs(x = "Longitude\n(c)", y="Latitude") + 
  coord_fixed(xlim=extent(sim_r)[1:2], ylim = extent(sim_r)[3:4], ratio = 1)

g12

lemon::grid_arrange_shared_legend(g10,g11,g12,position = "right")

ggsave(
  "png/prob_map_precip.png",
  lemon::grid_arrange_shared_legend(g10,g11,g12,position = "right"),
  width = 9*1.5,
  height = 2.5*1.5,
  dpi = 350
)


###################################
# Figure 8: split up precip probs #
###################################

# par(mfrow=c(1,2))
# plot(rasterFromXYZ(cbind(coords, prob_map_vals)), main="all rain probs")
# hist(prob_map_vals, main="all rain probs")
png("png/prob_map_split_1x4.png", res = 350, width = 2400, height = 600)
par(mfrow=c(1,4), mar = c(3.5,2,1,0))
#, main="rain > 2\" in observed\n plotted by prob of > 2\"")
plot(rasterFromXYZ(cbind(coords, prob_map_vals*two_in_rain_oNA)), cex.axis=0.7, cex.lab=0.7)
par(mar=c(3.5,4,1,1))
#, main="prob of > 2\" associated with pts \nthat actually were over 2\"")
hist(prob_map_vals*two_in_rain_oNA, main=NULL, xlab="",ylab="", prob=T, cex.axis=0.7, cex.lab=0.7)
mtext("Density", side=2, line=2, cex=0.5)
mtext("P(rain > 50.8 mm)", side=1, line=2, cex=0.5)
par(mar=c(3.5,2,1,1))
#, main="rain < 2\" in observed\n plotted by prob of < 2\"")
plot(rasterFromXYZ(cbind(coords, prob_map_vals*two_in_rain_less2)), cex.axis=0.7, cex.lab=0.7)
par(mar=c(3.5,4,1,1))
#, main="prob of < 2\" associated with pts \nthat actually were under 2\"")
hist(prob_map_vals*two_in_rain_less2, main=NULL,xlab="",ylab="", prob=T, cex.axis=0.7, cex.lab=0.7)
mtext("Density", side=2, line=2, cex=0.5)
mtext("P(rain > 50.8 mm)", side=1, line=2, cex=0.5)
dev.off()

# 2x2 plot

png("png/prob_map_split.png", res = 350, width = 1600, height = 1600)
par(mfrow=c(2,2), mar = c(3,2,1,0))
plot(rasterFromXYZ(cbind(coords, prob_map_vals*two_in_rain_oNA)), cex.axis=0.7, cex.lab=0.7)
par(mar=c(3,4,1,1))
hist(prob_map_vals*two_in_rain_oNA, main=NULL, xlab="",ylab="", prob=T, cex.axis=0.7, cex.lab=0.7)
mtext("Density", side=2, line=2, cex=0.66)
mtext("P(rain > 50.8 mm)", side=1, line=2, cex=0.66)
par(mar = c(3,2,1,0))
plot(rasterFromXYZ(cbind(coords, prob_map_vals*two_in_rain_less2)), cex.axis=0.7, cex.lab=0.7)
par(mar=c(3,4,1,1))
hist(prob_map_vals*two_in_rain_less2, main=NULL,xlab="",ylab="", prob=T, cex.axis=0.7, cex.lab=0.7)
mtext("Density", side=2, line=2, cex=0.66)
mtext("P(rain > 50.8 mm)", side=1, line=2, cex=0.66)
dev.off()

nlon <- length(unique(coords[,1]))
nlat <- length(unique(coords[,2]))

r <- rasterFromXYZ(cbind(coords, two_in_rain_o))
m <- matrix(data = values(r), nrow = nlon, ncol = nlat)
image(m[,nrow(m):1])

count <- 0
border <- matrix(0, nlon, nlat)
for (i in 2:(nlon-1)) {
  for (j in 2:(nlat-1)) {
    if(is.na(m[i,j]) | m[i,j] == 0){next}
    if(sum(m[i+1,j] == 0, m[i,j+1] == 0, 
           m[i-1,j] == 0, m[i,j-1] == 0,
           m[i+1,j+1] == 0, m[i-1,j+1] == 0, 
           m[i-1,j-1] == 0, m[i+1,j-1] == 0,
           na.rm = T) > 1){
      count <- count + 1
      border[i,j] <- 2
    }
  }
}

image(m[,nrow(m):1])
image(border[,nrow(border):1])

precip.max <- max(ST4_pred$value)
storm

g13= ggplot(aes(x=x,y=y,fill=value),data=NAM_pred) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max+max(ub_rain)),
                       name = expression(sqrt(mm))) + 
  labs(x = "Longitude\n(a)", y="Latitude") + #title =paste(paste0(str_to_title(name),": NAM Forecast")),
  coord_fixed(xlim=range(NAM_pred$x), ylim=range(NAM_pred$y), ratio = 1)

g14= ggplot(aes(x=x,y=y,fill=value),data=NAM_pred+cbind(0,0,ub_rain)) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max+max(ub_rain)),
                       name = expression(sqrt(mm))) + 
  labs(x = "Longitude\n(b)", y="Latitude") + #title =paste(paste0(str_to_title(name),": NAM + 95% UB")),
  coord_fixed(xlim=range(NAM_pred$x), ylim=range(NAM_pred$y), ratio = 1)

g15= ggplot(aes(x=x,y=y,fill=value),data=ST4_pred) + 
  geom_tile() + theme_classic() + 
  geom_polygon(data=subset(map_data("state"), region %in% regions), 
               aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max+max(ub_rain)),
                       name = expression(sqrt(mm))) + 
  labs(x = "Longitude\n(c)", y="Latitude") +#title =paste(paste0(str_to_title(name),": Stage IV")),
  coord_fixed(xlim=range(ST4_pred$x), ylim=range(ST4_pred$y), ratio = 1)

grid_arrange_shared_legend(g13,g14,g15,position = "right")


ggsave(
  "png/pred_NAM_us_ST4.png",
  grid_arrange_shared_legend(g13,g14,g15,position = "right"),
  width = 8.5*1.5,
  height = 2.5*1.5,
  dpi = 350
)
