library(raster) #raster, extent
library(fields) #US()
library(ggplot2)

black <- rgb(0,0,0,1)
orange <- rgb(230/255,159/255,0,1)
skyblue <- rgb(86/255,180/255,233/255,1)
bluegreen <- rgb(0,158/255,115/255,1)
yellow <- rgb(240/255,228/255,66/255,1)
blue <- rgb(0,114/255,178/255,1)
verm <- rgb(213/255, 94/255, 0,1)
redpurp <- rgb(204/255, 121/255, 167/255, 1)
plot(1:10, col=c(black,yellow,blue, redpurp,verm))

if(!dir.exists("png/")) dir.create("png/", recursive = T)

###################################
# Figure 1: landfalls by location #
###################################

regions<- c("alabama", "alaska", "arizona", "arkansas", "california", "colorado", 
            "connecticut", "delaware", "district of columbia", "florida", "georgia", 
            "hawaii", "idaho", "illinois", "indiana", "iowa", "kansas", "kentucky", 
            "louisiana", "maine", "maryland", "massachusetts", "michigan", "minnesota",
            "mississippi", "missouri", "montana", "nebraska", "nevada", "new hampshire", 
            "new jersey", "new mexico", "new york", "north carolina", "north dakota", 
            "ohio", "oklahoma", "oregon", "pennsylvania", "rhode island", "south carolina", 
            "south dakota", "tennessee", "texas", "utah", "vermont", "virginia", 
            "washington", "west virginia", "wisconsin", "wyoming")

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
  geom_segment(aes(x = -88, y = 26, xend = -88, yend = 32), lty=2)+
  geom_segment(aes(x = -78, y = 32, xend = -84, yend = 32), lty=2)+
  labs(x = "Longitude", y="Latitude") + 
  theme(text = element_text(size=11)) +
  coord_fixed(xlim=c(-98, -70), ylim = c(24.5,42), ratio = 1)+
  scale_colour_manual(name="Intensity", values = c(skyblue,yellow,verm),
                      guide = guide_legend(override.aes = list(shape = c(15,16,18), 
                                                               color = c(skyblue,yellow,verm))))
  
g_byint

png(paste0("png/g_byint_part.png"),width=1775, height=1000, res=350)
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


png(paste0("png/NAM_ST4_error_", storm_year, storm_name,".png"),
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

png(paste0("png/NAM_ST4_error_", storm_year, storm_name,"_notitle.png"),
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

png("png/error_counts.png", height = 500, width = 700)
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



###############################
# Figure 5: post vs MLE plots #
###############################

# from GibbsSamplerHurrRegrNA.R
load("RData/Gibbs_sqrt.RData")

my_cols <- c(skyblue,yellow,verm)

# Lower and upper bound extra space to see end of distbns
up_b <- c(0.07,1.45); lw_b <- c(0,0.25)
smush <- c(0.02, 0.15)
## Plot MLE vs posterior estimates for the 47 storms
png("png/NMV_compareMLEpost.png", width = 1400, height = 700)
par(mar=c(4.5,1.5,1.5,1.5), mfrow=c(1,2))
for (i in 1:P) {
  plot(0, 0, col = "white", ylab = "", 
       xlim=c(min(theta_hat[,i], emp_thetaMED[,i])-lw_b[i],
              max(theta_hat[,i], emp_thetaMED[,i])+up_b[i]),
       ylim=c(0,2), yaxt='n', cex = 3, cex.axis =2, cex.lab = 2,
       xlab = "")#bquote(theta[.(i)]~"bottom and"~hat(theta)[.(i)]~"top"))
  if(i==1){  mtext("(a)", side=1, line=3, cex=2)} 
  if(i==2){  mtext("(b)", side=1, line=3, cex=2)} 
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

########################################
# Figure 6, for every prediction storm #
########################################
for (ste in 1:6) {
  load(paste0("RData/prediction/prediction",ste,"nopw.RData"))
  precip.max <- max(ST4_pred$value, NAM_pred$value)
  
  g13= ggplot(aes(x=x,y=y,fill=value),data=NAM_pred) + 
    geom_tile() + theme_classic() + 
    geom_polygon(data=subset(map_data("state"), region %in% regions), 
                 aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
    scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max+max(ub_rain)),
                         name = expression(sqrt(mm))) + 
    labs(x = "Longitude\n(a)", y="Latitude") + #title =paste(paste0(str_to_title(name),": NAM Forecast")),
    coord_fixed(xlim=range(NAM_pred$x), ylim=range(NAM_pred$y), ratio = 1)
  
  g14= ggplot(aes(x=x,y=y,fill=value),data=NAM_pred[,-1]+cbind(ub_rain,0,0)) + 
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
    paste0("png/pred_NAM_us_ST4_storm",ste,".png"),
    grid_arrange_shared_legend(g13,g14,g15,position = "right"),
    width = 8.5*1.5,
    height = 2.5*1.5,
    dpi = 350
  )

}
#################################
# Table 1: Coverage percentages #
#################################

LM1_files <- grep(list.files(path="csv/prediction", full.names = T), pattern='LM', invert=TRUE, value=TRUE)
pred_res <- do.call(rbind, lapply(LM1_files, read.csv))
LM1_res <- 1-colMeans(pred_res)

pred_files <- (list.files("csv/prediction", pattern = "LM2", full.names = T))
pred_res <- do.call(rbind, lapply(pred_files, read.csv))
LM2_res <- 1-colMeans(pred_res)

pred_files <- (list.files("csv/prediction", pattern = "LM3", full.names = T))
pred_res <- do.call(rbind, lapply(pred_files, read.csv))
LM3_res <- 1-colMeans(pred_res)

rbind(LM1_res, LM2_res, LM3_res)[,3:4]

####################################
# Figure for hydrology application #
####################################

# run  FL_basin_mm.R first, with storm ste=4 and basin=4
png("png/hydro_app_FL_Michael_basin1.png", res=350, height=1000, width=3000)
par(mfrow=c(1,3))
par(mar=c(4.5,5,1,2))
plot(NA, xlim=c(-85.5,-83.5), ylim=c(29,32), asp=1,xlab= "Longitude", ylab="Latitude", cex.axis=1.3, cex.lab=1.3)
plot(SC_mask, add=T, legend = F); US(add=T, res=1)
par(mar=c(4.5,5,1,4))
plot(NA, xlim=c(-85.5,-83.5), ylim=c(29,32), asp=1,xlab= "Longitude", ylab="Latitude", cex.axis=1.3,cex.lab=1.3)
plot(rasterFromXYZ(rasterToPoints(mask_SC * NAM_r)),add=T, col=precipcolors, legend.width=1.5,
     legend.args = list(text = expression(sqrt(mm)), side = 3, 
                        font = 2, line = 1, cex = 0.8, width = 2))
US(add=T, res=1)
par(mar=c(4.5,5,1,2))
m <- 2
my_dens <- density(sum_sq_rains[,m], n = 2048)
pred_val_dens <- my_dens$y[ which(abs(my_dens$x - ST4_l_sq) == min(abs(my_dens$x - ST4_l_sq)))]
plot(my_dens$x, my_dens$y, type="l", xlab="mm", ylab="Density",
     xlim = c(0, 25000), cex.axis=1.3,cex.lab=1.3)
abline(v=NAM_l_sq, col="green")
abline(v=ST4_l_sq, col="blue", lty=2)
dev.off()
