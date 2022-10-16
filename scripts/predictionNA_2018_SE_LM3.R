##  Prediction for the future storms (validate with 2018 and 2019 method)
##  Establish regressors for an incoming storm (eg: $\boldsymbol{x}_p^T = [1,1,0]$ for a Florida storm)
##  Use x_p to generate theta_i^*=Bx_p + omega_i for i in 1:1000
##    From these theta_i^*, generate 1000 corresponding error fields
##  Add these generated error fields to the available NAM forecast 
##  Check to see if pointwise 95% percentile is greater than $\approx$ 95\% of pixels for the corresponding ST4. 
##  Compute probability maps that rain totals will surpass a particular threshold, etc.

library(LaplacesDemon)
suppressMessages(library(geoR))
suppressMessages(library(raster))
suppressMessages(library(RandomFields))
library(stringr)

if(!dir.exists("RData/prediction/")){
  dir.create("RData/prediction/", recursive = T)
}
if(!dir.exists("pdf/prediction/")){
  dir.create("pdf/prediction/", recursive = T)
}
if(!dir.exists("csv/prediction/")){
  dir.create("csv/prediction/", recursive = T)
}

# number of synthetic precipitation fields
Ngen <- 1000

# storm to evaluate
ste <- 4 #NAM_pred, ST4_pred, and x_pred, name of PDF change based on s
#change pwmean, sum(post)_cov_mtx, csv, load call below

# subtract the pointwise mean?
subtractPWmean <- F

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

PWstamp <- ifelse(subtractPWmean, "subtractpw", "nopw")
before <- Sys.time()

load(paste0("RData/Gibbs_sqrt_LM3",
            if(subtractPWmean){"_subtractPWmean"},".RData"))

path <- "csv/prediction_sqrt"
pred_dirs <- list.dirs(path, recursive = F, full.names = F)

# for (ste in 1:length(pred_dirs)) {

pred_dir <- pred_dirs[ste]

year <- substr(pred_dir, 1, 4)
name <- substr(pred_dir, 5, nchar(pred_dir))
print(paste0(year,name))

NAM_pred <- read.csv(paste0("csv/prediction_sqrt/",
                            year, name, "/", year, name,"_NAMdf.csv"))
ST4_pred <- read.csv(paste0("csv/prediction_sqrt/",
                            year, name, "/", year, name,"_ST4df.csv"))


## Uncertainty corresponding to pointwise mean
# PWM1_df <- read.csv("csv/PWM1_df.csv", row.names = 1)
# load("RData/sum_cov_mtx.RData") #post
# 
# ind <- c() 
# for (i in 1:nrow(PWM1_df)) {
#     if(length(which(abs(PWM1_df$x - NAM_pred$x[i]) < 1e-3 & abs(PWM1_df$y - NAM_pred$y[i]) < 1e-3))>0){
#       ind[i] <- which(abs(PWM1_df$x - NAM_pred$x[i]) < 1e-3 & abs(PWM1_df$y - NAM_pred$y[i]) < 1e-3)
#       next
#     }
# }
# 
# this_PWM <- PWM1_df[ind,]
# this_post_cov <- sum_cov_mtx[ind,ind] #post

# Alberto 2018: GULF storm
# First is always 1, ATL baseline; next are indicators for FL and GULF respectively
# 2018: Alberto FL, Florence ATL, Gordon GULF, Michael FL
# 2019: Barry GULF, Dorian ATL
# if(name %in% c("alberto", "michael"))  x_pred <- c(1,1,0)
# if(name %in% c("florence", "dorian"))  x_pred <- c(1,0,0)
# if(name %in% c("gordon", "barry"))     x_pred <- c(1,0,1)

pdf(paste0("pdf/prediction/prediction_sqrt_",name,
           year,"_GIS_GHiG_NA_flatPWmean_",PWstamp, Ngen,"_LM3.pdf"))

mask <- raster("lsmask.nc")
mask[mask==-1]  <- NA
extent(mask)[1] <- extent(mask)[1]-360
extent(mask)[2] <- extent(mask)[2]-360
mask.regrid <- raster::resample(mask, projectRaster(raster(
  "nam_218_20170826_0000_012.grb2"),
  crs = "+proj=longlat +datum=WGS84"), method='ngb') 

# run GibbsSamplerHurrRegr first for B and as.square

theta_pred <- mu_burn[10*(1:Ngen),]

# (neg_rows <- which(theta_pred < 0))
# for (i in 1:length(neg_rows)) {
#   col2chg <- 1
#   count <- Ngen
#   # print(c(col2chg,count))
#   while(neg_rows[i] > count){
#     count <- count + Ngen
#     col2chg <- col2chg + 1
#     # print(c(count, col2chg))
#   }
#   theta_pred[neg_rows[i]] <- median(theta_pred[,col2chg])
# }
# (neg_rows <- which(theta_pred < 0))

NAM_r <- rasterFromXYZ(NAM_pred[,c(3,4,2)])
ST4_r <- rasterFromXYZ(ST4_pred[,c(3,4,2)]) #row order messed up from logprecipvariogram write.csv

## Find the extra points from the files and remove them
## When NAM and ST4 have different amounts of pixels

## if NAM has more rows...
if(nrow(rasterToPoints(NAM_r)) > nrow(rasterToPoints(ST4_r))){
  ST40 <- ST4_r
  values(ST40)[!is.na(values(ST40))] <- 0
  # plot(ST40)
  NAM_r <- NAM_r - ST40
  NAM_pred <- as.data.frame(rasterToPoints(NAM_r))
  colnames(NAM_pred) <- c("x","y","value")
}

## if ST4 has more rows
if(nrow(rasterToPoints(NAM_r)) < nrow(rasterToPoints(ST4_r))){
  NAM0 <- NAM_r
  values(NAM0)[!is.na(values(NAM0))] <- 0
  # plot(NAM0)
  ST4_r <- ST4_r - NAM0
  ST4_pred <- as.data.frame(rasterToPoints(ST4_r))
  colnames(ST4_pred) <- c("x","y","value")
}

par(mfrow=c(1,2))
plot(NAM_r, main = "NAM forecast")
plot(ST4_r, main = "ST4 forecast")

if(subtractPWmean){
  PW_mean <- raster("error_rasters_summary_sqrt/PW_mean.grd")*mask.regrid
  NAM_r <- NAM_r - PW_mean
  NAM_pred <- as.data.frame(rasterToPoints(NAM_r))
  colnames(NAM_pred) <- c("x","y","value")
}

coords <- cbind(NAM_pred$x, NAM_pred$y)
coordsST4 <- cbind(ST4_pred$x, ST4_pred$y)

simvals <- simPWM <- matrix(NA, nrow = nrow(NAM_pred), ncol = Ngen) 
RFoptions(spConform=F) #for faster simulation
for (g in 1:Ngen) {
  if(g %% 50 == 0) cat(g,"\n")
  
  # Uncertainty in PWmean
  # simPWM[,g] <- t(chol(this_post_cov)) %*% rnorm(nrow(this_PWM))
  
  
  simvals[,g] <- RFsimulate(model = RMexp(var = exp(theta_pred[g,2]), 
                                          scale = exp(theta_pred[g,2]-theta_pred[g,1])), 
                            # err.model = RMnugget(var = 0),
                            x = cbind(NAM_pred$x, NAM_pred$y))
}

# example error field
sim <- cbind(coords, simvals[,1])
sim_r <-rasterFromXYZ(sim)

plot(NAM_r-ST4_r, main="Actual Error Field")
plot(sim_r, main="Simulated Error Field")

# calculate 95%ile for each point
ub_rain      <- apply(simvals, 1, function(x) quantile(x, .95))
# ub_rain_PW   <- apply(simvals + simPWM, 1, function(x) quantile(x, .95))
ub_rain99    <- apply(simvals, 1, function(x) quantile(x, .99))
# ub_rain99_PW <- apply(simvals + simPWM, 1, function(x) quantile(x, .99))
ub_rain100   <- apply(simvals, 1, max)
# ub_rain100_PW<- apply(simvals + simPWM, 1, max)
plot(rasterFromXYZ(cbind(coords, ub_rain+NAM_pred$value)), main = "NAM + 95% UB on EFs")

# compare NAM with ST4, and NAM+UB with ST4
off_base <- ifelse(NAM_pred$value < ST4_pred$value, 1, 0)
mean(off_base)
plot(rasterFromXYZ(cbind(coords,off_base)), main= bquote("NAM < ST4"~.(round(mean(off_base),4))))

par(mfrow=c(2,2))
off_est <- ifelse(ub_rain+NAM_pred$value < ST4_pred$value, 1, 0)
mean(off_est)
plot(rasterFromXYZ(cbind(coords,off_est)), main=bquote("NAM + 95% UB < ST4"~.(round(mean(off_est),4))))

# off_est_PW <- ifelse(ub_rain_PW+NAM_pred$value < ST4_pred$value, 1, 0)
# mean(off_est_PW)
# plot(rasterFromXYZ(cbind(coords,off_est_PW)), main=bquote("NAM + 95% UB < ST4"~.(round(mean(off_est_PW),4))))

off_est99 <- ifelse(ub_rain99+NAM_pred$value < ST4_pred$value, 1, 0)
mean(off_est99)
plot(rasterFromXYZ(cbind(coords,off_est99)), main=bquote("NAM + 99% UB < ST4"~.(round(mean(off_est99),4))))

# off_est99_PW <- ifelse(ub_rain99_PW+NAM_pred$value < ST4_pred$value, 1, 0)
# mean(off_est99_PW)
# plot(rasterFromXYZ(cbind(coords,off_est99_PW)), main=bquote("NAM + 99% UB < ST4"~.(round(mean(off_est99_PW),4))))

off_est100 <- ifelse(ub_rain100+NAM_pred$value < ST4_pred$value, 1, 0)
mean(off_est100)
plot(rasterFromXYZ(cbind(coords,off_est100)), main=bquote("NAM + 100% UB < ST4"~.(round(mean(off_est100),4))))

# off_est100_PW <- ifelse(ub_rain100_PW+NAM_pred$value < ST4_pred$value, 1, 0)
# mean(off_est100_PW)
# plot(rasterFromXYZ(cbind(coords,off_est100_PW)), main=bquote("NAM + 100% UB < ST4"~.(round(mean(off_est100_PW),4))))

ests <- cbind(mean(off_base), mean(off_est), mean(off_est99), mean(off_est100))
# ests_PW<-cbind(mean(off_base), mean(off_est_PW), mean(off_est99_PW), mean(off_est100_PW))

# write.csv(rbind(ests,ests_PW), paste0(file = "csv/prediction/", pred_dir,"_sqrt.csv"))
write.csv(ests, paste0(file = "csv/prediction/", pred_dir,"_sqrt_LM3_",
                       PWstamp,".csv"))

#2in = 50.8mm
sim_prob <- matrix(NA, nrow = nrow(simvals), ncol = Ngen)
for (i in 1:Ngen) {
  sim_prob[,i] <- ifelse(simvals[,i]+NAM_pred$value > sqrt(50.8), 1, 0)
}

two_in_rain_f <- ifelse(NAM_pred$value > sqrt(50.8), 1, 0)
two_in_rain_o <- ifelse(ST4_pred$value > sqrt(50.8), 1, 0)
two_in_rain_oNA <- ifelse(ST4_pred$value > sqrt(50.8), 1, NA)
two_in_rain_less2 <- ifelse(ST4_pred$value < sqrt(50.8), 1, NA)
prob_map_vals <- apply(sim_prob, 1, mean)

par(mfrow=c(1,3))
plot(rasterFromXYZ(cbind(coords, two_in_rain_f)), main="rain > 2\" in forecast", cex.main=1.5, legend=F, cex.axis=1.5)
plot(rasterFromXYZ(cbind(coords,prob_map_vals)), main = "prob map of rain > 2\"",
     cex.main=1.5, zlim=c(0,1), cex.axis=1.5)
# addRasterLegend(rasterFromXYZ(cbind(coords, prob_map_vals)),zlim=c(0,1),cex.axis = 1.4)
plot(rasterFromXYZ(cbind(coords, two_in_rain_o)), main="rain > 2\" in observed", cex.main=1.5, legend=F, cex.axis=1.5)

par(mfrow=c(1,2))
plot(rasterFromXYZ(cbind(coords, prob_map_vals)), main="all rain probs")
hist(prob_map_vals, main="all rain probs")

par(mfrow=c(2,2))
plot(rasterFromXYZ(cbind(coords, prob_map_vals*two_in_rain_oNA)))#, main="rain > 2\" in observed\n plotted by prob of > 2\"")
hist(prob_map_vals*two_in_rain_oNA, main=NULL, xlab=NULL)#, main="prob of > 2\" associated with pts \nthat actually were over 2\"")

# par(mfrow=c(1,2))
plot(rasterFromXYZ(cbind(coords, prob_map_vals*two_in_rain_less2)))#, main="rain < 2\" in observed\n plotted by prob of < 2\"")
hist(prob_map_vals*two_in_rain_less2, main=NULL,xlab=NULL)#, main="prob of < 2\" associated with pts \nthat actually were under 2\"")

# dev.off()

after <- Sys.time()
after - before

## Dave plots for water folks, from 4/28/20 notes
# pdf("pdf/prediction/NAM_ST4_with_error_sims.pdf")

min_p <- min(values(NAM_r), values(ST4_r), values(NAM_r), na.rm = T)
max_p <- max(values(NAM_r), values(ST4_r), values(NAM_r), na.rm = T)
for (i in (1:3)+1) {
  sim <- cbind(coords, simvals[,i])
  sim_r <-rasterFromXYZ(sim)
  min_p <- min(min_p, values(NAM_r+sim_r), na.rm = T)
  max_p <- max(max_p, values(NAM_r+sim_r), na.rm = T)
}
par(mfrow=c(2,3))
plot(NAM_r, zlim=c(min_p,max_p), main="NAM forecast"); #US(add=T, col="lightgray")
plot(ST4_r, zlim=c(min_p,max_p), main="ST4 observed"); #US(add=T, col="lightgray")
plot(NAM_r, zlim=c(min_p,max_p), main="NAM bias adj"); #US(add=T, col="lightgray")
for (i in (1:3)+1) {
  sim <- cbind(coords, simvals[,i])
  sim_r <-rasterFromXYZ(sim)
  plot(NAM_r+sim_r, main=paste("NAM bias adj + error field",i),zlim=c(min_p,max_p))
  #US(add=T, col="lightgray")
}

#plot with c(0,7 bounds)
par(mfrow=c(2,3))
plot(NAM_r, zlim=c(0,20), main="NAM forecast"); #US(add=T, col="lightgray")
plot(ST4_r, zlim=c(0,20), main="ST4 observed"); #US(add=T, col="lightgray")
plot(NAM_r, zlim=c(0,20), main="NAM bias adj"); #US(add=T, col="lightgray")
for (i in (1:3)+1) {
  sim <- cbind(coords, simvals[,i])
  sim_r <-rasterFromXYZ(sim)
  temp <- NAM_r + sim_r
  values(temp)[values(temp)>=20] = 20
  plot(temp, main=paste("NAM bias adj + error field",i),zlim=c(0,20))
  #US(add=T, col="lightgray")
}

# rm(list=setdiff(ls(), c("simvals", "NAM_pred", "ST4_pred", "s")))
save.image(paste0("RData/prediction/prediction",ste,"_LM3_",PWstamp,".RData"))
dev.off()

# }

rgs <- apply(simvals, 2, function(x) diff(range(x)))
sds <- apply(simvals, 2, function(x) sd(x))
hist(rgs)
hist(sds)

rgs <- apply(simvals, 1, function(x) diff(range(x)))
sds <- apply(simvals, 1, function(x) sd(x))
hist(rgs)
hist(sds)
