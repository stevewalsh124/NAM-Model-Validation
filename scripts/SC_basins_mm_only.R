# South Carolina river basins from the SC website: 
# https://sc-department-of-health-and-environmental-control-gis-sc-dhec.hub.arcgis.com/datasets/sc-major-river-basins
# PERIMETER	    HUC	      EXTHUC	  BASIN	              FEATURE
# 2002276.84709	03130011	99000000	APALACHICOLA RIVER	STREAM

# Read SHAPEFILE.shp from the current working directory (".")
# dsn="directory where the shapefile, projection file, etc are located" 
# layer="name of the file without .shp extention"

require(rgdal)
require(sf)
require(raster)
library(fields)
library(pracma)
library(scoringRules)
library(plgp)
library(LaplacesDemon)

make.pdf <- T

#storm(s) to eval
ste <- 3

# number of models compared
nM <- 9
# Ngen <- 1000 is loaded in the prediction RData

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

# Best states for pred storms 1:6 are AL, NC, AL/MS, FL/GA/AL, LA, NC/VA
# pick state based on the storm
if(ste==1) state <- "AL"
if(ste==2) state <- "NC"
if(ste==3) state <- "AL"#"MS"#
if(ste==4) state <- "FL"#"AL"#
if(ste==5) state <- "LA"
if(ste==6) state <- "NC"

if(make.pdf) pdf(paste0("~/NAM-Model-Validation/pdf/basin/",state,"_basins_storm",ste,".pdf"))

if(state == "AL") basins <- 1:11
if(state == "MS") basins <- 1:20
if(state == "NC") basins <- (1:18)[-14] #14 is NA
if(state == "NC" & ste == 6) basins <- c(1:10,12:13) #14 is NA, 11 and 15:18 don't intersect
if(state == "SC") basins <- 1:8 
if(state == "LA") basins <- 1:12
if(state == "FL") basins <- 1:8

log_scores <- CRPSs <- matrix(NA, max(basins), nM)       # for scoringRules version (condl on theta)
log_scores_marf <- log_scores_sum_marf <- matrix(NA, max(basins), nM)  # average out thetas, fn below
bas_locs <- coastals <- c() # keep track of the landfall region for each basin, coastal or not

log_score <- function(thet, dist_mtx, obsvn, mean_fn, G=(length(thet)/2)){
  thet <- matrix(thet, length(thet)/2, 2)
  log_pred_dens <- CRPSs <- c()
  for(g in 1:G){
    # if(g %% 500 ==0) print(g)
    sig2 = exp(thet[g,2])
    phi = exp(thet[g,2]-thet[g,1])
    my_cov_mtx <- sig2 * exp(-0.5 * my_dist_mtx / phi)
    log_pred_dens[g] <-  dmvnorm(x = obsvn, mean = mean_fn, 
                                 sigma = my_cov_mtx, checkSymmetry = F, log = T)
    # CRPSs[g] <- crps_norm(y = obsvn, mean = mean_fn, sd = chol(my_cov_mtx))
    if(length(thet)==2) break
  }
  return(log_pred_dens)
}

log_score_straw <- function(thet, obsvn, mean_fn, G=length(thet)){
  thet <- matrix(thet, length(thet), 1)
  log_pred_dens <- CRPSs <- c()
  for(g in 1:G){
    # if(g %% 500 ==0) print(g)
    sig2 = exp(thet[g])
    my_cov_mtx <- diag(sig2, length(mean_fn))
    log_pred_dens[g] <-  dmvnorm(x = obsvn, mean = mean_fn, 
                                 sigma = my_cov_mtx, checkSymmetry = F, log = T)
    # CRPSs[g] <- crps_norm(y = obsvn, mean = mean_fn, sd = chol(my_cov_mtx))
    if(length(thet)==1) break
  }
  return(log_pred_dens)
}

meanOrig <- function(X){log(mean(exp(X - max(X)))) + max(X)}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

for (basin in basins) {
  print(paste0("basin ",basin))
  SC_mask <- raster(paste0("~/NAM-Model-Validation/basin/",state,"mask",basin,".grd"))
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",ste,"nopw.RData"))
  
  if(state != "FL") mask_SC <- projectRaster(SC_mask, crs = "+proj=longlat +datum=WGS84", method = "ngb")
  if(state == "FL") mask_SC <- SC_mask
  
  par(mfrow=c(1,1))
  mask_SC <- raster::resample(mask_SC, rasterFromXYZ(cbind(coords, simvals[,1])), method="ngb")
  print(sum(!is.na(values(mask_SC))))
  if(sum(!is.na(values(mask_SC))) < 30) next
  plot(mask_SC, main = paste("basin #",basin, bas_locs[basin],coastals[basin]), 
       xlim = extent(mask_SC)[1:2] + c(-1,1), 
       ylim = extent(mask_SC)[3:4] + c(-1,1))
  US(add=T)
  
  # get basin-specific dist matrix
  mask_locs <- rasterToPoints(NAM_r * mask_SC)[,1:2]
  my_dist_mtx <- sqrt(plgp::distance(X1 = mask_locs))
  
  # determine location of each basin (by majority vote)
  basin_df <- data.frame(mask_locs)
  bas_loc <- ifelse(basin_df$x < -88, "G", ifelse(basin_df$y > 32, "A", "F"))
  # print(table(bas_loc))
  bas_locs[basin] <- Mode(bas_loc)
  
  # determine coastal or not for the basin
  if(basin == basins[1]){
    values(mask) <- ifelse(is.na(values(mask)),0,1)
    masky <- data.frame(rasterToPoints(mask)) # CONUS mask
    masky <- masky[masky$y < 14/27*masky$x+87,] # remove top left of CONUS
    masky <- masky[masky$y > 25,] # remove Cuba
    masky <- masky[masky$Land.Sea.Mask == 0,] # only look at the ocean values
  }
  coast_dist_mtx <- sqrt(plgp::distance(X1 = mask_locs, X2=masky[,c("x","y")]))
  print(ifelse(min(coast_dist_mtx) < 1 ,"C","N"))
  coastals[basin] <- ifelse(min(coast_dist_mtx) < 1,"C","N") # pick if coastal or not
  
  Ngen <- ncol(simvals)
  sum_sq_rains <- matrix(NA, Ngen, nM)
  
  for (storm in ste) {
    for (mm in 1:nM) {
      print(paste0("storm ", storm, ", model ", mm))
      if(mm==1) load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",storm,"subtractpw.RData"))
      if(mm==2) load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",storm,"_LM2_subtractpw.RData"))
      if(mm==3) load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",storm,"_LM3_subtractpw.RData"))
      if(mm==4) load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",storm,"_LM4_subtractpw.RData"))
      if(mm==5) load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",storm,"nopw.RData"))
      if(mm==6) load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",storm,"_LM2_nopw.RData"))
      if(mm==7) load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",storm,"_LM3_nopw.RData"))
      if(mm==8) load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",storm,"_LM4_nopw.RData"))
      if(mm==9) load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",storm,"nosp_nopw"))
      NAM_pred <- data.frame(NAM_pred$x, NAM_pred$y, NAM_pred$value)
      ST4_pred <- cbind(ST4_pred$x, ST4_pred$y, ST4_pred$value)
      colnames(NAM_pred) <- c("x","y","value")
      sim1 <- rasterFromXYZ(cbind(coords, NAM_pred[,3]+simvals[,1])) * mask_SC
      
      if(mm==1){
        par(mfrow=c(1,1), mar = c(5,4,4,2)+0.1)
        plot(sim1, main=paste0("storm ", storm, ": ", name,", basin ", basin)); US(add=T)
      }
      sum(unique(values(sim1)),na.rm = T)
      
      
      
      for (i in 1:Ngen) {
        # if((i %% 5000)==0) print(c(mm,i,storm))
        sim1 <- rasterFromXYZ(cbind(coords, NAM_pred[,3]+simvals[,i])) * mask_SC
        values(sim1)[which(values(sim1) < 0 )] <- 0
        
        sum_sq_rains[i,mm] <- sum(values(sim1)^2,na.rm = T)
      }
      
      # values(NAM_r)[which(values(NAM_r)<0)] <- 0
      nam_df <- rasterToPoints(NAM_r * mask_SC)[,3]
      obs_df <- rasterToPoints(ST4_r * mask_SC)[,3]
      
      # Obtain log score from each of the 4 sampling schemes for theta
      if(mm==9) all_scores <- log_score_straw(theta_pred, obs_df, nam_df)
      if(mm!=9) all_scores <- log_score(theta_pred, my_dist_mtx, obs_df, nam_df)
      log_scores_marf[basin,mm] <- meanOrig(all_scores)
      log_scores_sum_marf[basin,mm] <- sum(all_scores)
      
      NAM_l <- sum(values(rasterFromXYZ(NAM_pred) * mask_SC),na.rm = T)
      ST4_l <- sum(values(rasterFromXYZ(ST4_pred) * mask_SC),na.rm = T)
      NAM_l_sq <- sum(values(rasterFromXYZ(NAM_pred) * mask_SC)^2,na.rm = T)
      ST4_l_sq <- sum(values(rasterFromXYZ(ST4_pred) * mask_SC)^2,na.rm = T)
      
    }
    
    
    # get bandwidths and y-limits
    bws <- ylims <- probs_i <- probs_s <- probs_t <- c()
    for (mm in 1:nM) {
      
      bws[mm] <- (4*sd(sum_sq_rains[,mm])^5/(3*Ngen))^0.2
      ylims[mm] <-  max(density(sum_sq_rains[,mm], bw = bws[mm])$y)
      
      xx <- density(sum_sq_rains[,mm], bw = bws[mm])$x
      yy <- density(sum_sq_rains[,mm], bw = bws[mm])$y
      
      if(sum(xx <= 0) >0){
        neg_precip <- which(xx <= 0)
        xx <- xx[-neg_precip]
        yy <- yy[-neg_precip]
      }
      
      f <- approxfun(xx, yy)
      C <- integrate(f, min(xx), max(xx))$value #max(xx)*.95
      p.unscaled <- ifelse(ST4_l_sq > max(xx) | ST4_l_sq < min(xx),
                           0,
                           integrate(f, min(xx), ST4_l_sq)$value)
      probs_i[mm] <- round(p.unscaled / C, 4)
      probs_s[mm] <- mean(sum_sq_rains[,mm] <= ST4_l_sq)
      probs_t[mm] <- ifelse(ST4_l_sq > max(xx) | ST4_l_sq < min(xx),
                            0,
                            round(trapzfun(f, min(xx), ST4_l_sq)$value/trapzfun(f, min(xx), max(xx))$value,4))
        
    }
    
    # Make histograms and kernel density estimation plots
    par(mfrow=c(2,3), mar = c(2,2,6,1)+0.1)
    
    # plot the histograms
    for (mm in 1:nM) {
      hist(sum_sq_rains[,mm], main=paste0("total mm precip,\n model ",mm," storm ", storm, ": ", name,
                                          "\n probs (int & sum & trap): \n", probs_i[mm], ", ", probs_s[mm],", ", probs_t[mm]),
           xlim = range(c(sum_sq_rains, ST4_l_sq+100), na.rm = T), prob = T, ylim = c(0, max(ylims)))
      # lines(density(sum_sq_rains[,m], bw = bws[m]))
      abline(v=NAM_l_sq, col="green")
      abline(v=ST4_l_sq, col="blue")
      log_scores[basin,mm] <- logs_sample(ST4_l_sq, sum_sq_rains[,mm])
      CRPSs[basin,mm] <- crps_sample(ST4_l_sq, sum_sq_rains[,mm])
    }
    
    # plot the kernel density estimation
    par(mfrow=c(2,3))
    for (mm in 1:nM) {
      my_dens <- density(sum_sq_rains[,mm], bw = bws[mm], n = 2048)
      pred_val_dens <- my_dens$y[ which(abs(my_dens$x - ST4_l_sq) == min(abs(my_dens$x - ST4_l_sq)))]
      plot(my_dens$x, my_dens$y, type="l", main=paste0("total mm precip,\n model ",mm," storm ", storm, ": ", name,
                                                       "\n probs (int & sum & trap): \n", probs_i[mm], ", ", probs_s[mm],", ", probs_t[mm],
                                                       "\n pred dens val =", round(pred_val_dens,6)),
           xlim = range(c(sum_sq_rains, ST4_l_sq+100), na.rm = T), ylim = c(0, max(ylims)))
      abline(v=NAM_l_sq, col="green")
      abline(v=ST4_l_sq, col="blue")
      points(x = c(ST4_l_sq), y= pred_val_dens, type = "p")
    }
    
  }
  
}

write.csv(cbind(log_scores,bas_locs,coastals), file = paste0("~/NAM-Model-Validation/csv/scores/basinwide/by_region/",state,
                                    "_basins_storm",ste,".csv"))
write.csv(cbind(CRPSs,bas_locs,coastals), file = paste0("~/NAM-Model-Validation/csv/scores/basinwide/by_region/",state,
                                    "_basins_storm",ste,"_CRPSs.csv"))
write.csv(cbind(log_scores_marf,bas_locs,coastals), file = paste0("~/NAM-Model-Validation/csv/scores/basinwide/by_region/alt/",state,
                                         "_basins_storm",ste,"_marf.csv"))
write.csv(cbind(log_scores_sum_marf,bas_locs,coastals), file = paste0("~/NAM-Model-Validation/csv/scores/basinwide/by_region/altsum/",state,
                                         "_basins_storm",ste,"_sum_marf.csv"))

if(make.pdf) dev.off()
