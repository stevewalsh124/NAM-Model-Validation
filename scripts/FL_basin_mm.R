# Florida river basins from the FL website: 
# https://geodata.dep.state.fl.us/datasets/6ceca90e86c84aeda5ee2b9ac9298518_0

# PERIMETER	    HUC	      EXTHUC	  BASIN	              FEATURE
# 2002276.84709	03130011	99000000	APALACHICOLA RIVER	STREAM

# Read SHAPEFILE.shp from the current working directory (".")
# dsn="directory where the shapefile, projection file, etc are located" 
# layer="name of the file without .shp extention"

require(rgdal)
require(sf)
require(raster)
library(fields)

make.pdf <- F

basin <- 1

# number of models compared
nM <- 3
# Ngen <- 1000 is loaded in the prediction RData

if(make.pdf) pdf(paste0("pdf/FL_basin",basin,"_mm_only.pdf"))
FL_mask <- raster(paste0("basin/FLmask",basin,".grd"))
plot(NA, xlim = extent(FL_mask)[1:2] + c(-1,1), ylim = extent(FL_mask)[3:4] + c(-1,1), type="n", asp=1)
plot(FL_mask, add=T)
US(add=T)
load(paste0("RData/prediction4"))
mask_FL <- raster::resample(FL_mask, rasterFromXYZ(cbind(coords, simvals[,1])), method="ngb")

sum_sq_rains <- matrix(NA, Ngen, nM)
sum_sq_straw <- c()

for (k in 4) {
  for (m in 1:nM) {
    print(paste0("storm ", k, ", model ", m))
    if(m==1) load(paste0("RData/prediction",k))
    if(m==2) load(paste0("RData/prediction",k,"_LM2"))
    if(m==3) load(paste0("RData/prediction",k,"_LM3"))
    
    sim1 <- rasterFromXYZ(cbind(coords, NAM_pred[,3]+simvals[,1])) * mask_FL
    
    if(m==1){
      par(mfrow=c(1,1), mar = c(5,4,4,2)+0.1)
      plot(sim1, main=paste0("storm ", k, ": ", name)); US(add=T)
    }
    sum(unique(values(sim1)),na.rm = T)
    
    for (i in 1:Ngen) {
      # spatially correlated errors
      sim1 <- rasterFromXYZ(cbind(coords, NAM_pred[,3]+simvals[,i])) * mask_FL
      values(sim1)[which(values(sim1) < 0 )] <- 0
      sum_sq_rains[i,m] <- sum(values(sim1)^2,na.rm = T)
    }
    
    if(m==1){
      #strawman: white noise errors; no spatial correlation
      for(i in 1:Ngen){
        sim1 <- rasterFromXYZ(cbind(coords, NAM_pred[,3] + 
                                      rnorm(n = nrow(NAM_pred), mean = 0, 
                                            sd = sqrt(mean(exp(theta_pred[,2])))))) * mask_FL
        values(sim1)[which(values(sim1) < 0 )] <- 0
        sum_sq_straw[i] <- sum(values(sim1)^2,na.rm = T)
      }  
    }
    
    NAM_l <- sum(values(rasterFromXYZ(NAM_pred) * mask_FL),na.rm = T)
    ST4_l <- sum(values(rasterFromXYZ(ST4_pred[,c(3,4,2)]) * mask_FL),na.rm = T)
    NAM_l_sq <- sum(values(rasterFromXYZ(NAM_pred) * mask_FL)^2,na.rm = T)
    ST4_l_sq <- sum(values(rasterFromXYZ(ST4_pred[,c(3,4,2)]) * mask_FL)^2,na.rm = T)
    
  }
  
  # Make histograms and kernel density estimation plots
  par(mfrow=c(1,3), mar = c(2,2,6,1)+0.1)
  
  # get bandwidths and y-limits
  bws <- ylims <- probs_i <- probs_s <- c()
  for (m in 1:nM) {
    bws[m] <- (4*sd(sum_sq_rains[,m])^5/(3*Ngen))^0.2
    ylims[m] <-  max(density(sum_sq_rains[,m], bw = bws[m])$y)
    
    xx <- density(sum_sq_rains[,m], bw = bws[m])$x
    yy <- density(sum_sq_rains[,m], bw = bws[m])$y
    f <- approxfun(xx, yy)
    C <- integrate(f, min(xx), max(xx))$value
    p.unscaled <- integrate(f, min(xx), ST4_l_sq)$value
    probs_i[m] <- round(p.unscaled / C, 4)
    
    probs_s[m] <- mean(sum_sq_rains[,m] <= ST4_l_sq)
  }
  
  # (nM + 1)st is for the straw man; no spatial correlation in errors 
  bws[m+1] <- (4*sd(sum_sq_straw)^5/(3*Ngen))^0.2
  ylims[m+1] <-  max(density(sum_sq_straw, bw = bws[m+1])$y)
  
  xx <- density(sum_sq_straw, bw = bws[m+1])$x
  yy <- density(sum_sq_straw, bw = bws[m+1])$y
  f <- approxfun(xx, yy)
  C <- integrate(f, min(xx), max(xx))$value
  if(min(xx) - bws[m+1] < ST4_l_sq){
    p.unscaled <- integrate(f, min(xx), ST4_l_sq)$value
  } else { p.unscaled <- 0 }
  probs_i[m] <- round(p.unscaled / C, 4)
  
  probs_s[m] <- mean(sum_sq_rains[,m] <= ST4_l_sq)
  
  # plot the histograms
  for (m in 1:nM) {
    hist(sum_sq_rains[,m], main=paste0("total mm precip,\n model ",m," storm ", k, ": ", name,
                                       "\n probs (int & sum): ", probs_i[m], ", ", probs_s[m]),
         xlim = range(sum_sq_rains, na.rm = T), prob = T, ylim = c(0, max(ylims)))
    # lines(density(sum_sq_rains[,m], bw = bws[m]))
    abline(v=NAM_l_sq, col="green")
    abline(v=ST4_l_sq, col="blue")
  }
  
  # plot the kernel density estimation
  for (m in 1:nM) {
    my_dens <- density(sum_sq_rains[,m], bw = bws[m], n = 2048)
    pred_val_dens <- my_dens$y[ which(abs(my_dens$x - ST4_l_sq) == min(abs(my_dens$x - ST4_l_sq)))]
    plot(my_dens$x, my_dens$y, type="l", main=paste0("total mm precip,\n model ",m," storm ", k, ": ", name,
                                                     "\n probs (int & sum): ", probs_i[m], ", ", probs_s[m],
                                                     "\n pred dens val =", round(pred_val_dens,6)),
         xlim = range(sum_sq_rains, na.rm = T), ylim = c(0, max(ylims)))
    abline(v=NAM_l_sq, col="green")
    abline(v=ST4_l_sq, col="blue")
    points(x = c(ST4_l_sq), y= pred_val_dens, type = "p")
  }
  m <- m + 1
  my_dens <- density(sum_sq_straw, bw = bws[m], n = 2048)
  pred_val_dens <- my_dens$y[ which(abs(my_dens$x - ST4_l_sq) == min(abs(my_dens$x - ST4_l_sq)))]
  plot(my_dens$x, my_dens$y, type="l", main=paste0("total mm precip,\n model ",m," storm ", k, ": ", name,
                                                   "\n probs (int & sum): ", probs_i[m], ", ", probs_s[m],
                                                   "\n pred dens val =", round(pred_val_dens,6)),
       xlim = range(sum_sq_rains, na.rm = T), ylim = c(0, max(ylims)))
  abline(v=NAM_l_sq, col="green")
  abline(v=ST4_l_sq, col="blue")
  points(x = c(ST4_l_sq), y= pred_val_dens, type = "p")
  
}

if(make.pdf) dev.off()
