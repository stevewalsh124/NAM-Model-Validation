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

basin <- 1

make_shp2raster <- F

if(make_shp2raster){
  shape <- readOGR(dsn = "Florida_Drainage_Basins_1997", layer = "Florida_Drainage_Basins_1997")
  # plot(shape)
  if(basin == 1) watershed = subset(shape, HUC %in% c("03120001","03120002","03120003","03120004"))
  if(basin == 2) watershed = subset(shape, HUC %in% c("03130004","03130011","03130012","03130013","03130014"))
  if(basin == 3) watershed = subset(shape, HUC %in% c("03140101","03140102","03140103","03140104",
                                                      "03140105","03140106","03140107"))
  if(basin == 4) watershed = subset(shape, HUC %in% c("03140101","03140102","03140103","03140104","03140105",
                                                      "03140106","03140107","03140202","03140203","03140305"))
  ext <- extent(watershed)
  r <- raster(ext, res=.01)  
  r <- rasterize(watershed, r, field=1)
  plot(1:5,ylim=c(27,32), xlim=c(-87,-82),type="n", asp=1)
  plot(r,add=T)
  US(add=T)
  writeRaster(r, paste0("basin/FLmask",basin), overwrite=T)
}


pdf(paste0("pdf/FL_basin",basin,".pdf"))
FL_mask <- raster(paste0("basin/FLmask",basin,".grd"))
plot(NA, xlim = extent(FL_mask)[1:2] + c(-1,1), ylim = extent(FL_mask)[3:4] + c(-1,1), type="n", asp=1)
plot(FL_mask, add=T)
US(add=T)
load(paste0("RData/prediction4"))
mask_FL <- resample(FL_mask, rasterFromXYZ(cbind(coords, simvals[,1])), method="ngb")

for (k in 1:6) {
  for (m in 1:3) {
    print(paste0("storm ", k, ", model ", m))
    if(m==1) load(paste0("RData/prediction",k))
    if(m==2) load(paste0("RData/prediction",k,"_LM2"))
    if(m==3) load(paste0("RData/prediction",k,"_LM3"))
    
    par(mfrow=c(1,1), mar = c(5,4,4,2)+0.1)
    sim1 <- rasterFromXYZ(cbind(coords, NAM_pred[,3]+simvals[,1])) * mask_FL
    if(m==1) plot(sim1, main=paste0("storm ", k, ": ", name)); US(add=T)
    sum(unique(values(sim1)),na.rm = T)
    
    sum_rains <- sum_sq_rains <- c()
    for (i in 1:Ngen) {
      sim1 <- rasterFromXYZ(cbind(coords, NAM_pred[,3]+simvals[,i])) * mask_FL
      values(sim1)[which(values(sim1) < 0 )] <- 0
      sum_rains[i] <- sum(values(sim1),na.rm = T)
      sum_sq_rains[i] <- sum(values(sim1)^2,na.rm = T)
    }
    
    NAM_l <- sum(values(rasterFromXYZ(NAM_pred) * mask_FL),na.rm = T)
    ST4_l <- sum(values(rasterFromXYZ(ST4_pred[,c(3,4,2)]) * mask_FL),na.rm = T)
    NAM_l_sq <- sum(values(rasterFromXYZ(NAM_pred) * mask_FL)^2,na.rm = T)
    ST4_l_sq <- sum(values(rasterFromXYZ(ST4_pred[,c(3,4,2)]) * mask_FL)^2,na.rm = T)
    
    par(mfrow=c(2,2), mar = c(2,2,2,1)+0.1)
    hist(sum_rains, main=paste0("total sqrt(mm) precip,\n model ",m," storm ", k, ": ", name))
    abline(v=NAM_l, col="green")
    abline(v=ST4_l, col="blue")
    
    hist(sum_sq_rains, main=paste0("total mm precip,\n model ",m," storm ", k, ": ", name))
    abline(v=NAM_l_sq, col="green")
    abline(v=ST4_l_sq, col="blue")
    
    FL_pixs <- dim(rasterToPoints(sim1))[1]
    if(FL_pixs > 0){
      hist((sum_sq_rains*1e3*9834.98)/(FL_pixs), 
      main=paste0("total m^3 in basin, precip,\n model ",m," storm ", k, ": ", name))
      abline(v=(NAM_l_sq*1e3*9834.98)/(FL_pixs), col="green")
      abline(v=(ST4_l_sq*1e3*9834.98)/(FL_pixs), col="blue")
    }
  }
}

dev.off()

mean_sq_rain <- 11472.2/101 # mean mm value for Michael 2018 (s=4) (sum of mm / # pixels)
mean_km_rain <- mean_sq_rain/1e6 # .00011 km, about 4.33 inches in each pixel
total_vol_km3 <- mean_km_rain * 9834.98 # 1.117117 cubic km
total_vol_m3 <- total_vol_km3 * 1e9 # 1.12e9 cubic meters over the watershed
