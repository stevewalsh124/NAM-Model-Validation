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
shape <- readOGR(dsn = "Florida_Drainage_Basins_1997", layer = "Florida_Drainage_Basins_1997")
# plot(shape)
apala = subset(shape, HUC %in% c("03120001","03120002","03120003","03120004"))

ext <- extent(apala)
r <- raster(ext, res=.01)  
r <- rasterize(apala, r, field=1)
plot(1:5,ylim=c(27,32), xlim=c(-87,-82),type="n", asp=1)
plot(r,add=T)
US(add=T)

# writeRaster(r, "FLmask")

# shape2 <- read_sf(dsn = ".", layer = "SHAPEFILE")
# shape2 <- st_read("Florida_Drainage_Basins_1997/Florida_Drainage_Basins_1997.shp")
# gch <- gConvexHull(methods::as( object = apala, Class = "Spatial" ))
# gb <- gBuffer(methods::as( object = apala, Class = "Spatial" )) 
# rss <- readShapeSpatial("Florida_Drainage_Basins_1997/Florida_Drainage_Basins_1997.shp")
pdf("pdf/FL_basin.pdf")
FL_mask <- raster("FLmask.grd")
plot(NA, xlim = c(-87,-82), ylim = c(27,32), type="n", asp=1)
plot(FL_mask, add=T)
US(add=T)
mask_FL <- resample(FL_mask, rasterFromXYZ(cbind(coords, simvals[,1])))

for (k in 1:6) {
  print(k)
  load(paste0("RData/prediction",k))
  
  sim1 <- rasterFromXYZ(cbind(coords, NAM_pred[,3]+simvals[,1])) * mask_FL
  plot(sim1, main=paste0("storm ", k, ": ", name)); US(add=T)
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
  
  hist(sum_rains, main=paste0("total sqrt(mm) precip, storm ", k, ": ", name))
  abline(v=NAM_l, col="green")
  abline(v=ST4_l, col="blue")
  
  hist(sum_sq_rains, main=paste0("total mm precip, storm ", k, ": ", name))
  abline(v=NAM_l_sq, col="green")
  abline(v=ST4_l_sq, col="blue")
}
dev.off()
