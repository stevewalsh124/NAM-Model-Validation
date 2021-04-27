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

writeRaster(r, "FLmask")

# # shape2 <- read_sf(dsn = ".", layer = "SHAPEFILE")
# shape2 <- st_read("Florida_Drainage_Basins_1997/Florida_Drainage_Basins_1997.shp")
# plot(shape2)
#
# gch <- gConvexHull(methods::as( object = apala, Class = "Spatial" ))
# gb <- gBuffer(methods::as( object = apala, Class = "Spatial" )) 
#
# rss <- readShapeSpatial("Florida_Drainage_Basins_1997/Florida_Drainage_Basins_1997.shp")
