# Florida river basins from the FL website: 
# https://geodata.dep.state.fl.us/datasets/6ceca90e86c84aeda5ee2b9ac9298518_0

# PERIMETER	    HUC	      EXTHUC	  BASIN	              FEATURE
# 2002276.84709	03130011	99000000	APALACHICOLA RIVER	STREAM

# Read SHAPEFILE.shp from the current working directory (".")
# dsn="directory where the shapefile, projection file, etc are located" 
# layer="name of the file without .shp extention"

require(rgdal)
shape <- readOGR(dsn = "FL_basin", layer = "Florida_Drainage_Basins_1997")

plot(shape)

require(sf)
# shape2 <- read_sf(dsn = ".", layer = "SHAPEFILE")
shape2 <- st_read("FL_basin/Florida_Drainage_Basins_1997.shp")
plot(shape2)

apala = subset(shape2, BASIN=="APALACHICOLA RIVER")

fields::US()
plot(apala, max.plot=1, add=T)
