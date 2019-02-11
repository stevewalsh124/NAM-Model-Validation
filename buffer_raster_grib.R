library(raster)
library(dplyr)
library(sf) #st_geometry
#harvey landfall: 0000 on Aug 26 2017 [21+4, c(7,8)] #+4 on row to move down 12 hours after landfall (looking at 24hr precip)
#irene  landfall: 1200 on Aug 27 2011 [53+4, c(7,8)] #7lat, 8lon

ST4_irene_24hr <- raster("/Users/stephenwalsh/Desktop/HarveyIreneNAMST4/ST4.2011082718.06h")+
  raster("/Users/stephenwalsh/Desktop/HarveyIreneNAMST4/ST4.2011082800.06h")+
  raster("/Users/stephenwalsh/Desktop/HarveyIreneNAMST4/ST4.2011082806.06h")+
  raster("/Users/stephenwalsh/Desktop/HarveyIreneNAMST4/ST4.2011082812.06h")

str(ST4_irene_24hr)
which(subset(rasterToPoints(ST4_irene_24hr)== max(rasterToPoints(ST4_irene_24hr)))==TRUE)
dim(rasterToPoints(ST4_irene_24hr))

NAM_df <- rasterToPoints(ST4_irene_24hr)
prod(dim(ST4_irene_24hr))
ST4_irene_24hr <- ST4_irene_24hr %>% projectRaster(crs = "+proj=longlat +datum=WGS84")
ST4_irene <- (raster("/Users/stephenwalsh/Desktop/HarveyIreneNAMST4/ST4.2011082718.06h") %>% projectRaster(crs = "+proj=longlat +datum=WGS84"))
plot(ST4_irene_24hr);plot(st_geometry(spData::us_states), add=TRUE)

hurdat <- read.csv("/Users/stephenwalsh/Downloads/R/Research/HURDAT/2017_08HARVEY_12UTC_0823.csv", header = F)
hurdat2 <- read.csv("/Users/stephenwalsh/Downloads/R/Research/HURDAT/2011_09IRENE_0UTC_0821.csv", header = F)
colnames(hurdat) <- colnames(hurdat2) <- c("BTid", "SNum", "Year", "Month", "Day", "Hour", "Lat", "Lon", "Winds", "Pressure")

# harvey_eye <- hurdat[25,7:8]
irene_eye <- hurdat2[57,7:8]

radius    <- 500
centerLat <- as.numeric(irene_eye[1])
centerLon <- as.numeric(irene_eye[2])

# scale longtitude grid (since it's a fcn of latitude)
lonscale = cos(abs(centerLat*pi/180));
Ydist = (NAM_df[,2]-centerLat)*110.54; # convert to km [i,j]
Xdist = (NAM_df[,1]-centerLon)*111.32*lonscale; #[i,j]
radialdist = abs(sqrt(Ydist^2+Xdist^2));
buffer_precip <- ifelse(radialdist < radius, NAM_df[,3], NA)

NAM_df_buffer <- cbind(NAM_df[,1:2],buffer_precip)
raster_buffer <-rasterFromXYZ(NAM_df_buffer)
plot(raster_buffer);plot(st_geometry(spData::us_states), add=TRUE)


# Turn the matrix into a raster
rast <- raster(xy)
# Give it lat/lon coords for 36-37°E, 3-2°S
extent(rast) <- c(36,37,-3,-2)
# ... and assign a projection
projection(rast) <- CRS("+proj=longlat +datum=WGS84")
plot(rast)

plot(raster_buffer, xlim=c(centerLon-6, centerLon+6), ylim=c(centerLat-6,centerLat+6));plot(st_geometry(spData::us_states), add=TRUE)
plot(ST4_irene_24hr, xlim=c(centerLon-6, centerLon+6), ylim=c(centerLat-6,centerLat+6));plot(st_geometry(spData::us_states), add=TRUE)
