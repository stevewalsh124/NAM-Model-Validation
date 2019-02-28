library(raster)
library(dplyr)
library(sf) #st_geometry
#harvey landfall: 0000 on Aug 26 2017 [21+4, c(7,8)] #+4 on row to move down 12 hours after landfall (looking at 24hr precip)
#irene  landfall: 1200 on Aug 27 2011 [53+4, c(7,8)] #7lat, 8lon

ST4_irene_12hr_1 <- raster("/Users/stephenwalsh/Desktop/HarveyIreneNAMST4/ST4.2011082718.06h")+
  raster("/Users/stephenwalsh/Desktop/HarveyIreneNAMST4/ST4.2011082800.06h")
ST4_irene_12hr_2 <- raster("/Users/stephenwalsh/Desktop/HarveyIreneNAMST4/ST4.2011082806.06h")+
  raster("/Users/stephenwalsh/Desktop/HarveyIreneNAMST4/ST4.2011082812.06h")

# str(ST4_irene_24hr)
# which(subset(rasterToPoints(ST4_irene_24hr) == max(rasterToPoints(ST4_irene_24hr)))==TRUE)
# dim(rasterToPoints(ST4_irene_24hr))
ST4_irene_12hr_1 <- ST4_irene_12hr_1 %>% projectRaster(crs = "+proj=longlat +datum=WGS84")
ST4_irene_12hr_2 <- ST4_irene_12hr_2 %>% projectRaster(crs = "+proj=longlat +datum=WGS84")
ST4_irene_24hr <- ST4_irene_12hr_1+ST4_irene_12hr_2


NAM_df1 <- rasterToPoints(ST4_irene_12hr_1)
NAM_df2 <- rasterToPoints(ST4_irene_12hr_2)
NAM_df <- rasterToPoints(ST4_irene_12hr_1+ST4_irene_12hr_2)

hurdat <- read.csv("/Users/stephenwalsh/Downloads/R/Research/HURDAT/2017_08HARVEY_12UTC_0823.csv", header = F)
hurdat2 <- read.csv("/Users/stephenwalsh/Downloads/R/Research/HURDAT/2011_09IRENE_0UTC_0821.csv", header = F)
colnames(hurdat) <- colnames(hurdat2) <- c("BTid", "SNum", "Year", "Month", "Day", "Hour", "Lat", "Lon", "Winds", "Pressure")

# harvey_eye <- hurdat[25,7:8]
irene_eye <- hurdat2[55,7:8]
irene_eye2 <- hurdat2[59,7:8]
radius    <- 600
centerLat <- as.numeric(irene_eye[1])
centerLon <- as.numeric(irene_eye[2])
centerLat2 <- as.numeric(irene_eye2[1])
centerLon2 <- as.numeric(irene_eye2[2])

# scale longtitude grid (since it's a fcn of latitude)
lonscale = cos(abs(centerLat*pi/180));
Ydist = (NAM_df[,2]-centerLat)*110.54; # convert to km [i,j]
Xdist = (NAM_df[,1]-centerLon)*111.32*lonscale; #[i,j]
radialdist = abs(sqrt(Ydist^2+Xdist^2));
buffer_precip <- ifelse(radialdist < radius, NAM_df[,3], NA)

NAM_df_buffer <- cbind(NAM_df[,1:2],buffer_precip)
raster_buffer <- rasterFromXYZ(NAM_df_buffer)

plot(ST4_irene_24hr, xlim=c(centerLon-6, centerLon+6), ylim=c(centerLat-6,centerLat+6));plot(st_geometry(spData::us_states), add=TRUE)
plot(raster_buffer , xlim=c(centerLon-6, centerLon+6), ylim=c(centerLat-6,centerLat+6));plot(st_geometry(spData::us_states), add=TRUE)


################try to subset the file name to obtain the date/time and use this to access hurdat row
ST4file
#The first 4 ST4 are 24 hrs. Two buffers for each 12hr group, using the 6th and 18th hour as centers (1st and 3rd files in ST4 folder)
eye_1 <- list.files(ST4_folder,full.names = T)[1]
#eye_2 <- list.files(ST4_folder,full.names = T)[3] #just add 12 to eye_1 hour later, aka look 2 rows down

time_start <- tail(str_locate_all(pattern ='/', eye_1)[[1]],1)[,2]
storm_time <-  substr(eye_1,time_start[1]+5,nchar(eye_1))
storm_year <- substr(storm_time,1,4)
storm_month <- substr(storm_time,5,6)
storm_day <- substr(storm_time,7,8)
storm_hour <- substr(storm_time,9,10)

print(paste("Storm is", storm_name, "Year is",storm_year,"Month is",storm_month,"Day is",storm_day,"Hour is",storm_hour))
eye1_latlon <- hurdat[hurdat$Day==as.integer(storm_day)&hurdat$Hour==as.integer(storm_hour),7:8]
eye2_latlon <- hurdat[hurdat$Day==as.integer(storm_day)&hurdat$Hour==as.integer(storm_hour)+6,7:8]

storm_name <- substr(ST4_folder,name_start[1]+5,nchar(ST4_folder))
class(hurdat$Day)

date_start <- str_locate_all(pattern ='_218_', NAMfile)[[1]][1,2] + 1
storm_date <- as.numeric(substr(NAMfile,date_start,date_start+7))


#################ST4 buffer function
ST4bufferprecip <- function(NAM_df,LatLonVec,LatLonVec2,radius){
  centerLat <- as.numeric(LatLonVec[1])
  centerLon <- as.numeric(LatLonVec[2])
  centerLat2 <- as.numeric(LatLonVec2[1])
  centerLon2 <- as.numeric(LatLonVec2[2])
  
  # scale longtitude grid (since it's a fcn of latitude)
  lonscale = cos(abs(centerLat*pi/180));
  Ydist = (NAM_df[,2]-centerLat)*110.54; # convert to km [i,j]
  Xdist = (NAM_df[,1]-centerLon)*111.32*lonscale; #[i,j]
  radialdist = abs(sqrt(Ydist^2+Xdist^2));
  
  lonscale2 = cos(abs(centerLat2*pi/180));
  Ydist2 = (NAM_df[,2]-centerLat2)*110.54; # convert to km [i,j]
  Xdist2 = (NAM_df[,1]-centerLon2)*111.32*lonscale; #[i,j]
  radialdist2 = abs(sqrt(Ydist2^2+Xdist2^2));
  
  buffer_precip <- ifelse(radialdist < radius, NAM_df[,3], 0)
  buffer_precip2 <- ifelse(radialdist < radius|radialdist2<radius, buffer_precip, NA)
  
  ST4_df_buffer <- cbind(NAM_df[,1:2],buffer_precip2)
  raster_buffer <- rasterFromXYZ(ST4_df_buffer)
  return(raster_buffer)
}
dim(trythis)
dim(ST4_irene_12hr_1)
trythis <- ST4bufferprecip(NAM_df1, irene_eye,irene_eye2,500)
trythis2 <- ST4bufferprecip(NAM_df2, irene_eye2,irene_eye,500)

plot(trythis,xlim=c(min(centerLon,centerLon2)-5.5, max(centerLon,centerLon2)+5.5), 
     ylim=c(min(centerLat,centerLat2)-5.5,max(centerLat,centerLat2)+5.5)); plot(st_geometry(spData::us_states), add=TRUE)
plot(trythis2,xlim=c(min(centerLon,centerLon2)-5.5, max(centerLon,centerLon2)+5.5), 
     ylim=c(min(centerLat,centerLat2)-5.5,max(centerLat,centerLat2)+5.5)); plot(st_geometry(spData::us_states), add=TRUE)
plot(trythis+trythis2,xlim=c(min(centerLon,centerLon2)-5.5, max(centerLon,centerLon2)+5.5), 
     ylim=c(min(centerLat,centerLat2)-5.5,max(centerLat,centerLat2)+5.5)); plot(st_geometry(spData::us_states), add=TRUE)
