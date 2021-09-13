# Plotting 47 Storms' Square Root Precipitation for First 24 Hours After Landfall
# Steve Walsh

library(sp)
library(ncdf4)
library(sf)
library(dplyr)
library(stringr) #str_locate_all "_218_" in the NAM file names
library(rgdal)
library(colorspace)
library(ggplot2)
library(SpatialEpi) #latlong2grid()
library(geoR) #as.geodata, variog, variofit
library(fields) #US()
library(raster)

# Imagine that the two circular buffer areas create a Venn diagram shape.
# Do you want it so that A\B and B\A have 12hr of precip while A intersect B has 24hr? (choose F)
# Do you want the entire union of A and B to have 24 hrs of precip? (choose T)
addbothbuffers = T
pred <- T

# Do you want to write csvs and rasters, as well as make a pdf?
write.pdf <- T

# Do you want to subtract the pointwise error field mean 
# (made from all 47 storms) from each error field? (If yes choose T)
subtractPWmean = F
makePWmean = F
if(makePWmean & subtractPWmean) stop("makePWmean and subtractPWmean can't both be TRUE")

if(write.pdf){
pdf(paste0("~/NAM-Model-Validation/pdf/sqrtbuffer_47_ngb_",
           if(makePWmean){"makePWmean"},if(subtractPWmean){"subtractPWmean"},if(pred){"pred"},
           "_bothbuffers.pdf"))
}

if(pred){
  storm.dirs <- list.dirs("~/NAM-Model-Validation/prediction", recursive = F)
  if(length(storm.dirs)==0){
    storm.dirs <- list.dirs("/Volumes/LACIEHD/NAMandST4_prediction", recursive = F) }
} else {storm.dirs <- list.dirs("~/NAMandST4", recursive = F) }

storms.out.of.hurdat <- c()

# Change this file so that it corresponds to the new buffering
PW_mean <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_mean.grd")

#[-c(2,3,6,8)]#THE NON12HR STORMS
#[c(5,16,22,25,35,47)] #R session aborted when fitting Matern
#[-c(2,3,6,8,5,16,22,25,35,47)]# to avoid aborting R and storms w/o 12hr precip
#[-c(1,12:36)] #to avoid JPEG2000 driver issue on Mac

# these will be lists of lists (#lists=#storms),
# where each element in a list is rasters for a storm
NAMlist <- list()
ST4list <- list()

#collecting the parameter estimates
SVGparamest <- SVGparamboot <- list()

phivec    <- c()
sig2vec   <- c()
tau2vec   <- c()
prRangevec<- c()
namevec <- c()

# load, modify and project land/sea mask 
mask <- raster("~/NAM-Model-Validation/lsmask.nc")
mask[mask==-1]  <- NA #changed from 0 to NA because mismatch rows due to off-coast pts
extent(mask)[1] <- extent(mask)[1]-360
extent(mask)[2] <- extent(mask)[2]-360
mask.regrid <- raster::resample(mask, projectRaster(raster(
  "~/NAM-Model-Validation/nam_218_20170826_0000_012.grb2"),
  crs = "+proj=longlat +datum=WGS84"), method='ngb')  #/Volumes/LACIEHD/

# Multiple plot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Buffer function (two eye locations)
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
  
  if(addbothbuffers){buffer_precip2 <- ifelse(radialdist < radius|radialdist2<radius, NAM_df[,3], NA)}
  else{buffer_precip <- ifelse(radialdist < radius, NAM_df[,3], 0)
  buffer_precip2 <- ifelse(radialdist < radius|radialdist2<radius, buffer_precip, NA)}
  
  ST4_df_buffer <- cbind(NAM_df[,1:2],buffer_precip2)
  raster_buffer <- rasterFromXYZ(ST4_df_buffer)
  return(raster_buffer)
}


ST4bufferprecip_addboth <- function(NAM_df,LatLonVec,LatLonVec2,radius){
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
  
  # buffer_precip <- ifelse(radialdist < radius, NAM_df[,3], 0)
  buffer_precip2 <- ifelse(radialdist < radius|radialdist2<radius, NAM_df[,3], NA)
  
  ST4_df_buffer <- cbind(NAM_df[,1:2],buffer_precip2)
  raster_buffer <- rasterFromXYZ(ST4_df_buffer)
  return(raster_buffer)
}

#Combine multiple rasters in all locations, not just their intersection
mosaicList <- function(rasList){
  #Internal function to make a list of raster objects from list of files.
  ListRasters <- function(list_names) {
    raster_list <- list() # initialise the list of rasters
    for (i in 1:(length(list_names))){
      grd_name <- list_names[i] # list_names contains all the names of the images in .grd format
      raster_file <- raster::raster(grd_name)
    }
    raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
  }
  
  #convert every raster path to a raster object and create list of the results
  raster.list <-sapply(rasList, FUN = ListRasters)
  
  # edit settings of the raster list for use in do.call and mosaic
  names(raster.list) <- NULL
  #####This function deals with overlapping areas
  raster.list$fun <- sum
  
  #run do call to implement mosaic over the list of raster objects.
  mos <- do.call(raster::mosaic, raster.list)
  
  #set crs of output
  crs(mos) <- crs(x = raster(rasList[1]))
  return(mos)
}


for(i in 1:length(storm.dirs)){
  # i <- 1
  currentST4_NAM_folders <- list.dirs(storm.dirs[i])[-1] #[-1] avoids parent folder
  ST4_folder <- currentST4_NAM_folders[1]
  NAM_folder <- currentST4_NAM_folders[2] #uncomment the [2] as well
  storm_yearname <- list.dirs(storm.dirs[i], full.names = F)[-1][1]
  
  name_start <- tail(str_locate_all(pattern ='/', ST4_folder)[[1]],1)
  storm_name <- substr(ST4_folder,name_start[1]+5,nchar(ST4_folder))
  namevec[i] <- storm_name
  print(paste0("Storm #",i,": ",storm_name))
  
  hurdat_csv <- list.files(storm.dirs[i], pattern = ".csv")
  
  print(length(list.files(NAM_folder))) #should be 15 or 18, depending if 72hr timestep exists
  NAMlist[[i]] <- list()
  ST4list[[i]] <- list()
  
  #Identifying where in the NAM folder the first 24 hours after landfall files are (latest timestep)
  if(length(list.files(NAM_folder,full.names = T))==18){first24 <- c(7,8)} else {
    if(length(list.files(NAM_folder,full.names = T))==15){first24 <- c(5,6)} else {
      if(length(list.files(NAM_folder,full.names = T))==12&storm_name=="rita"){first24 <- c(5,6)} else {
        if(length(list.files(NAM_folder, full.names = T))==53){first24 <- c(13,25)} else {
          print("Total NAM files for this storm not 15 or 18 or 53; skipping"); next}}}} 
  
  for (k in 1:length(list.files(ST4_folder,full.names = T)[1:4])) {
    # k <-1
    ST4file <- list.files(ST4_folder,full.names = T)[k]
    ST4list[[i]][[k]] <- (raster(ST4file) %>% projectRaster(crs = "+proj=longlat +datum=WGS84", method = "ngb"))
    ST4list[[i]][[k]] <- raster::resample(ST4list[[i]][[k]], mask.regrid, method = "ngb")
    
    print(min(values(ST4list[[i]][[k]]), na.rm=T))
    # values(ST4list[[i]][[k]])[which(values(ST4list[[i]][[k]]<=1))] <- 1 
    #change nonpositive projected precipitation to 0.0000000001 to implement log precip wihtout -Inf
  }   
  
  
  for (j in 1:length(list.files(NAM_folder,full.names = T)[first24])) {
    # j <-1
    NAMfile <- list.files(NAM_folder,full.names = T)[first24[j]]
    
    #Establish correct NAM band for total precip, mostly depends on date of storm (3 w issues)
    date_start <- str_locate_all(pattern ='_218_', NAMfile)[[1]][1,2] + 1
    storm_date <- as.numeric(substr(NAMfile,date_start,date_start+7))
    if(storm_date < 20060701){NAM_band <- 10} else 
      if(storm_date > 20060701 && storm_date < 20170101){NAM_band <- 11} else {NAM_band <- 373}
    # if(storm_name=="charley"&&j==1&&storm_date<20050000){NAM_band <- 13} 
    # if(storm_name=="charley"&&j==2&&storm_date<20050000){NAM_band <- 14} 
    # if(storm_name=="ivan"   &&j==1&&storm_date<20050000){NAM_band <- 12} 
    # if(storm_name=="ivan"   &&j==2&&storm_date<20050000){NAM_band <- 14} # these were resolved
    # if(storm_name=="matthew"&&j==1&&storm_date<20050000){NAM_band <- 14} # by getting the "hi"
    # if(storm_name=="matthew"&&j==2&&storm_date<20050000){NAM_band <- 12} # files for 2004 NAM
    # if(storm_name=="frances"&&j==1&&storm_date<20050000){NAM_band <- 13} 
    # if(storm_name=="frances"&&j==2&&storm_date<20050000){NAM_band <- 14} 
    
    print(NAM_band)
    NAMlist[[i]][[j]] <- (raster(NAMfile, band=NAM_band) %>% projectRaster(crs = "+proj=longlat +datum=WGS84",
                                                                           method = "ngb"))
    print(min(values(NAMlist[[i]][[j]]), na.rm=T))
    # values(NAMlist[[i]][[j]])[which(values(NAMlist[[i]][[j]]<=1))] <- 1 
    #change nonpositive projected precipitation to 0.0000000001 to implement log precip without -Inf
  }
  
  
  ##############################################################################################################
  ST4_first12  <- ST4list[[i]][[1]]+ST4list[[i]][[2]]
  ST4_second12 <- ST4list[[i]][[3]]+ST4list[[i]][[4]]
  
  ST4_first12_logmask <- ((ST4_first12))*mask.regrid #removed log
  values(ST4_first12_logmask)[which(values(ST4_first12_logmask<=0))] <- 0 #log precip scale is nonnegative
  ST4_second12_logmask <- ((ST4_second12))*mask.regrid #removed log
  values(ST4_second12_logmask)[which(values(ST4_second12_logmask<=0))] <- 0 
  
  ST4_df_first12 <- rasterToPoints(ST4_first12_logmask)
  ST4_df_second12 <- rasterToPoints(ST4_second12_logmask)
  
  NAM_first12  <- NAMlist[[i]][[1]]
  NAM_second12 <- NAMlist[[i]][[2]]
  
  NAM_first12_logmask <- ((NAM_first12))*mask.regrid #removed log
  # values(NAM_first12_logmask)[which(values(NAM_first12_logmask<=0))] <- 0 #log precip scale is nonnegative
  NAM_second12_logmask <- ((NAM_second12))*mask.regrid #removed log
  # values(NAM_second12_logmask)[which(values(NAM_second12_logmask<=0))] <- 0 
  
  NAM_df_first12  <- rasterToPoints(NAM_first12_logmask)
  NAM_df_second12 <- rasterToPoints(NAM_second12_logmask)
  PW_mean_df <- rasterToPoints(PW_mean)
  
  #The first 4 ST4 are 24 hrs. Two buffers for each 12hr group, 
  #using the 6th and 18th hour as centers (1st and 3rd files in ST4 folder)
  eye_1 <- list.files(ST4_folder,full.names = T)[1]
  eye_2 <- list.files(ST4_folder,full.names = T)[3] # add 12 hours to eye_1, aka look 2 rows down
  
  if(pred) {
    best_track <- st_read(list.files(storm.dirs[i], full.names = T, pattern = "pts.shp")[1])
    eye1_latlon <- apply(st_coordinates(best_track)[1:2,],2,mean)[2:1]
    eye2_latlon <- apply(st_coordinates(best_track)[3:2,],2,mean)[2:1]
  } else {
    hurdat_file <- list.files(storm.dirs[i], pattern = ".csv", recursive = F, full.names = T)[1]
    hurdat <- read.csv(hurdat_file, header = F)
    colnames(hurdat) <- c("BTid", "SNum", "Year", "Month", "Day", "Hour", "Lat", "Lon", "Winds", "Pressure")
    
    time_start  <- tail(str_locate_all(pattern ='/', eye_1)[[1]],1)[,2]
    storm_time  <- substr(eye_1,time_start[1]+5,nchar(eye_1))
    storm_year  <- substr(storm_time,1,4)
    storm_month <- substr(storm_time,5,6)
    storm_day   <- substr(storm_time,7,8)
    storm_hour  <- substr(storm_time,9,10)
    
    time_start2 <- tail(str_locate_all(pattern ='/', eye_2)[[1]],1)[,2]
    storm_time2 <- substr(eye_2,time_start2[1]+5,nchar(eye_2))
    storm_year2 <- substr(storm_time2,1,4)
    storm_month2<- substr(storm_time2,5,6)
    storm_day2  <- substr(storm_time2,7,8)
    storm_hour2 <- substr(storm_time2,9,10)
    
    print(paste("Storm is", storm_name, "Year is",storm_year,"Month is",storm_month,
                "Day is",storm_day,"Hour is",storm_hour))
    eye1_latlon <- hurdat[hurdat$Day==as.integer(storm_day) &
                            hurdat$Hour==as.integer(storm_hour), 7:8]
    eye2_latlon <- hurdat[hurdat$Day==as.integer(storm_day2)&
                            hurdat$Hour==as.integer(storm_hour2),7:8]
    
    if(nrow(eye1_latlon)==0){print("eye 1 hurdat trouble; using last")
      eye1_latlon <- tail(hurdat,1)[,7:8]
      storms.out.of.hurdat <- c(storms.out.of.hurdat, i)}
    if(nrow(eye2_latlon)==0){print("eye 2 hurdat trouble; using last")
      eye2_latlon <- tail(hurdat,1)[,7:8]
      storms.out.of.hurdat <- c(storms.out.of.hurdat, i)}
  }
  
  date_start <- str_locate_all(pattern ='_218_', NAMfile)[[1]][1,2] + 1
  storm_date <- as.numeric(substr(NAMfile,date_start,date_start+7))
  
  #######################################################################
  # Original radius 700, changed to 250 and 500 for efficiency in trials runs with simulating error fields
  radius <- 700
  ST4trythis  <- ST4bufferprecip(ST4_df_first12, eye1_latlon, eye2_latlon, radius)
  ST4trythis2 <- ST4bufferprecip(ST4_df_second12, eye2_latlon, eye1_latlon, radius)
  
  NAMtrythis  <- ST4bufferprecip(NAM_df_first12, eye1_latlon, eye2_latlon, radius)
  NAMtrythis2 <- ST4bufferprecip(NAM_df_second12, eye2_latlon, eye1_latlon, radius)
  
  PW_mean_buff1 <- ST4bufferprecip(PW_mean_df, eye1_latlon, eye2_latlon, radius)
  PW_mean_buff2 <- ST4bufferprecip(PW_mean_df, eye2_latlon, eye1_latlon, radius)
  PW_mean_buff <- PW_mean_buff1 + PW_mean_buff2
  
  NAM_plotter <- sqrt(NAMtrythis+NAMtrythis2)
  # values(NAM_plotter)[which(values(NAM_plotter<=0))] <- 0 
  
  ST4_plotter <- sqrt(ST4trythis+ST4trythis2)
  # values(ST4_plotter)[which(values(ST4_plotter<=0))] <- 0   
  
  ## Find the extra points from the files and remove them
  ## When NAM and ST4 have different amounts of pixels
  ## if NAM has more rows...
  if(nrow(rasterToPoints(NAM_plotter)) > nrow(rasterToPoints(ST4_plotter))){
    ST40 <- ST4_plotter
    values(ST40)[!is.na(values(ST40))] <- 0
    # plot(ST40)
    NAM_plotter <- NAM_plotter - ST40
  }
  
  ## if ST4 has more rows
  if(nrow(rasterToPoints(NAM_plotter)) < nrow(rasterToPoints(ST4_plotter))){
    NAM0 <- NAM_plotter
    values(NAM0)[!is.na(values(NAM0))] <- 0
    # plot(NAM0)
    ST4_plotter <- ST4_plotter - NAM0
  }
  
  if(nrow(rasterToPoints(NAM_plotter)) > nrow(rasterToPoints(ST4_plotter))){
    ST40 <- ST4_plotter
    values(ST40)[!is.na(values(ST40))] <- 0
    # plot(ST40)
    NAM_plotter <- NAM_plotter - ST40
  }
  
  # 24 HOUR PRECIP ANALYSIS
  # Convert the three rasters to data frames for ggplot
  NAM_spdf <- as(NAM_plotter, "SpatialPixelsDataFrame")
  NAM_df <- as.data.frame(NAM_spdf)
  colnames(NAM_df) <- c("value", "x", "y")
  
  ST4_spdf <- as(ST4_plotter, "SpatialPixelsDataFrame")
  ST4_df <- as.data.frame(ST4_spdf)
  colnames(ST4_df) <- c("value", "x", "y")
  
  # Use these to load into simulate_error_field.R, prediction, etc
  if(write.pdf){
  if(pred){
    if(!dir.exists(paste0("~/NAM-Model-Validation/csv/prediction_sqrt/",storm_yearname,"/"))) {
      dir.create(paste0("~/NAM-Model-Validation/csv/prediction_sqrt/",storm_yearname,"/"), recursive = T)
    }
    write.csv(NAM_df, paste0("~/NAM-Model-Validation/csv/prediction_sqrt/",storm_yearname,"/",storm_yearname,"_NAMdf.csv"))
    write.csv(ST4_df, paste0("~/NAM-Model-Validation/csv/prediction_sqrt/",storm_yearname,"/",storm_yearname,"_ST4df.csv"))
  } else {
    if(!dir.exists("~/NAM-Model-Validation/csv/nam_df_sqrt/")) {
      dir.create("~/NAM-Model-Validation/csv/nam_df_sqrt/", recursive = T)
    }
    if(!dir.exists("~/NAM-Model-Validation/csv/st4_df_sqrt/")) {
      dir.create("~/NAM-Model-Validation/csv/st4_df_sqrt/", recursive = T)
    }
    write.csv(NAM_df, paste0("~/NAM-Model-Validation/csv/nam_df_sqrt/namdf_",storm_year,storm_name,"_",radius,"km.csv"))
    write.csv(ST4_df, paste0("~/NAM-Model-Validation/csv/st4_df_sqrt/st4df_",storm_year,storm_name,"_",radius,"km.csv"))
  }
  }
  
  if(subtractPWmean){
    error <- ST4_plotter - NAM_plotter - PW_mean_buff
    
  } else {
    error <- ST4_plotter - NAM_plotter
  }
  

  if(makePWmean){
    if(!dir.exists("~/NAM-Model-Validation/error_rasters_sqrt/")) {
      dir.create("~/NAM-Model-Validation/error_rasters_sqrt/", recursive = T)
    }
    if(!dir.exists("~/NAM-Model-Validation/error_rasters_squared_sqrt/")) {
      dir.create("~/NAM-Model-Validation/error_rasters_squared_sqrt/", recursive = T)
    }
    if(!dir.exists("~/NAM-Model-Validation/error_rasters_counts_sqrt/")) {
      dir.create("~/NAM-Model-Validation/error_rasters_counts_sqrt/", recursive = T)
    }
    if(write.pdf){
    #write rasters from the errors to combine all on common map
    writeRaster(error, paste0("~/NAM-Model-Validation/error_rasters_sqrt/",storm_year,storm_name), overwrite=T) 
    writeRaster(error*error, paste0("~/NAM-Model-Validation/error_rasters_squared_sqrt/",storm_year,storm_name), overwrite=T)
    }
    # in the loop for I in the storms
    error_count <- error
    error_count[!is.na(error_count)] <- 1
    error_count[is.na(error_count)] <- 0
    if(write.pdf){
    writeRaster(error_count, paste0("~/NAM-Model-Validation/error_rasters_counts_sqrt/",storm_year,storm_name), overwrite=T)
    }
  }
  
  error_spdf <- as((error), "SpatialPixelsDataFrame")
  error_df <- as.data.frame(error_spdf)
  colnames(error_df) <- c("value", "x", "y")
  
  if(write.pdf){
  if(!pred){
    if(subtractPWmean){
      if(!dir.exists("~/NAM-Model-Validation/csv/error_df_sqrt/subtractPWmeanT_flat/")) {
        dir.create("~/NAM-Model-Validation/csv/error_df_sqrt/subtractPWmeanT_flat/", recursive = T)
      }
      write.csv(error_df, paste0("~/NAM-Model-Validation/csv/error_df_sqrt/subtractPWmeanT_flat/errordf_PW_",
                                 storm_year,storm_name,"_",radius,"deg.csv"))
    } else {
      if(!dir.exists("~/NAM-Model-Validation/csv/error_df_sqrt/subtractPWmeanF/")) {
        dir.create("~/NAM-Model-Validation/csv/error_df_sqrt/subtractPWmeanF/", recursive = T)
      }
      write.csv(error_df, paste0("~/NAM-Model-Validation/csv/error_df_sqrt/subtractPWmeanF/errordf_",
                                 storm_year,storm_name,"_",radius,"deg.csv"))
    }
  }
  }
  
  error.max <- max(abs(error_df$value))
  precip.max <- max(c(ST4_df$value,NAM_df$value))
  precip.min <- min(c(ST4_df$value,NAM_df$value))
  
  brk=seq(-1*ceiling(error.max),ceiling(error.max),ceiling(error.max/5))
  my.limits = c(-1*ceiling(error.max),ceiling(error.max))
  plot.edge <-7.2  #this is approx the number of degrees (100km) from the storm's eye
  precipcolors <- c("#FFFFFF", "#E0F8E0", "#81F781",
                    "#2EFE2E", "#C8FE2E", "#FFFF00",
                    "#FACC2E", "#FFBF00", "#FF8000", "#FF4000","#FF0000")
  regions<- c("alabama", "alaska", "arizona", "arkansas", "california", "colorado", 
              "connecticut", "delaware", "district of columbia", "florida", "georgia", 
              "hawaii", "idaho", "illinois", "indiana", "iowa", "kansas", "kentucky", 
              "louisiana", "maine", "maryland", "massachusetts", "michigan", "minnesota",
              "mississippi", "missouri", "montana", "nebraska", "nevada", "new hampshire", 
              "new jersey", "new mexico", "new york", "north carolina", "north dakota", 
              "ohio", "oklahoma", "oregon", "pennsylvania", "rhode island", "south carolina", 
              "south dakota", "tennessee", "texas", "utah", "vermont", "virginia", 
              "washington", "west virginia", "wisconsin", "wyoming")
  
  
  g1= ggplot(aes(x=x,y=y,fill=value),data=NAM_df) + 
    geom_tile() + theme_classic() + 
    geom_polygon(data=subset(map_data("state"), region %in% regions), 
                 aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
    scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max)) +
    labs(title =paste(storm_name, storm_date,"NAM"), x = "Longitude", y="Latitude") + 
    coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                       max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
                ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                       max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)
  
  
  g2= ggplot(aes(x=x,y=y,fill=value),data=ST4_df) + 
    geom_tile() + theme_classic() + 
    geom_polygon(data=subset(map_data("state"), region %in% regions), 
                 aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
    scale_fill_gradientn(colors = precipcolors ,na.value = "white",limits=c(0,precip.max)) + 
    labs(title =paste(paste(storm_name, storm_date,"ST4")),x = "Longitude", y="Latitude") + 
    coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                       max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
                ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                       max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)
  
  g3= ggplot(aes(x=x,y=y,fill=value),data=error_df) + 
    geom_tile() + theme_classic() + 
    geom_polygon(data=subset(map_data("state"), region %in% regions), 
                 aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
    scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white") +
    labs(title =paste(paste(storm_name, storm_date,"Error: Stage IV - NAM")),
         x = "Longitude", y="Latitude")+
    coord_fixed(xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                       max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
                ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
                       max(eye1_latlon[1],eye2_latlon[1])+plot.edge), ratio = 1)
  
  
  simul <- cbind(error_df[,2],error_df[,3],error_df[,1])
  simul.geo <- as.geodata(simul)
  
  simul.modvar <- variog(simul.geo, estimator.type="modulus")
  
  initial.values <- expand.grid(seq(0.01, 6, by=.02), seq(0.01, 4, by=.02))
  # initial.values <- expand.grid(seq(0, 6, by=.01), seq(0, 4, by=.01))
  # initial.values <- c(2.96,1.34)
  sim.mod.cressie     <- variofit(simul.modvar, ini.cov.pars = initial.values, 
                                  cov.model = "matern", fix.nugget = F, 
                                  weights = "cressie", fix.kappa = T, kappa = 1)
  
  SVGparamest[[i]] <- sim.mod.cressie
  
  print(sim.mod.cressie)
  tau2vec[i]  <- sim.mod.cressie$nugget
  sig2vec[i]  <- sim.mod.cressie$cov.pars[1]
  phivec[i]   <- sim.mod.cressie$cov.pars[2]
  prRangevec[i]<-sim.mod.cressie$practicalRange
  
  #   boot.mod.cressie <- boot.variofit(geodata = simul.geo, 
  #                                     obj.variog = simul.modvar,
  #                                     model.pars = sim.mod.cressie)
  #   
  #   SVGparamboot[[i]] <- boot.mod.cressie
  
  if(!pred){
    plot(simul.modvar, main=paste0("SVG for ", storm_year," ",storm_name),
         ylim=c(0,10));lines(sim.mod.cressie,col="orange")#;lines(sim.mod.cressie.nug,col="red")
  }
  
  multiplot(g1, g2, g3, cols=1)
}

svg.param.ests.error <- cbind(namevec, phivec, prRangevec, tau2vec, sig2vec)
if(write.pdf){
if(subtractPWmean) {
  write.csv(svg.param.ests.error, 
            paste0("~/NAM-Model-Validation/csv/svg.param.ests.error_deg_subtractPWmean_",
            if(pred){"pred_"},"sqrt.csv"))
} else {
  write.csv(svg.param.ests.error, 
            paste0("~/NAM-Model-Validation/csv/svg.param.ests.error_deg_noPWmean_",
            if(pred){"pred_"},"sqrt.csv"))
}
}

# save.image("LogPrecipVariograms.wks")

if(write.pdf){
  dev.off()
}

# Make the new raster summary files for pointwise (PW) mean, variance, std errors
if(makePWmean){
  raster_files <- list.files("~/NAM-Model-Validation/error_rasters_sqrt/", pattern = ".grd", full.names = T)
  error_sum <- mosaicList(raster_files)
  
  raster_files_sq <- list.files("~/NAM-Model-Validation/error_rasters_squared_sqrt/", pattern = ".grd", full.names = T)
  error_sum_sq <- mosaicList(raster_files_sq)
  
  error_files <- list.files("~/NAM-Model-Validation/error_rasters_counts_sqrt/", pattern = ".grd", full.names = T)
  error_counts <- mosaicList(error_files)
  
  if(write.pdf){
  if(!dir.exists("~/NAM-Model-Validation/error_rasters_summary_sqrt")) {
    dir.create("~/NAM-Model-Validation/error_rasters_summary_sqrt", recursive = T)
  }
  writeRaster(error_counts, "~/NAM-Model-Validation/error_rasters_summary_sqrt/error_counts", overwrite=T)
  writeRaster(error_sum, "~/NAM-Model-Validation/error_rasters_summary_sqrt/error_sum", overwrite=T)
  writeRaster(error_sum_sq, "~/NAM-Model-Validation/error_rasters_summary_sqrt/error_sum_sq", overwrite = T)
  }
  
  PW_mean <- error_sum/error_counts
  if(write.pdf){
  writeRaster(PW_mean, "~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_mean", overwrite=T)
  }
  plot(PW_mean, main="Pointwise Mean")
  PW_mean_spdf <- as((PW_mean), "SpatialPixelsDataFrame")
  PW_mean_df <- as.data.frame(PW_mean_spdf)
  colnames(PW_mean_df) <- c("value", "x", "y")
  
  # g4= ggplot(aes(x=x,y=y,fill=value),data=PW_mean_df) +
  #   geom_tile() + theme_classic() +
  #   geom_polygon(data=subset(map_data("state"), region %in% regions),
  #                aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0) +
  #   scale_fill_gradient2(low = "blue",mid = "white", high = "red",na.value = "white") +
  #   labs(title =paste("Pointwise Mean"))+
  #   coord_fixed(xlim=c(-110,-65),ylim=c(25,50), ratio = 1)
  # plot(g4) + scale_colour_manual(values = c("red","green","blue","black"), limits = c("4", "6", "8", "10"))
  
  plot(PW_mean, main="Pointwise Mean")
  plot(abs(PW_mean),main="Absolute Value of Pointwise Mean")
  S2 <- (error_sum_sq - (error_sum * error_sum)/error_counts)/(error_counts - 1)  
  if(write.pdf){
    writeRaster(S2, "~/NAM-Model-Validation/error_rasters_summary_sqrt/S2", overwrite=T)
  }
  plot(S2, main="Pointwise Variance Map")
  
  plot(PW_mean/sqrt(S2))
} else {
  # if not making the raster summaries, load them
  error_counts <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/error_counts.grd")
  error_sum    <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/error_sum.grd")
  error_sum_sq <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/error_sum_sq.grd")
}


# Plot pointwise means, variances and standard errors
par(mfrow=c(2,3))
pwms <- vars <- stds <- list()
for (n in c(1,6)) { 
  par(mar=c(5,4,4,2)+.1)
  # par(mfrow=c(1,2))
  Nmap <- error_counts
  Nmap[Nmap < n] <- NA
  # Nmap[Nmap >= n] <- 1
  
  pwms[[n]] <- PW_mean_geq5 <- error_sum/Nmap
  
  val <- max(abs(c(floor(4*minValue(pwms[[n]]))/4, ceiling(4*maxValue(pwms[[n]]))/4)))
  plot(pwms[[n]], main = paste0("PW mean when n >= ",n),
       col=c("blue",cm.colors(length(seq(-1.5,1.5,.25)[abs(seq(-1.5,1.5,.25)) < val ])-1),"red"),
       breaks=c(-val,seq(-1.5,1.5,.25)[abs(seq(-1.5,1.5,.25)) < val ],val))
  US(add=T, col="gray75")
  
  vars[[n]] <- S2_geq5 <- ((error_sum_sq - (error_sum*error_sum)/Nmap)/(Nmap-1))*mask.regrid
  plot(S2_geq5, main=paste0("Pointwise Variance when n >= ", n));US(add=T, col="gray49")
  
  stds[[n]] <- std_error_geq5 <- PW_mean_geq5/(sqrt(S2_geq5/Nmap))
  # png("stdgeq5.png",width = 1000, height = 1000) #cex.main = 2.5, cex.axis = 2
  val <- max(abs(c(floor(10*minValue(stds[[n]]))/10, ceiling(10*maxValue(stds[[n]]))/10)))
  plot(stds[[n]], main = paste0("Std err when n >= ",n),
       col=c("blue",cm.colors(length(seq(-3,3,.25))-1),"red"),
       breaks=c(-val,seq(-3,3,.25),val))
  US(add=T, col="gray75")
}

# for (i in 1:n) {
#   plot(pwms[[i]], main="Pointwise Mean",
#        col=c("blue",cm.colors(5),"red"), breaks=c(-2,-.8,-.4,-.2,.2,.4,.8,2))
#   US(add=T, col="gray75")
# }