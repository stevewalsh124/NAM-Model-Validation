library(dplyr)
track_files <- list.files("NAM_NHC_track/NAMbesttrack2015/", full.names = T)
#5th column is model we want "  NAM" (note the spaces)
#11th column we want " HU" or " TS" for hurricane/tropical storm (note the space)

HU_TS_files <- c()
NAM_files <- c()
NAM_HU_files <- c()

allNAM <- data.frame()
everyNAMNHC <- data.frame()

#loop through all files
all_track_folders <- list.dirs("NAM_NHC_track")
for(f in 2:7){#length(all_track_folders)){
  print(all_track_folders[f])
  for (i in 1:length(list.files(all_track_folders[f], full.names = T))) {
    # print(i)
    track_file <- read.csv(list.files(all_track_folders[f], full.names = T)[i], header = F)
    print(unique(track_file$V5))
    # # if(ncol(track_file)==30){
    # #   colnames(track_file) <- c("BASIN","CY","DATE","TECHNUM/MIN","TECH","TAU","LatN/S","LonE/W",
    # #                             "VMAX","MSLP","TY","RAD","WINDCODE","RAD1","RAD2","RAD3","RAD4",
    # #                             "POUTER","ROUTER","RMW","GUSTS","EYE","SUBREGION","MAXSEAS","INITIALS",
    # #                             "DIR","SPEED","STORMNAME","DEPTH","SEAS")
    # #   
    # # } 
    # # if(ncol(track_file)==40){
    # #   colnames(track_file) <- c("BASIN","CY","DATE","TECHNUM/MIN","TECH","TAU","LatN/S","LonE/W",
    # #                             "VMAX","MSLP","TY","RAD","WINDCODE","RAD1","RAD2","RAD3","RAD4",
    # #                             "POUTER","ROUTER","RMW","GUSTS","EYE","SUBREGION","MAXSEAS","INITIALS",
    # #                             "DIR","SPEED","STORMNAME","DEPTH","SEAS"  ,"SEASCODE","SEAS1","SEAS2",
    # #                             "SEAS3","SEAS4","USER1","USER2","USER3","USER4","USER5")
    # #   
    # # }else{print("the number of columns is not 30 or 40")}
    # 
    # #indicates if a file contained NAM info
    # if("  NAM" %in% levels(track_file$V5)){NAM_files[length(NAM_files)+1] <- 1}
    # else{NAM_files[length(NAM_files)+1] <- 0}
    # 
    # #indicates if a file contained HU/TS info
    # if(" HU" %in% levels(track_file$V11) | " TS" %in% levels(track_file$V11)){HU_TS_files[length(HU_TS_files)+1] <- 1}
    # else{HU_TS_files[length(HU_TS_files)+1] <- 0}
    # 
    # #checks to see if there are any lines where the info is "  NAM" *AND* (" HU" *OR* " TS")
    # track_NAM_HU <- (track_file %>%
    #                    # filter(V11 %in% c(" HU"," TD"," TS")) %>% #NAM info not marked TS/HU appropriately
    #                    filter(V5=="  NAM"))
    # 
    # if(nrow(track_NAM_HU)!=0){
    #   print(i)
    #   print(head(track_NAM_HU))
    #   NAM_HU_files[length(NAM_HU_files)+1] <- 1
    #   allNAM <- merge(allNAM, track_NAM_HU, all = T)
    # }
    # else{print("no lines with NAM and HUR info"); NAM_HU_files[length(NAM_HU_files)+1] <- 0}
    everyNAMNHC <- merge(everyNAMNHC,track_file,all=T)
  }
  
}

cbind(HU_TS_files,NAM_files, NAM_HU_files)

which(track_file$V5=="  NAM")
which(track_file$V11==" HU")
which(track_file$V5=="  NAM" & track_file$V11==" HU")

allNAM
sort(unique(allNAM$V3))

substr(track_file$V3,1,6)
grepl(track_file$V5, "NAM")

NAM_rows <- c()
for (i in 1:length(track_file$V5)) {
  if(grepl("NAM",track_file$V5[i])){
      NAM_rows[i] <- track_file$V5[i] #find out if ith row model ID contains substr "NAM"
  }
}

sum(NAM_rows, na.rm = T)/length(track_file$V5) #what percent of rows contain "NAM" in model ID
d <- levels(track_file$V5)[unique(NAM_rows)]   #get the unique model IDs containing "NAM"
d[!is.na(d)]


all.storm.folders <- list.dirs("/Volumes/LACIEHD/NAMandST4", recursive = F) #/Volumes/LACIEHD/
most.storm.folders <- list.dirs("/Volumes/LACIEHD/NAMandST4", recursive = F)[[c(1,10,40)]] #/Volumes/LACIEHD/

for(i in 1:length(most.storm.folders)){
  # i <- 1
  currentST4_NAM_folders <- list.dirs(most.storm.folders[i])[-1] #[-1] avoids parent folder
  ST4_folder <- currentST4_NAM_folders[1]
  NAM_folder <- currentST4_NAM_folders[2] #uncomment the [2] as well
  
  
  for (j in 1:length(list.files(NAM_folder,full.names = T)[first24])) {
    # j <-1
    NAMfile <- list.files(NAM_folder,full.names = T)[first24[j]]
    
    
    date_start <- str_locate_all(pattern ='_218_', NAMfile)[[1]][2,2] + 1
    storm_date <- as.numeric(substr(NAMfile,date_start,date_start+7))
    forec_date <- as.numeric(substr(NAMfile,date_start,date_start+7))
    forec_step <- as.numeric(substr(NAMfile,date_start+9,date_start+10))
    forec_time <- as.integer(paste0(forec_date,forec_step))
    print(forec_time %in% allNAM$V3)
  }
}
