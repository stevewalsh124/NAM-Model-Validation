library(dplyr) # %>%, filter
library(stringr) #str_replace_all
#5th column is model we want "  NAM" (note the spaces)
#11th column we want " HU" or " TS" for hurricane/tropical storm (note the space)

# # if(ncol(track_file)==40){
# #   colnames(track_file) <- c("BASIN","CY","DATE","TECHNUM/MIN","TECH","TAU","LatN/S","LonE/W",
# #                             "VMAX","MSLP","TY","RAD","WINDCODE","RAD1","RAD2","RAD3","RAD4",
# #                             "POUTER","ROUTER","RMW","GUSTS","EYE","SUBREGION","MAXSEAS","INITIALS",
# #                             "DIR","SPEED","STORMNAME","DEPTH","SEAS", #uptohere for 30 cols
# #                             "SEASCODE","SEAS1","SEAS2",
# #                             "SEAS3","SEAS4","USER1","USER2","USER3","USER4","USER5")
# #   
# # }else{print("the number of columns is not 30 or 40")}


NAM_files <- c() #tells which file numbers had NAM info

allNAM <- data.frame() #contains all lines from files where MODEL ID is "  NAM"
allNAM2<- data.frame()
everyNAMNHC <- data.frame()

## GET NHC INFO
#loop through all files
all_track_folders <- list.dirs("NAM_NHC_track")[-1]
for(f in 1:length(all_track_folders)){
  print(all_track_folders[f])
  for (i in 1:length(list.files(all_track_folders[f], full.names = T))) {
    # print(i)
    track_file <- read.csv(list.files(all_track_folders[f], full.names = T)[i], colClasses = "character" ,header = F)

    #indicates if a file contained NAM info
    ifelse("  NAM" %in% levels(track_file$V5),
           NAM_files[length(NAM_files)+1] <- 1,
           NAM_files[length(NAM_files)+1] <- 0)

    #checks to see if there are any lines where Model ID is "  NAM" in $V5
    track_NAM <- (track_file %>% filter(V5=="  NAM"))

    if(nrow(track_NAM)!=0){
      print(i)
      print(head(track_NAM[,1:15]))
      allNAM <- merge(allNAM, track_NAM, all = T)
      if(length(allNAM2)==0){allNAM2 <- track_NAM}
      allNAM2 <- full_join(allNAM2,track_NAM)
    }else{print("no lines with NAM info")}
  }
}

#resulting allNAM has 43025 rows; 43048 rows
dim(allNAM);dim(allNAM2)

which(track_file$V5=="  NAM")
allNAM_12_24_step <- allNAM[which(allNAM$V5=="  NAM" &
                                    (grepl("12",allNAM$V6)|grepl("24",allNAM$V6))&
                                    !grepl("6",allNAM$V6)&!grepl("0",allNAM$V6)&
                                    !grepl("3",allNAM$V6)&!grepl("-",allNAM$V6)),]

allNAM2_12_24_step <- allNAM2[which(allNAM2$V5=="  NAM" &
                                    (grepl("12",allNAM2$V6)|grepl("24",allNAM2$V6))&
                                    !grepl("6",allNAM2$V6)&!grepl("0",allNAM2$V6)&
                                    !grepl("3",allNAM2$V6)&!grepl("-",allNAM2$V6)),]
#another apporach to getting all 
all.equal(allNAM_12_24_step,allNAM2_12_24_step)
all.equal(allNAM_12_24_step[,6], filter(allNAM,(V6==24|V6==12|V6=="  12"|V6=="  24"))[,6])
all.equal(allNAM2_12_24_step[,6], filter(allNAM2,(V6==24|V6==12|V6=="  12"|V6=="  24"))[,6])

allNAM_12_24_step  #only 6633 out of the 43025 has timestep 12 or 24
allNAM2_12_24_step #only 6639 out of the 43048 has timestep 12 or 24

NHC_info <- c() #collects date time and timestep for each line of NHC info
for (i in 1:dim(allNAM_12_24_step)[1]) {
  ifelse(nchar(str_remove_all(as.character(allNAM_12_24_step$V6[i]), " "))==1, 
         temp <- paste0(0,allNAM_12_24_step$V6[i]),
         temp <- allNAM_12_24_step$V6[i])
  NHC_info[i] <- paste0(allNAM_12_24_step$V3[i],temp)
}
NHC_info <- str_replace_all(NHC_info," ","") #remove random spaces
length(NHC_info) #6633 (6639) lines of files where model ID was "  NAM" and forecast timestep was 12 or 24

length(unique(NHC_info)) #4947 (4953)







# NAM_rows <- c()
# for (i in 1:length(track_file$V5)) {
#   if(grepl("NAM",track_file$V5[i])){
#       NAM_rows[i] <- track_file$V5[i] #find out if ith row model ID contains substr "NAM"
#   }
# }
# 
# sum(NAM_rows, na.rm = T)/length(track_file$V5) #what percent of rows contain "NAM" in model ID
# model_IDS_wNAM <- levels(track_file$V5)[unique(NAM_rows)]   #get the unique model IDs containing "NAM"
# model_IDS_wNAM[!is.na(model_IDS_wNAM)] # remove NA








all.storm.folders <- list.dirs("NAMandST4", recursive = F) #/Volumes/LACIEHD/
most.storm.folders <- list.dirs("NAMandST4", recursive = F)[1:47] #/Volumes/LACIEHD/

NAM_info <- matrix(nrow = 2*length(most.storm.folders), ncol = 3)
NAM_NHC_lines <- c()

for(i in 1:length(most.storm.folders)){
  # i <- 30
  currentST4_NAM_folders <- list.dirs(most.storm.folders[i])[-1] #[-1] avoids parent folder
  ST4_folder <- currentST4_NAM_folders[1]
  NAM_folder <- currentST4_NAM_folders[2] #uncomment the [2] as well
  name_start <- tail(str_locate_all(pattern ='/', ST4_folder)[[1]],1)
  storm_name <- substr(ST4_folder,name_start[1]+5,nchar(ST4_folder))
  
  #Identifying where in the NAM folder the first 24 hours after landfall files are (latest timestep)
  if(length(list.files(NAM_folder,full.names = T))==18){first24 <- c(7,8)} else {
    if(length(list.files(NAM_folder,full.names = T))==15|storm_name=="rita"){first24 <- c(5,6)} else {
      if(length(list.files(NAM_folder,full.names = T))==26){first24 <- c(5,9)} #eg: 2004frances, storm3
      else {print("Total NAM files for this storm not 15 or 18; skipping"); next}}} #storm 3 for example
  
  
  for (j in 1:length(list.files(NAM_folder,full.names = T)[first24])) {
    print(c(i,j))# j <-1
    NAMfile <- list.files(NAM_folder,full.names = T)[first24[j]]
    
    date_start <- str_locate_all(pattern ='_218_', NAMfile)[[1]][2,2] + 1
    grib_start <- str_locate_all(pattern ='.gr', NAMfile)[[1]][1,1]
    forec_date <- substr(NAMfile,date_start,date_start+7)
    forec_time <- substr(NAMfile,date_start+9,date_start+10)
    forec_step <- substr(NAMfile,grib_start-2,grib_start-1)
    forec_info <- paste0(forec_date,forec_time, forec_step)
    print(c(2*(i-1)+j,storm_name, forec_info, forec_info %in% NHC_info))
    NAM_info[2*(i-1)+j,] <- c(storm_name, forec_info, forec_info %in% NHC_info)
    if(forec_info %in% NHC_info){
      if(length(NAM_NHC_lines)==0){NAM_NHC_lines <- allNAM_12_24_step[grepl(forec_info, NHC_info),]}
      NAM_NHC_lines <- full_join(NAM_NHC_lines, 
                             allNAM_12_24_step[grepl(forec_info, NHC_info),])
    }
  }
}

dim(NAM_NHC_lines)
NAM_NHC_lines[order(NAM_NHC_lines$V3),1:8]
# write.csv(NAM_NHC_lines,file = "NAM_NHC_lines_best_track.csv")
colnames(NAM_info) <- c("storm_name","'YYYYMMDDHHTS'","info_in_NHC?")
head(NAM_info)

new_false <- NAM_info[which(NAM_info[,3]=="FALSE")]
new_true <- NAM_info[which(NAM_info[,3]=="TRUE")]


if(NAM_info[i,3]=="TRUE"){}#use points to add centers to NAM plotter}

#ran logprecipvariograms.Rmd, then 

plot(NAM_plotter,xlim=c(min(eye1_latlon[2],eye2_latlon[2])-plot.edge, 
                       max(eye1_latlon[2],eye2_latlon[2])+plot.edge),
    ylim=c(min(eye1_latlon[1],eye2_latlon[1])-plot.edge, 
           max(eye1_latlon[1],eye2_latlon[1])+plot.edge))

d <-NAM_BT[grepl(paste0(storm_date,storm_hour_NAM),as.character(NAM_BT[,4])),]
for (i in 1:nrow(d)) {
  points(as.numeric(substr(d[i,9],1,nchar(as.character(d[i,9]))-1))/-10, 
         as.numeric(substr(d[i,8],1,nchar(as.character(d[i,8]))-1))/ 10, pch=16)
}
d
