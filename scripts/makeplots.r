# some code to make some image plots in R with minimal whitespace so
# we can jam alot of them on a single plot.

 # storm type
landfall = read.table('csv/landfall_locns.csv',sep=',',header=T,stringsAsFactors=FALSE)$x
lf = paste(1:47,substr(landfall,1,1))

 # read in data
load('RData/TC_array.RData')
dim(TC_array)

 # read in the coordinate data
load('RData/coords.RData')

 # get the uniqe coordinates for lat and lon
longitude = sort(unique(coords[,1]))
latitude = sort(unique(coords[,2]))

 # determine coordinates for the image
ilon = match(coords[,1],longitude)
jlat = match(coords[,2],latitude)

 # make a 47x3x286x214 array so the last to dimensions give a matrix 
tcdata = array(TC_array,c(47,3,286,214))

 # now, loop thru the images and fill it with Steve's data
for(i in 1:47) for(j in 1:3){
 tcdata[i,j,,] = NA
 for(k in 1:nrow(coords)) tcdata[i,j,ilon[k],jlat[k]] = TC_array[i,j,k]
}

 # make sure this looks ok
image(longitude,latitude,tcdata[1,1,,])
plot(ilon,jlat,pch='.')
points(ilon,jlat,pch='.',col=ifelse(is.na(TC_array[1,1,]),'yellow','red'))

 # make 0's in rainfall NA for plotting purposes
for(i in 1:47){
  tcdata[i,1,,] = ifelse(tcdata[i,1,,] <= 0,NA,tcdata[i,1,,])
  tcdata[i,2,,] = ifelse(tcdata[i,2,,] <= 0,NA,tcdata[i,2,,])
  x = ifelse(is.na(tcdata[i,2,,]) | is.na(tcdata[i,1,,]),NA,tcdata[i,3,,])
  tcdata[i,3,,] = ifelse(is.na(tcdata[i,2,,]) | is.na(tcdata[i,1,,]) |
                         (tcdata[i,2,,]==0 & tcdata[i,1,,]==0),NA,tcdata[i,3,,])
}

 # add in the maps library
library(maps)
library(RColorBrewer)
data(usaMapEnv)

 # mess with color maps for the image
imcols = terrain.colors(20)
imcols = topo.colors(20)[20:1]
imcols = heat.colors(20)[20:1]
imcols = brewer.pal(9, 'Blues')
brewer.pal(20, 'BuPu')
imcols = colorRampPalette(c('white','blue','yellow','gold'))(35)
 # here's my choice
whiteBlueCols = colorRampPalette(c('white','blue'))(20)
blueYellowCols = colorRampPalette(c('blue','yellow'))(10)
imcols = c(whiteBlueCols,blueYellowCols[-1])
 # for differences
blueYellowCols2 = colorRampPalette(c('blue','yellow'))(16)
yellowRedCols2 = colorRampPalette(c('yellow','red'))(16)
imcolsDiff = c(blueYellowCols2,yellowRedCols2[-1])

 # make a plot of the 42 observations
pdf('pdf/allStorms.pdf',width=8,height=6)
par(mfrow=c(13,12),oma=c(4,4,1,1),mar=c(0,0,0,0))
for(i in 1:47){
  image(longitude,latitude,tcdata[i,1,,],zlim=range(tcdata[1:5,1:2,,],na.rm=T),axes=F,col=imcols)
  box(); text(-97,42,lf[i])
  map('usa',add=T)
  image(longitude,latitude,tcdata[i,2,,],zlim=range(tcdata[1:5,1:2,,],na.rm=T),axes=F,col=imcols)
  box(); text(-97,42,lf[i])
  map('usa',add=T)
  image(longitude,latitude,sqrt(abs(tcdata[i,3,,]))*sign(tcdata[i,3,,]),zlim=c(-3.7,3.7),axes=F,col=imcolsDiff)
  box(); text(-97,42,lf[i])
  map('usa',add=T)
}


 # add in some legends
par(cex=0.7)  #, mai=c(0.1,0.1,0.2,0.1))
par(fig=c(0.15,0.35,0.00,0.03), new=TRUE)
precipRange = range(tcdata[1:5,1:2,,],na.rm=T)
image(x=seq(0,18.1,length=31),y=0:1,matrix(seq(0,18.1,length=30),nrow=30,ncol=1),
      col=imcols,axes=F)
box(); axis(1,at=c(0,5,10,15))
mtext('precip (mm)',side=1,line=2.5)

par(fig=c(0.65,0.85,0.00,0.03), new=TRUE)
precipRange = range(tcdata[1:5,1:2,,],na.rm=T)
image(x=seq(-3.7,3.7,length=31),y=0,matrix(seq(-3.7,3.7,length=31),nrow=31,ncol=1),
      col=imcolsDiff,axes=F)
box(); axis(1,at=c(-sqrt(10),-sqrt(4),0,sqrt(4),sqrt(10)),labels=c(-10,-4,0,4,10))
mtext('precip (mm)',side=1,line=2.5)
dev.off()

 # now a starter for the paper figure - let's try 6 of each:
 # pick gulf storms to show
igulf = c(7,9,13,20,40,47)
ifla = c(1,2,14,19,29,46)
iatl = c(4,27,31,36,38,43)

 # make the plots
pdf('pdf/someStorms.pdf',width=6,height=6.5)
par(mfrow=c(10,6),oma=c(1.5,2,1.5,1),mar=c(0,0,0,0))
for(ii in 1:6){
  i = igulf[ii]
  image(longitude,latitude,tcdata[i,1,,],zlim=range(tcdata[1:5,1:2,,],na.rm=T),axes=F,col=imcols)
  box(); #text(-97,42,lf[i])
  if(ii<=2) mtext('obseration',side=3,line=.3,cex=.7)
  map('usa',add=T)
  image(longitude,latitude,tcdata[i,2,,],zlim=range(tcdata[1:5,1:2,,],na.rm=T),axes=F,col=imcols)
  box(); #text(-97,42,lf[i])
  if(ii<=2) mtext('prediction',side=3,line=.3,cex=.7)
  map('usa',add=T)
  image(longitude,latitude,sqrt(abs(tcdata[i,3,,]))*sign(tcdata[i,3,,]),zlim=c(-3.7,3.7),axes=F,col=imcolsDiff)
  box(); #text(-97,42,lf[i])
  if(ii<=2) mtext('obs - pred',side=3,line=.3,cex=.7)
  map('usa',add=T)
}
for(ii in 1:6){
  i = ifla[ii]
  image(longitude,latitude,tcdata[i,1,,],zlim=range(tcdata[1:5,1:2,,],na.rm=T),axes=F,col=imcols)
  box(); #text(-97,42,lf[i])
  map('usa',add=T)
  image(longitude,latitude,tcdata[i,2,,],zlim=range(tcdata[1:5,1:2,,],na.rm=T),axes=F,col=imcols)
  box(); #text(-97,42,lf[i])
  map('usa',add=T)
  image(longitude,latitude,sqrt(abs(tcdata[i,3,,]))*sign(tcdata[i,3,,]),zlim=c(-3.7,3.7),axes=F,col=imcolsDiff)
  box(); #text(-97,42,lf[i])
  map('usa',add=T)
}
for(ii in 1:6){
  i = iatl[ii]
  image(longitude,latitude,tcdata[i,1,,],zlim=range(tcdata[1:5,1:2,,],na.rm=T),axes=F,col=imcols)
  box(); #text(-97,42,lf[i])
  map('usa',add=T)
  image(longitude,latitude,tcdata[i,2,,],zlim=range(tcdata[1:5,1:2,,],na.rm=T),axes=F,col=imcols)
  box(); #text(-97,42,lf[i])
  map('usa',add=T)
  image(longitude,latitude,sqrt(abs(tcdata[i,3,,]))*sign(tcdata[i,3,,]),zlim=c(-3.7,3.7),axes=F,col=imcolsDiff)
  box(); #text(-97,42,lf[i])
  map('usa',add=T)
}
mtext('Gulf landfall',side=2,line=.5,outer=T,at=.85)
mtext('Florida landfall',side=2,line=.5,outer=T,at=.55)
mtext('Atlantic landfall',side=2,line=.5,outer=T,at=.25)


 # add in some legends
par(cex=0.7)  #, mai=c(0.1,0.1,0.2,0.1))
par(fig=c(0.0,0.2,0.06,0.09), new=TRUE)
precipRange = range(tcdata[1:5,1:2,,],na.rm=T)
image(x=seq(0,18.1,length=31),y=0:1,matrix(seq(0,18.1,length=30),nrow=30,ncol=1),
      col=imcols,axes=F)
box(); axis(1,at=c(0,5,10,15),cex=.5)
mtext('precip (mm)',side=1,line=2.0,cex=.7)

par(fig=c(0.33,.53,0.06,0.09), new=TRUE)
precipRange = range(tcdata[1:5,1:2,,],na.rm=T)
image(x=seq(-3.7,3.7,length=31),y=0,matrix(seq(-3.7,3.7,length=31),nrow=31,ncol=1),
      col=imcolsDiff,axes=F)
box(); axis(1,at=c(-sqrt(12),-sqrt(8),-sqrt(4),0,sqrt(4),sqrt(8),sqrt(12)),labels=c('-12','','-4','0','4','8','12'),cex=.5)
mtext(paste(expression(Delta),'precip (mm)'),side=1,line=2.0,cex=.7)

dev.off()


