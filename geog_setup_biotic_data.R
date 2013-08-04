setwd('~/maxent/geog/data')

library(sp)
library(rgdal)

## PART I - setup biotic spatial coordinates
## read in biotic coordinates in and aggregate into a list
## in data col1 is lat and col2 is long, switch this order on read in

datasets = c('bbs','cbc','fia','gentry','mcdb','naba')
coordFileEnding = '_lat_long.csv'
parmFileEnding = '_dist_test.csv'

meteData = list()

for(i in seq_along(datasets)){
  columns = 2:1
  if(datasets[i] == 'naba'){
    meteData[[i]] = cbind(
      read.table(paste(datasets[i],coordFileEnding,sep=''),
                 header=FALSE,sep=',')[,columns],
      read.table(paste(datasets[i],parmFileEnding,sep=''),
                 header=FALSE,sep=',',quote=''))  # requires the quote argument
  }
  else{
    meteData[[i]] = cbind(
      read.table(paste(datasets[i],coordFileEnding,sep=''),
                 header=FALSE,sep=',')[,columns],
      read.table(paste(datasets[i],parmFileEnding,sep=''),
                header=FALSE,sep=','))
  }
  colnames(meteData[[i]]) = c('Longitude','Latitude','SiteID','year','S','N',
                               'p1','p2','p3','p4')
}

names(meteData) = datasets

## for the time being fix this datapoint which appears to have lost its decimal palce
meteData$mcdb[16,1] = -114.97

## examine spatial domain of North American datasets
range(meteData$bbs$Longitude,meteData$cbc$Longitude,meteData$fia$Longitude,
      meteData$naba$Longitude)
range(meteData$bbs$Latitude,meteData$cbc$Latitude,meteData$fia$Latitude,
      meteData$nabaLatitude)

studyExtentXvals = c(-167,-167,-50,-50,-167)
studyExtentYvals = c(19,70,70,19,19)

meteData = lapply(meteData, function(x){
             x = x[as.logical(point.in.polygon(x$Longitude,x$Latitude,
                 studyExtentXvals,studyExtentYvals)),] 
             x
             }
           )
 
## convert meteData into a list of SpatialPointsDataFrames
for(i in seq_along(datasets)){
  meteData[[i]] = SpatialPointsDataFrame(coords = meteData[[i]][,1:2],
                  data = meteData[[i]],
                  proj4string = CRS('+proj=longlat +ellps=WGS84'))
}

gc()

## write data product to file
save(meteData,file="meteData.Rdata")

## Optional: take a look at geographic distribution of biotic data
cls = c('green3', 'dodgerblue', 'red', 'purple', 'orange', 'yellow')
names_ordered = c('fia', 'bbs', 'cbc', 'naba', 'mcdb', 'gentry')
library(maps)
pdf('../figs/maps_of_data.pdf')
par(mfrow=c(1,1))
map('world',c('canada','usa','mexico'),xlim=c(-170,-50))
for(i in 1:length(meteData))
  points(meteData[[names_ordered[i]]],pch=19, cex=.5, col=cls[i])
## bounding box corners
#points(c(-157,-52),c(24,67),col='red',pch=19,cex=3) 
for(i in 1:length(meteData)) {
  map('world',c('canada','usa','mexico'),xlim=c(-170,-50))
  for(j in 1:i) {
    if (j < i)
      points(meteData[[names_ordered[j]]],pch=19,cex=.5,col='grey')
    else
      points(meteData[[names_ordered[j]]],pch=19,cex=.5,col=cls[j])
  }
}
dev.off()


