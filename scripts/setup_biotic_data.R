
library(sp)
library(rgdal)
library(maps)

## PART I - setup biotic spatial coordinates
## read in biotic coordinates in and aggregate into a list
## in data col1 is lat and col2 is long, switch this order on read in

datasets = c('bbs_2012', 'bbs_2008_2012', 'cbc','fia','gentry','naba')
siteFileEnding = '_coords.csv'
spFileEnding = '_spab.csv'

mete_data = vector('list', length(datasets))
names(mete_data) = datasets

for (i in seq_along(datasets)) {
    if (grepl('bbs', datasets[i]))
        coords = read.csv('./data/raw/bbs_coords.csv')
    else
        coords = read.csv(paste('./data/raw/', datasets[i],
                                '_coords.csv', sep=''))
    spab = read.csv(paste('./data/raw/', datasets[i], '_spab.csv', sep=''),
                    na.string='\\N')
    ## remove rows with NA in abundance column
    spab = spab[!is.na(spab$ab), ]
    sr = tapply(spab$sp_id, spab$site_id, function(x) length(unique(x)))
    ab = tapply(spab$ab, spab$site_id, sum)
    ## match up the site ids so rows of coords match
    coords = coords[match(names(sr), coords$site_id), ]
    mete_data[[i]] = cbind(coords, sr, ab)
    colnames(mete_data[[i]]) = c(names(coords), 'sr', 'ab')
}

pdf('./figs/white_et_al_2012_data_map.pdf')
cls = c('green3', 'dodgerblue', 'red', 'mediumpurple', 'darkgreen')
names_ordered = c('fia', 'bbs_2012', 'cbc', 'naba', 'gentry')
map('world', interior=F)
for(i in seq_along(mete_data))
  points(mete_data[[names_ordered[i]]]$longitude,
         mete_data[[names_ordered[i]]]$latitude,
         pch=19, cex=.5, col=cls[i])
dev.off()

## examine spatial domain of North American datasets
range(mete_data$bbs_2012$longitude,mete_data$cbc$longitude,mete_data$fia$longitude,
      mete_data$naba$longitude)
range(mete_data$bbs_2012$latitude,mete_data$cbc$latitude,mete_data$fia$latitude,
      mete_data$nabalatitude)

extent_x = c(-167,-167,-40,-40,-167)
extent_y = c(-30,70,70,-30,-30)

mete_data = lapply(mete_data, function(x){
                x = x[as.logical(point.in.polygon(
                                 x$longitude, x$latitude,
                                 extent_x, extent_y)), ] 
                return(x)
             }
           )
 
## convert mete_data into a list of SpatialPointsDataFrames
for(i in seq_along(datasets)) {
    mete_data[[i]] = SpatialPointsDataFrame(
                         coords = mete_data[[i]][ , c('longitude', 'latitude')],
                         data = mete_data[[i]], 
                         proj4string = CRS('+proj=longlat +ellps=WGS84'))
}

gc()

## write data product to file
save(mete_data, file="./data/mete_data.Rdata")
#load('./data/mete_data.Rdata')

for (i in seq_along(mete_data)) {
    write.csv(mete_data[[i]]@data[ , c('site_id', 'longitude','latitude')],
              file=paste('./data/', names(mete_data)[i], '_coords.csv',
                         sep=''), row.names=F)
}

## Optional: take a look at geographic distribution of biotic data
library(maps)
pdf('./figs/maps_of_data.pdf')
par(mfrow=c(1,1))
cls = c('green3', 'dodgerblue', 'dodgerblue', 'red', 'mediumpurple', 'darkgreen')
names_ordered = c('fia', 'bbs_2012', 'bbs_2008_2012',
                  'cbc', 'naba', 'gentry')
map('world', xlim=c(-170,-40), ylim=c(-40,70))
for(i in 1:length(mete_data))
  points(mete_data[[names_ordered[i]]]$longitude,
         mete_data[[names_ordered[i]]]$latitude,
         pch=19, cex=.5, col=cls[i])
## bounding box corners
#points(c(-157,-52),c(24,67),col='red',pch=19,cex=3) 
for(i in 1:length(mete_data)) {
    map('world',xlim=c(-170,-40), ylim=c(-40,70))
    for(j in 1:i) {
        if (j < i)
            points(mete_data[[names_ordered[i]]]$longitude,
                   mete_data[[names_ordered[i]]]$latitude,
                   pch=19,cex=.5,col='grey')
        else
            points(mete_data[[names_ordered[i]]]$longitude,
                   mete_data[[names_ordered[i]]]$latitude,
                   pch=19,cex=.5,col=cls[i])
    }
}
dev.off()


