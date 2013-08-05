## this file should be run from the ./maxent/geog/ directory

setwd('~/maxent/geog/')

i = as.integer(commandArgs()[length(commandArgs())])

library(raster)
load('./data/meteData.Rdata')
load('./data/sampleCircles.Rdata')
load('~/plant_evol/gbif/wwfeco.Rdata')
enviLayerslaea = stack('./data/enviLayerslaea.grd')
bioDat = meteData[[i]]

## Sample mean and variance inside circles
enviAvg = extract(enviLayerslaea, sampleCircles[[i]], fun=mean, na.rm=T,
           layer=1, nl=nlayers(enviLayerslaea), progress='text')
colnames(enviAvg) = sapply(colnames(enviAvg),
                     function(x) paste(x, '.mean', sep=''))
enviVar = extract(enviLayerslaea,sampleCircles[[i]], fun=var, na.rm=T,
           layer=1, nl=nlayers(enviLayerslaea), progress='text')
colnames(enviVar) = sapply(colnames(enviVar),
                     function(x) paste(x, '.var', sep=''))

## Sample biome at the coordinates of the atlas sample
biome = overlay(wwfeco,SpatialPoints(bioDat))$eco_code

## Output data product
nameDat = names(meteData)[i]
outDat = data.frame(bioDat, enviAvg, enviVar, biome)
write.csv(outDat, file=paste('~/maxent/geog/data/',
                             nameDat, '_envidat.csv', sep=''), row.names=F)

