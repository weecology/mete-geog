i = as.integer(commandArgs()[6])
sink(file = paste('geog_calc_envi_out_',i,'.txt',sep=''))

library(raster)
load('meteData.Rdata')
load('sampleCircles.Rdata')
load('/home/danmcglinn/GIS/WWFecoregionbiome/wwfeco.Rdata')
enviLayerslaea = stack('enviLayerslaea.grd')
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
biome = overlay(wwfeco,SpatialPoints(bioDat))[,19]

## Output data product
nameDat = names(meteData)[i]
outDat = data.frame(bioDat, enviAvg, enviVar, biome)
write.csv(outDat, file=paste(nameDat, '_envidat2.csv', sep=''), row.names=F)

sink()
unlink(paste('geog_calc_envi_out_',i,'.txt',sep=''))
