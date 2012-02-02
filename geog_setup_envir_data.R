library(sp)
library(raster)
library(rgdal)

startDir = '/home/danmcglinn/maxent/geog/'
setwd(startDir)
## PART II - query climate and ecosystem data

## set study exent
box = c(-167,-50,19,70) 
myExtent = extent(box)
myRes = 0.1
wgsPrjString <- "+proj=longlat +ellps=WGS84 "
laeaPrjString <- "+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"

rasTemplate = raster(myExtent,crs=wgsPrjString)
res(rasTemplate) = myRes

## Read in raster datasets (Climate and NDVI) and crop to myExtent
## load bioclim datasets
## 5 arc-minutes resolution
setwd('/home/danmcglinn/GIS/WorldClimData')
#setwd('C:/Users/White Lab/Documents/Lab data/GIS/WorldClim data')
load('bioclim_5m.Rdata') 
bioStack = crop(bioStack, myExtent)
## load worldclim elev data
## 5 arc-minutes resolution
load('alt_5m.Rdata')
alt = crop(alt,myExtent)

## load NDVI data
setwd('/home/danmcglinn/GIS/GlobalVegData/Mean&Std')
#setwd('C:/Users/White Lab/Documents/Lab data/GIS/GlobalVegData/Mean&Std')
ndviJun = stack('JUNAV18.rst') 
ndviDec = stack('decav18.rst') 
ndvi = stack(ndviJun,ndviDec)
ndvi = crop(ndvi,myExtent)
fixNDVI = function(x){
  x2 = ifelse(x ==1,NA,x)
  (x2-128)/128
}
ndvi = calc(ndvi,fixNDVI)

## resample WorldClim layers to resolution of NDVI before stacking them
bioTemp = resample(bioStack, rasTemplate, method='bilinear',progress='text')
bioStack = bioTemp
altTemp = resample(alt, rasTemplate, method='bilinear',progress='text')
alt = altTemp
rm(altTemp, bioTemp)
gc()

## stack all the layers
enviLayers = stack(alt,ndvi,bioStack)
enviLayers@layernames = c('alt','ndviJun','ndviDec',
                        bioStack@layernames)

## reproject to lambert equal area
temp = enviLayers[[1]]
temp.ext = projectExtent(temp, CRS(laeaPrjString))
enviLayerslaea = projectRaster(from=enviLayers, to=temp.ext,
                 crs=CRS(laeaPrjString),progress='text')
enviLayerslaea@layernames = enviLayers@layernames

setwd(startDir)
writeRaster(enviLayerslaea,file='./data/enviLayerslaea.grd',bandorder='BIL')
#enviLayerslaea = stack('enviLayerslaea.grd')

rm(temp, temp.ext, alt, bioStack, enviLayers, ndviDec, ndviJun, wwfeco)
gc()

## Draw sample circles around atlas locations
makeCir = function(p,r){
  # A function that draws a circle of radius r around a point: p (x,y)
  points = c()
  for(i in 1:360){
    theta = i*2*pi/360
    y = p[2] + r*cos(theta)
    x = p[1] + r*sin(theta)
    points = rbind(points,c(x,y))
  }
  points = rbind(points,points[1,])
  circle = Polygon(points,hole=F)
  circle
}

circs.sp = list()
radii = rep(50,6)

makeMeteCir = function(dataList,radii,laeaProjString){
  # given a list of spatial points dataframes it generates a list of
  # spatial polygons of circles of given a given radius in 
  # lambert equal area projection 
  for(i in seq_along(dataList)){
    atlaslaea = spTransform(meteData[[i]], CRS(laeaProjString))
    if(radii[i] > 0){
      circs = sapply(1:nrow(atlaslaea), 
                function(x){
                  circ = makeCir(atlaslaea@coords[x,],radii[i])
                  circ = Polygons(list(circ),ID=atlaslaea@data$SiteID[x])
                }
              )
      circs.sp[[i]] = SpatialPolygons(circs, proj4string=CRS(laeaProjString))
    }
    else{
     circs.sp[[i]] = SpatialPoints(coordinates(atlaslaea),
                     proj4string=CRS(laeaProjString))
    }
  }
 names(circs.sp) = names(dataList)
 circs.sp
}

setwd(paste(startDir,'/data/',sep=''))
load('meteData.Rdata')

sampleCircles <- makeMeteCir(meteData,radii,laeaPrjString)

save(sampleCircles,file='sampleCircles.Rdata')