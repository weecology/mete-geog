library(sp)
library(raster)
library(rgdal)

## PART II - query climate and ecosystem data

extent_x = c(-167,-167,-40,-40,-167)
extent_y = c(-30,70,70,-30,-30)

## set study exent
box = c(-167, -40, -30, 70) 
my_extent = extent(box)
my_res = 0.1

## several different equal area projections to consider
## ref: http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
wgs_prj = "+proj=longlat +ellps=WGS84 "
# lambert amizmial equal area
laea_prj = "+proj=laea +lat_0=20 +lon_0=-103.5 +units=km"
# sinusodial
sinu_prj = "+proj=sinu +lat_0=20 +lon_0=-103.5 +units=km"
# cyindrical
cea_prj = "+proj=cea"


ras_template = raster(my_extent, crs=wgs_prj)
res(ras_template) = my_res

## Read in raster datasets (Climate and NDVI) and crop to myExtent
## load bioclim datasets
## 5 arc-minutes resolution
load('~/gis/WorldClimData/bioclim_5m.Rdata') 
bio_stack = crop(bioStack, my_extent)
## load worldclim elev data
## 5 arc-minutes resolution
load('~/gis/WorldClimData/alt_5m.Rdata')
alt = crop(alt, my_extent)

## load NDVI data
ndvi_jun = stack('~/gis/GlobalVegData/Mean&Std/JUNAV18.rst') 
ndvi_dec = stack('~/gis/GlobalVegData/Mean&Std/decav18.rst') 
ndvi = stack(ndvi_jun, ndvi_dec)
ndvi = crop(ndvi, my_extent)
fix_ndvi = function(x) {
    x2 = ifelse(x == 1, NA, x)
    (x2 - 128) / 128
}
ndvi = calc(ndvi, fix_ndvi)

## resample WorldClim layers to resolution of NDVI before stacking them
bio_temp = resample(bio_stack, ras_template, method='bilinear',
                   progress='text')
bio_stack = bio_temp
alt_temp = resample(alt, ras_template, method='bilinear',
                   progress='text')
alt = alt_temp
rm(alt_temp, bio_temp)
gc()

## stack all the layers
envi_layers = stack(alt, ndvi, bio_stack)
names(envi_layers) = c('alt','ndvi_june','ndvi_dec',
                        names(bio_stack))

## reproject to sinusodial equal area
temp = envi_layers[[1]]
temp.ext = projectExtent(temp, CRS(sinu_prj))
envi_layers_sinu = projectRaster(from=envi_layers, to=temp.ext,
                                crs=CRS(sinu_prj),
                                progress='text')

writeRaster(envi_layers_sinu, file='./data/envi_layers.grd',
            bandorder='BIL')
#envi_Layers = stack('./data/envi_layers.grd')

rm(temp, temp.ext, alt, bioStack, bio_stack, envi_layers,
   ndvi_dec, ndvi_jun, wwfeco)
gc()

## Draw sample circles around atlas locations
make_cir = function(p, r) {
    # A function that draws a circle of radius r around a point: p (x,y)
    points = c()
    for (i in 1:360) {
        theta = i * 2 * pi / 360
        y = p[2] + r * cos(theta)
        x = p[1] + r * sin(theta)
        points = rbind(points, c(x, y))
    }
    points = rbind(points, points[1, ])
    circle = Polygon(points, hole=F)
    return(circle)
}

make_mete_cir = function(mete_data, radius, proj_string){
    # given a list of spatial points dataframes it generates a list of
    # spatial polygons of circles of given a radius in 
    # a supplied projection
    circs_sp = vector('list', length(mete_data))
    for (i in seq_along(mete_data)) {
        atlas = spTransform(mete_data[[i]], CRS(proj_string))
        if (radius > 0) {
            circs = sapply(1:nrow(atlas), function(x) {
                           circ = make_cir(atlas@coords[x, ], radius)
                           circ = Polygons(list(circ), ID=atlas$site_id[x])})   
            circs_sp[[i]] = SpatialPolygons(circs, proj4string=
                                            CRS(proj_string))
        }
        else {
            circs_sp[[i]] = SpatialPoints(coordinates(atlas),
                                          proj4string=CRS(proj_string))
        }
    }
    names(circs_sp) = names(mete_data)
    return(circs_sp)
}

load('./data/mete_data.Rdata')

radii = c(50, 25, 10)

library(snowfall)
library(foreach)

npar = 3
sfInit( parallel=TRUE, cpus=npar, type="SOCK")
sfLibrary(sp)
sfLibrary(rgdal)
sfLibrary(raster)
sfExport("mete_data", "radii", "sinu_prj", "make_cir", "make_mete_cir")
foreach (i=1:length(radii)) %dopar% {
    samp_circles = make_mete_cir(mete_data, radii[i], sinu_prj)
    save(samp_circles, 
         file=paste('./data/samp_circles_', radii[i],'.Rdata', sep='')) 
}
sfStop()
