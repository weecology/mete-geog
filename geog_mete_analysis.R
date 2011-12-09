## Author: Dan McGlinn
## Contact: daniel.mcglinn@usu.edu
## Date Created: 08/03/2011
## Date Modified: N/A
## Purpose: To link the coordinates of the biotic surveys with environmental data

bio1	Annual Mean Temperature
bio2	Mean Diurnal Range (Mean of monthly (max temp - min temp))
bio3	Isothermality (BIO2/BIO7) (* 100)
bio4	Temperature Seasonality (standard deviation *100)
bio5	Max Temperature of Warmest Month
bio6	Min Temperature of Coldest Month
bio7	Temperature Annual Range (BIO5-BIO6)
bio8	Mean Temperature of Wettest Quarter
bio9	Mean Temperature of Driest Quarter
bio10	Mean Temperature of Warmest Quarter
bio11	Mean Temperature of Coldest Quarter
bio12	Annual Precipitation
bio13	Precipitation of Wettest Month
bio14	Precipitation of Driest Month
bio15	Precipitation Seasonality (Coefficient of Variation)
bio16	Precipitation of Wettest Quarter
bio17	Precipitation of Driest Quarter
bio18	Precipitation of Warmest Quarter
bio19	Precipitation of Coldest Quarter

## PART III - Examine Regression models of S & N 

#read in data for each dataset
setwd('C:/Users/White Lab/Documents/Lab data/MaxEnt geog/data')
file.names<-dir()
file.names<-file.names[grep('_envidat.csv',file.names)]
dat<-list()
for(i in 1:length(file.names)){
 dat[[i]]<-read.table(file.names[i],header=T,sep=',')
}
names(dat)<-unlist(strsplit(file.names,'_dat.csv'))

#################

bbsS.ndvi<-lm(dat$bbs$S~dat$bbs$ndvi)
summary(bbsS.ndvi)

##try to recreate patterns from Hurlbert 2004
par(mfrow=c(3,1))
plot(dat$bbs$ndvi,dat$bbs$S)
plot(dat$bbs$ndvi,dat$bbs$logN)
plot(dat$bbs$logN,dat$bbs$S)

mod.objs<-list()
mod.sums<-array(NA,dim=c(6,4,2),dimnames=list(names(dat),c('r.sq','adj.r.sq','num.vars','num.data'),c('S','N')))
for(i in 1:length(dat)){
 ##drop rows with NA values for now
 dat[[i]]<-dat[[i]][apply(dat[[i]],1,function(x)sum(is.na(x))==0),]
 S.full.mod<-lm(dat[[i]]$S~dat[[i]]$Longitude+dat[[i]]$Latitude+dat[[i]]$bio1+dat[[i]]$bio2+dat[[i]]$bio3+dat[[i]]$bio4+dat[[i]]$bio5+
                         dat[[i]]$bio6+dat[[i]]$bio7+dat[[i]]$bio8+dat[[i]]$bio9+dat[[i]]$bio10+
                         dat[[i]]$bio11+dat[[i]]$bio12+dat[[i]]$bio13+dat[[i]]$bio14+dat[[i]]$bio15+
                         dat[[i]]$bio16+dat[[i]]$bio17+dat[[i]]$bio18+dat[[i]]$ndvi+dat[[i]]$biome)
 N.full.mod<-lm(dat[[i]]$logN~dat[[i]]$Longitude+dat[[i]]$Latitude+dat[[i]]$bio1+dat[[i]]$bio2+dat[[i]]$bio3+dat[[i]]$bio4+dat[[i]]$bio5+
                         dat[[i]]$bio6+dat[[i]]$bio7+dat[[i]]$bio8+dat[[i]]$bio9+dat[[i]]$bio10+
                         dat[[i]]$bio11+dat[[i]]$bio12+dat[[i]]$bio13+dat[[i]]$bio14+dat[[i]]$bio15+
                         dat[[i]]$bio16+dat[[i]]$bio17+dat[[i]]$bio18+dat[[i]]$ndvi+dat[[i]]$biome)

 S.sub.mod<-step(S.full.mod,trace=F)
 N.sub.mod<-step(N.full.mod,trace=F)
 mod.objs[[i]]<-list(S.mod=S.sub.mod,N.mod=N.sub.mod)
 mod.sums[i,,]<-round(cbind(c(summary(S.sub.mod)$r.squared,summary(S.sub.mod)$adj.r.squared,summary(S.sub.mod)$df[1],nrow(dat[[i]])),
                      c(summary(N.sub.mod)$r.squared,summary(N.sub.mod)$adj.r.squared,summary(N.sub.mod)$df[1],nrow(dat[[i]]))),2)
}
names(mod.objs)<-names(dat)

mod.sums
, , S

       r.sq adj.r.sq num.vars num.data
bbs    0.55     0.55       25     2697
cbc    0.73     0.72       28     1874
fia    0.06     0.05       16    10196
gentry 0.09     0.07        7      208
mcdb   0.49     0.36       22      100
naba   0.13     0.11        9      382

, , N

       r.sq adj.r.sq num.vars num.data
bbs    0.27     0.27       25     2697
cbc    0.36     0.35       26     1874
fia    0.16     0.16       23    10196
gentry 0.06     0.03        7      208
mcdb   0.52     0.41       19      100
naba   0.04     0.03        7      382

setwd('C:/Users/White Lab/Documents/Lab data/MaxEnt geog/analysis')
#pdf('geog_prelim_r2.pdf')
plot(mod.sums[,2,1],mod.sums[,2,2],xlab='S adjRsq',ylab='N adjRsq',type='n')
text(mod.sums[,2,1],mod.sums[,2,2],labels=names(dat))
#dev.off()

##geographic variation in residuals
library(maps)
setwd('C:/Users/White Lab/Documents/Lab data/MaxEnt geog/analysis')

S.resid<-residuals(mod.objs$bbs$S.mod)
S.resid.std<-scale(S.resid)
S.resid.pro<-S.resid.std/
#
N.resid<-residuals(mod.objs$bbs$N.mod)
N.resid.std<-scale(N.resid)
pdf('geog_bbs_resid.pdf')
i=1
 par(mfrow=c(2,1))
 map('world',c('usa','canada'),xlim=c(-160,-50),ylim=c(25,55))
 map('state',xlim=c(-160,-50),ylim=c(25,55),add=TRUE)
 mtext(side=3,paste('S residuals: ',names(dat)[i],sep=''))
# points(dat[[i]]$Longitude,dat[[i]]$Latitude,pch='.')
 points(dat[[i]]$Longitude,dat[[i]]$Latitude,cex=.5*S.resid.std,pch=1,col='red')
 points(dat[[i]]$Longitude,dat[[i]]$Latitude,cex=.5*-1*S.resid.std,pch=1,col='blue')
 map('world',c('usa','canada'),xlim=c(-160,-50),ylim=c(25,55))
 map('state',xlim=c(-160,-50),ylim=c(25,55),add=TRUE)
 mtext(side=3,paste('N residuals: ',names(dat)[i],sep=''))
# points(dat[[i]]$Longitude,dat[[i]]$Latitude,pch='.')
 points(dat[[i]]$Longitude,dat[[i]]$Latitude,cex=.5*N.resid.std,pch=1,col='red')
 points(dat[[i]]$Longitude,dat[[i]]$Latitude,cex=.5*-1*N.resid.std,pch=1,col='blue')
dev.off()


pdf('geog_resid_1.pdf')
par(mfrow=c(3,2))
#
for(i in 1:3){
 if(names(dat)[i]=='fia')
  map('state',xlim=c(-100,-70),ylim=c(25,55))
 else{
  map('world',c('usa','canada'),xlim=c(-160,-50),ylim=c(25,55))
  map('state',xlim=c(-160,-50),ylim=c(25,55),add=TRUE)
 }
 mtext(side=3,paste('S residuals: ',names(dat)[i],sep=''))
 points(dat[[i]]$Longitude,dat[[i]]$Latitude,cex=4*residuals(mod.objs[[i]]$S.mod)/max(residuals(mod.objs[[i]]$S.mod)),pch=1,col='blue')
 if(names(dat)[i]=='fia')
  map('state',xlim=c(-100,-70),ylim=c(25,55))
 else{
  map('world',c('usa','canada'),xlim=c(-160,-50),ylim=c(25,55))
  map('state',xlim=c(-160,-50),ylim=c(25,55),add=TRUE)
 }
 mtext(side=3,paste('N residuals: ',names(dat)[i],sep=''))
 points(dat[[i]]$Longitude,dat[[i]]$Latitude,cex=4*residuals(mod.objs[[i]]$N.mod)/max(residuals(mod.objs[[i]]$N.mod)),pch=1,col='red')
}
dev.off()

###
setwd('C:/Users/White Lab/Documents/Lab data/MaxEnt geog/analysis')
pdf('geog_resid_2.pdf')
par(mfrow=c(3,2))
#
for(i in 4:6){
 if(names(dat)[i]=='gentry')
  map('world')
 else{
  map('world',c('usa','canada'),xlim=c(-160,-50),ylim=c(25,55))
  map('state',xlim=c(-160,-50),ylim=c(25,55),add=TRUE)
 }
 mtext(side=3,paste('S residuals: ',names(dat)[i],sep=''))
 points(dat[[i]]$Longitude,dat[[i]]$Latitude,cex=4*residuals(mod.objs[[i]]$S.mod)/max(residuals(mod.objs[[i]]$S.mod)),pch=1,col='blue')
 if(names(dat)[i]=='gentry')
  map('world')
 else{
  map('world',c('usa','canada'),xlim=c(-160,-50),ylim=c(25,55))
  map('state',xlim=c(-160,-50),ylim=c(25,55),add=TRUE)
 }
 mtext(side=3,paste('N residuals: ',names(dat)[i],sep=''))
 points(dat[[i]]$Longitude,dat[[i]]$Latitude,cex=4*residuals(mod.objs[[i]]$N.mod)/max(residuals(mod.objs[[i]]$N.mod)),pch=1,col='red')
}
dev.off()






###################
library(maps)
library(maptools)
library(classInt)
library(raster)
library(rgdal)
library(gstat)
##bbs example
bbs.dat<-cbind(dat$bbs,residuals(mod.objs$bbs$S.mod),residuals(mod.objs$bbs$N.mod))
names(bbs.dat)<-c(names(bbs.dat)[-(33:34)],'Sresid','Nresid')
coordinates(bbs.dat)<-c('Longitude','Latitude')
proj4string(bbs.dat)<-CRS('+proj=longlat')

#Map point data
STATES<-map(database='world',c('usa','canada'),xlim=c(-165,-50),fill=T,plot=F) 
STATES_sp<-map2SpatialPolygons(STATES, IDs=STATES$names, CRS('+proj=longlat'))
states.outline<-list('sp.polygons',STATES_sp)
Nintervals<-classIntervals(c(bbs.dat$Nresid),n=15,style='equal')
Sintervals<-classIntervals(c(bbs.dat$Sresid),n=15,style='equal')
bird.layout<-list(states.outline)

spplot(bbs.dat,'Sresid',sp.layout=bird.layout, cuts=Sintervals$brks,pch=19,
	col.regions=colorRampPalette(c('Blue','Light Blue','Dark Green','Yellow','Red'))(20),
	key.space=list(space='right',border=TRUE))  
windows()
spplot(bbs.dat,'Nresid',sp.layout=bird.layout, cuts=Nintervals$brks,pch=19,
	col.regions=colorRampPalette(c('Blue','Light Blue','Dark Green','Yellow','Red'))(20),
	key.space=list(space='right',border=TRUE))  

usa<-map('usa',fill=T,plot=F)
usa_sp<-map2SpatialPolygons(usa, IDs=usa$names, CRS('+proj=longlat'))

#Make a grid that we will fill 
bb = bbox(STATES_sp) #Defines bounding box
cs = c(0.5,0.5) #Defines cell size in lon/lat
cc = bb[,1] + (cs/2) #Defines cell center offset (location of lower-left corner)
cd = ceiling(diff(t(bb))/cs) #Defines number of cells in each dimension
usa_grid = GridTopology(cc,cs,cd) #(GridTopology(cellcentre.offset,cellsize,cells.dim)
usa_grid
usa_sg = SpatialGrid(usa_grid,proj4string=CRS(proj4string(bbs.dat)))
summary(usa_sg)

#Make variogram
Sresid.vg = variogram(resid~1, bbs.dat, alpha=c(45,135)); plot(Sresid.vg)
Nresid.vg = variogram(resid~1, bbs.dat, alpha=c(45,135)); plot(Nresid.vg)

Sresid.fit = fit.variogram(Sresid.vg,vgm(150,'Sph',700,80))
plot(Sresid.vg,Sresid.fit)

Sresid.ok = krige(resid~1, bbs.dat, usa_sg, Sresid.fit, maxdist=800, nmin=3)

#Remove grid cells whose centers are not inside the polygon of the usa
cells2plot = overlay(Sresid.ok,usa_sp)
cells2plot = ifelse(is.na(cells2plot),NA,1)
Sresid.ok.inner = Sresid.ok
Sresid.ok.inner@data = Sresid.ok.inner@data*cells2plot

ok.data = data.frame(Sresid.pred.mean=Sresid.ok.inner$var1.pred, 
	Sresid.pred.var=Sresid.ok.inner$var1.var)

birds_sgdf = SpatialGridDataFrame(usa_sg,ok.data,proj4string=CRS(proj4string(usa_sg)))
birds_sgdf@data$Sresid.pred.mean = Sresid.ok.inner$var1.pred


###Make a plot
route.locs = list('sp.points',bbs.dat,fill=F,pch=1,col=1) #places points at route locations
scale = list('SpatialPolygonsRescale',layout.scale.bar(),offset=c(-120,30),
	scale=4,fill=c('transparent','black')) #makes scale bar
start = list('sp.text', c(-120,29), '0', cex=0.7)#text for scale bar
stop = list('sp.text', c(-116,29), '4', cex=0.7)
mylayout<-list(route.locs, scale, start, stop)#combine layout attributes into a list

spplot(birds_sgdf,c('Sresid.pred.mean',
	col.regions=colorRampPalette(c('Blue','Light Blue','Dark Green','Yellow','Red'))(21),
	sp.layout = mylayout)



