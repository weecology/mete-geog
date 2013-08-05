## Author: Dan McGlinn
## Contact: daniel.mcglinn@usu.edu
## Date Created: 08/03/2011
## Date Modified: N/A
## Purpose: To link the coordinates of the biotic surveys with environmental data
## PART III - Examine Regression models of S & N 

#read in data for each dataset
setwd('~/maxent/geog/')
file.names<-dir('./data/')
file.names<-paste('./data/', file.names[grep('_envidat.csv',file.names)], sep='')
dat<-list()
for(i in 1:length(file.names)){
 dat[[i]]<-read.table(file.names[i],header=T,sep=',')
 dat[[i]]$logS = log10(dat[[i]]$S)
 dat[[i]]$logN = log10(dat[[i]]$N)
}
names(dat)<-unlist(strsplit(file.names,'_envidat.csv'))

bbsS.ndvi<-lm(dat$bbs$S~dat$bbs$ndviJun.mean)
summary(bbsS.ndvi)

##try to recreate patterns from Hurlbert 2004
par(mfrow=c(3,1))
plot(dat$bbs$ndviJun.mean,dat$bbs$S)
plot(dat$bbs$ndviJun.mean,dat$bbs$N)
plot(dat$bbs$N,dat$bbs$S)

mod.objs<-vector('list',6)
mod.sums<-array(NA,dim=c(6,4,4),dimnames=list(names(dat),
                c('r.sq','adj.r.sq','num.vars','num.data'),c('S','N','logS','logN')))
for(i in 1:length(dat)){
 ##drop rows with NA values for now
 dat[[i]]<-dat[[i]][apply(dat[[i]],1,function(x)sum(is.na(x))==0),]
 S.full.mod<-lm(S~Longitude+Latitude+alt.mean+ndviJun.mean+ndviDec.mean+   
                mdr.mean+iso.mean+tseas.mean+tmax.mean+tmin.mean+tar.mean+
                twetq.mean+tdryq.mean+twarmq.mean+tcoldq.mean+ap.mean+pwet.mean+
                pdry.mean+pseas.mean+pwetq.mean+pdryq.mean+pwarmq.mean+pcoldq.mean+
                mat.mean+alt.var+ndviJun.var+ndviDec.var+mdr.var+iso.var+
                tseas.var+tmax.var+tmin.var+tar.var+twetq.var+tdryq.var+twarmq.var+
                tcoldq.var+ap.var+pwet.var+pdry.var+pseas.var+pwetq.var+pdryq.var+
                pwarmq.var+pcoldq.var+mat.var+biome,data=dat[[i]])
 N.full.mod<-lm(N~Longitude+Latitude+alt.mean+ndviJun.mean+ndviDec.mean+   
                mdr.mean+iso.mean+tseas.mean+tmax.mean+tmin.mean+tar.mean+
                twetq.mean+tdryq.mean+twarmq.mean+tcoldq.mean+ap.mean+pwet.mean+
                pdry.mean+pseas.mean+pwetq.mean+pdryq.mean+pwarmq.mean+pcoldq.mean+
                mat.mean+alt.var+ndviJun.var+ndviDec.var+mdr.var+iso.var+
                tseas.var+tmax.var+tmin.var+tar.var+twetq.var+tdryq.var+twarmq.var+
                tcoldq.var+ap.var+pwet.var+pdry.var+pseas.var+pwetq.var+pdryq.var+
                pwarmq.var+pcoldq.var+mat.var+biome,data=dat[[i]])
 logS.full.mod<-lm(logS~Longitude+Latitude+alt.mean+ndviJun.mean+ndviDec.mean+   
                mdr.mean+iso.mean+tseas.mean+tmax.mean+tmin.mean+tar.mean+
                twetq.mean+tdryq.mean+twarmq.mean+tcoldq.mean+ap.mean+pwet.mean+
                pdry.mean+pseas.mean+pwetq.mean+pdryq.mean+pwarmq.mean+pcoldq.mean+
                mat.mean+alt.var+ndviJun.var+ndviDec.var+mdr.var+iso.var+
                tseas.var+tmax.var+tmin.var+tar.var+twetq.var+tdryq.var+twarmq.var+
                tcoldq.var+ap.var+pwet.var+pdry.var+pseas.var+pwetq.var+pdryq.var+
                pwarmq.var+pcoldq.var+mat.var+biome,data=dat[[i]])
 logN.full.mod<-lm(logN~Longitude+Latitude+alt.mean+ndviJun.mean+ndviDec.mean+   
                mdr.mean+iso.mean+tseas.mean+tmax.mean+tmin.mean+tar.mean+
                twetq.mean+tdryq.mean+twarmq.mean+tcoldq.mean+ap.mean+pwet.mean+
                pdry.mean+pseas.mean+pwetq.mean+pdryq.mean+pwarmq.mean+pcoldq.mean+
                mat.mean+alt.var+ndviJun.var+ndviDec.var+mdr.var+iso.var+
                tseas.var+tmax.var+tmin.var+tar.var+twetq.var+tdryq.var+twarmq.var+
                tcoldq.var+ap.var+pwet.var+pdry.var+pseas.var+pwetq.var+pdryq.var+
                pwarmq.var+pcoldq.var+mat.var+biome,data=dat[[i]])
 S.sub.mod<-step(S.full.mod,trace=F)
 N.sub.mod<-step(N.full.mod,trace=F)
 logS.sub.mod<-step(logS.full.mod,trace=F)
 logN.sub.mod<-step(logN.full.mod,trace=F)
 mod.objs[[i]]<-list(S.mod=S.sub.mod,N.mod=N.sub.mod,logS.mod=logS.sub.mod,
                     logN.mod=logN.sub.mod)
 mod.sums[i,,]<-round(cbind(c(summary(S.sub.mod)$r.squared,
                              summary(S.sub.mod)$adj.r.squared,
                              summary(S.sub.mod)$df[1],nrow(dat[[i]])),
                            c(summary(N.sub.mod)$r.squared,
                              summary(N.sub.mod)$adj.r.squared,
                              summary(N.sub.mod)$df[1],nrow(dat[[i]])),
                            c(summary(logS.sub.mod)$r.squared,
                              summary(logS.sub.mod)$adj.r.squared,
                              summary(logS.sub.mod)$df[1],nrow(dat[[i]])),
                            c(summary(logN.sub.mod)$r.squared,
                              summary(logN.sub.mod)$adj.r.squared,
                              summary(logN.sub.mod)$df[1],nrow(dat[[i]]))),2)
}
names(mod.objs)<-names(dat)

save(dat,mod.objs,mod.sums,file='geog_stepwise_input_&_output.Rdata')
#load('geog_stepwise_input_&_output.Rdata')
mod.sums
, , S

       r.sq adj.r.sq num.vars num.data
bbs    0.56     0.55       35     2672
cbc    0.73     0.72       44     1893
fia    0.06     0.06       26    10350
gentry 0.81     0.22       26       34
mcdb   0.87     0.62       45       69
nabc   0.24     0.20       20      388

, , N

       r.sq adj.r.sq num.vars num.data
bbs    0.24     0.23       38     2672
cbc    0.02     0.01       12     1893
fia    0.19     0.19       36    10350
gentry 0.86     0.42       26       34
mcdb   0.98     0.95       45       69
nabc   0.05     0.04        7      388

, , logS

       r.sq adj.r.sq num.vars num.data
bbs    0.58     0.57       38     2672
cbc    0.71     0.70       43     1893
fia    0.06     0.06       27    10350
gentry 0.86     0.44       26       34
mcdb   0.86     0.62       45       69
nabc   0.20     0.16       18      388

, , logN

       r.sq adj.r.sq num.vars num.data
bbs    0.31     0.30       42     2672
cbc    0.33     0.32       33     1893
fia    0.18     0.18       36    10350
gentry 0.89     0.53       26       34
mcdb   0.95     0.87       44       69
nabc   0.11     0.07       17      388


pdf('./figs/prelim_r2adj.pdf')
plot(mod.sums[,2,1],mod.sums[,2,2],xlab='S adjRsq',ylab='N adjRsq',type='n')
text(mod.sums[,2,1],mod.sums[,2,2],labels=names(dat))
dev.off()

pdf('./figs/prelim_r2_both.pdf')
plot(mod.sums[,2,1],mod.sums[,2,2],xlab='S Rsq',ylab='N Rsq',type='n',ylim=c(0,1),
     xlim=c(0,1))
text(mod.sums[,1,1],mod.sums[,1,2],labels=names(dat),col='red')
text(mod.sums[,2,1],mod.sums[,2,2],labels=names(dat))
legend('topleft',c('R sqr','adjusted R sqr'),col=2:1,lty=1,bty='n')
abline(a=0,b=1,lty=2)
dev.off()

pdf('./figs/prelim_r2_both_loglog.pdf')
plot(mod.sums[,2,3],mod.sums[,2,4],xlab='logS Rsq',ylab='logN Rsq',type='n',ylim=c(0,1),
     xlim=c(0,1))
text(mod.sums[,1,3],mod.sums[,1,4],labels=names(dat),col='red')
text(mod.sums[,2,3],mod.sums[,2,4],labels=names(dat))
legend('topleft',c('R sqr','adjusted R sqr'),col=2:1,lty=1,bty='n')
abline(a=0,b=1,lty=2)
dev.off()

##geographic variation in residuals
library(maps)

S.resid<-residuals(mod.objs$bbs$S.mod)
S.resid.std<-scale(S.resid)
S.resid.pro<-S.resid.std/
#
N.resid<-residuals(mod.objs$bbs$N.mod)
N.resid.std<-scale(N.resid)
pdf('./figs/bbs_resid.pdf')
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


pdf('./figs/geog_resid_1.pdf')
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
pdf('./figs/geog_resid_2.pdf')
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



