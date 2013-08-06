setwd('~/maxent/geog')
load('./results/stepwise_input_&_output.Rdata')

sub_names = c('bbs', 'cbc', 'nabc', 'fia')
true = names(mod.objs) %in% sub_names
mod.sums = mod.sums[true, , ]
mod.objs = mod.objs[true]
dat = dat[true]

##drop predictors at random from gentry and mcdb
#set.seed(1)
#for(i in 4:5){
#  for(j in 1:4){
#    mod.split = strsplit(as.character(formula(mod.objs[[i]][[j]])),' ')
#    theCoeffs = mod.split[[3]][mod.split[[3]]!='+']
#    ncoef = floor(length(theCoeffs) / 2)
#    mod.objs[[i]][[j]] = lm(paste(mod.split[[2]],'~',
#                         paste(sample(theCoeffs,ncoef),sep='',collapse='+'),sep=''),
#                         data=dat[[i]])
#    mod.sum = summary(mod.objs[[i]][[j]])
#    mod.sums[i,1,j] = round(mod.sum$r.squared,2)
#    mod.sums[i,2,j] = round(mod.sum$adj.r.squared,2)
#    mod.sums[i,3,j] = round(ncoef,2)                   
#  }
#} 
  
## output predictions to files
fields = c("Longitude","Latitude","SiteID","year","S","N")

for(i in seq_along(dat)){
  mod = mod.objs[[i]]
  out = dat[[i]][,fields]
  out = data.frame(out,logS=log10(out$S),logN=log10(out$N))
  out = data.frame(out,dat[[i]][ ,c('p1', 'p2', 'p3', 'p4')])
  out = data.frame(out,predS=predict(mod$S),predN=predict(mod$N),
                   predlogS=predict(mod$logS),predlogN=predict(mod$logN))
  write.csv(out,paste('./data/',names(dat)[i],'_out.csv',sep=''),row.names=F)
}

  