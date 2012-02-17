setwd('/home/danmcglinn/maxent/geog')
load('geog_stepwise_input_&_output.Rdata')

##drop predictors at random from gentry and mcdb
set.seed(1)
for(i in 4:5){
  for(j in 1:4){
    mod.split = strsplit(as.character(formula(mod.objs[[i]][[j]])),' ')
    theCoeffs = mod.split[[3]][mod.split[[3]]!='+']
    ncoef = floor(length(theCoeffs) / 2)
    mod.objs[[i]][[j]] = lm(paste(mod.split[[2]],'~',
                         paste(sample(theCoeffs,ncoef),sep='',collapse='+'),sep=''),
                         data=dat[[i]])
    mod.sum = summary(mod.objs[[i]][[j]])
    mod.sums[i,1,j] = round(mod.sum$r.squared,2)
    mod.sums[i,2,j] = round(mod.sum$adj.r.squared,2)
    mod.sums[i,3,j] = round(ncoef,2)                   
  }
} 
  
## output predictions to files
for(i in seq_along(dat)){
  mod = mod.objs[[i]]
  out = dat[[i]][,1:6]
  out = data.frame(out,logS=log10(out$S),logN=log10(out$N))
  out = data.frame(out,dat[[i]][,7:10])
  out = data.frame(out,predS=predict(mod$S),predN=predict(mod$N),
                   predlogS=predict(mod$logS),predlogN=predict(mod$logN))
  write.csv(out,paste('./data/',names(dat)[i],'_out.csv',sep=''),row.names=F)
}

mod.sums
, , S

       r.sq adj.r.sq num.vars num.data
bbs    0.56     0.55       35     2672
cbc    0.73     0.72       44     1893
fia    0.06     0.06       26    10350
gentry 0.21    -0.24       12       34
mcdb   0.37     0.10       21       69
naba   0.24     0.20       20      388

, , N

       r.sq adj.r.sq num.vars num.data
bbs    0.24     0.23       38     2672
cbc    0.02     0.01       12     1893
fia    0.19     0.19       36    10350
gentry 0.42     0.10       12       34
mcdb   0.48     0.25       21       69
naba   0.05     0.04        7      388

, , logS

       r.sq adj.r.sq num.vars num.data
bbs    0.58     0.57       38     2672
cbc    0.71     0.70       43     1893
fia    0.06     0.06       27    10350
gentry 0.22    -0.22       12       34
mcdb   0.45     0.20       21       69
naba   0.20     0.16       18      388

, , logN

       r.sq adj.r.sq num.vars num.data
bbs    0.31     0.30       42     2672
cbc    0.33     0.32       33     1893
fia    0.18     0.18       36    10350
gentry 0.36     0.00       12       34
mcdb   0.53     0.31       21       69
naba   0.11     0.07       17      388



  