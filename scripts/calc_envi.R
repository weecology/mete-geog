library(raster)
library(rgdal)
library(foreach)
library(doSNOW)
library(snowfall)

radii = c(50, 25, 15)

load('./data/mete_data.Rdata')
sinu_prj = "+proj=sinu +lat_0=20 +lon_0=-103.5 +units=km"
mete_data = lapply(mete_data, function(x) spTransform(x, CRS(sinu_prj)))
envi_layers = stack('./data/envi_layers_sinu.grd')

npar = 8
sfInit(parallel=TRUE, cpus=npar, type="SOCK")
sfLibrary(sp)
sfLibrary(raster)
registerDoSNOW(sfGetCluster())
sfExport("envi_layers", "mete_data")

foreach (r=radii, .inorder=FALSE) %dopar% {
    foreach (i=1:length(mete_data), .inorder=FALSE) %dopar% {
        ## Sample mean and variance inside circles
        envi_avg = extract(envi_layers, coordinates(mete_data[[i]]),
                           fun=mean, na.rm=T, buffer=r,
                           layer=1, nl=nlayers(envi_layers), progress='text')
        
        colnames(envi_avg) = sapply(colnames(envi_avg), function(x) {
                                    paste(x, '.mean', sep='')})
        envi_var = extract(envi_layers, coordinates(mete_data[[i]]),
                           fun=var, na.rm=T, buffer=r,
                           layer=1, nl=nlayers(envi_layers), progress='text')
        
        colnames(envi_var) = sapply(colnames(envi_var), function(x) {
                                    paste(x, '.var', sep='')})
        ## Output data product
        name_dat = names(mete_data)[i]
        out_dat = data.frame(mete_data[[i]], envi_avg, envi_var)
        cols_to_drop = c('longitude.1', 'latitude.1')
        out_dat = out_dat[ , !(names(out_dat) %in% cols_to_drop)]
        filename = paste('./data/', name_dat, '_envi_mean_var_', r, '.csv',
                         sep='')
        write.csv(out_dat, file=filename, row.names=F)
        print(paste('r is', r, 'i is', i))
    }
}

sfStop()
