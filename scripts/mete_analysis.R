## Author: Dan McGlinn
## Contact: daniel.mcglinn@usu.edu
## Date Created: 08/03/2011
## Date Modified: N/A
## Purpose: To link the coordinates of the biotic surveys with environmental data
## PART III - Examine Regression models of S & N 

library(glmnet)
library(dummies)

library(foreach)
library(snowfall)
library(doSNOW)

npar = 6
sfInit(parallel=TRUE, cpus=npar, type="SOCK")
sfLibrary(dummies)
sfLibrary(glmnet)
registerDoSNOW(sfGetCluster())

## set smallest acceptable value of richness
minS = 5

#read in data for each dataset
radii = 50
for (r in radii) {
    load(paste('./data/mete_data_biotic_and_envi_', r, '.Rdata', sep=''))
    sfExport('mete_data')
    foreach (i= 1:6, .inorder=FALSE) %dopar% { 
        datname = names(mete_data)[i]
        dat = mete_data[[i]]
        
        if (datname == 'gentry') {
            hab_dum = dummy(dat$biome)
        }
        else {
            hab_dum = dummy(dat$landuse)
            #hab_dum = hab_dum[ , !(colnames(hab_dum) == 'landuseNA')]
        }
        dat = cbind(dat, hab_dum)
        ## if sampling effort is zero change to NA
        ## then log transform 
        if (datname == 'cbc') {
            dat$duration_hrs[dat$duration_hrs == 0] = NA
            dat$duration_hrs = log10(dat$duration_hrs)
        }
        if (datname == 'naba') {
            dat$party_hours = log10(dat$party_hours)
        }
        
        lo_s_rows = dat$sr < minS      
        dat = dat[!lo_s_rows, ]
        
        cols_to_drop = c('site_id', 'sr', 'ab', 'logS', 'logN', 
                         'cn', 'statecd', 'unitcd', 'countrycd',
                         'plot', 'year', 'elev', 'biome', 'landuse')
        cols_to_drop = c(cols_to_drop, paste(cols_to_drop, 1, sep='.'))
        
        col_indices = match(cols_to_drop, names(dat))
        col_indices = col_indices[!is.na(col_indices)]
        x = dat[ , -col_indices]
        
        na_rows = apply(x, 1, function(x) any(is.na(x)))
        x = x[!na_rows, ]
        
        # pull out 20% for cross-validation
        n_test_pts = round(nrow(x) * .2)
        test_rows = sample(nrow(x), size=n_test_pts)
        x_test = as.matrix(x[test_rows, ])
        x = as.matrix(x[-test_rows, ])
        
        # track site_ids
        site_ids = dat$site_id[!na_rows]
        site_ids_test = site_ids[test_rows]
        site_ids_fit = site_ids[-test_rows]
        
        yvars = c('sr', 'lsr', 'ab', 'lab', 'ns', 'lns')
        
        mod_objs = pred_out = vector('list', length(yvars))
        names(mod_objs) = names(pred_out) = yvars
        
        for (j in seq_along(yvars)) {
            if (yvars[j] %in% c('sr', 'ab')) {
                y = dat[ , yvars[j]]
                y = y[!na_rows]
                y_test = as.matrix(y[test_rows])
                y = as.matrix(y[-test_rows])
                mod_objs[[j]] = cv.glmnet(x, y, family='poisson', alpha=0)
            }
            else {
                if (yvars[j] == 'lsr')
                    y = log10(dat$sr)
                if (yvars[j] == 'lab')
                    y = log10(dat$'ab')
                if (yvars[j] == 'ns')
                    y = dat$ab / dat$sr
                if (yvars[j] == 'lns')
                    y = log10(dat$ab / dat$sr)
                y = y[!na_rows]
                y_test = as.matrix(y[test_rows])
                y = as.matrix(y[-test_rows])
                mod_objs[[j]] = cv.glmnet(x, y, family='gaussian', alpha=0)
            }
            pred_fit = predict(mod_objs[[j]], newx=x, s='lambda.min',
                               type='response')
            pred_test = predict(mod_objs[[j]], newx=x_test, s='lambda.min',
                                type='response')
            pred_out[[j]]$fit = data.frame(site_id = site_ids_fit, 
                                           obs = y, pred=pred_fit[,1])
            pred_out[[j]]$test = data.frame(site_id = site_ids_test,
                                            obs = y_test, pred=pred_test[,1])
            print(paste('completed', yvars[j], 'for', datname, 'at', r))
        }
        ## save model results
        save(mod_objs, pred_out, 
             file=paste('./results/', datname, '_mod_pred_', r, '.Rdata',
                        sep=''))
    }
}

sfStop()


