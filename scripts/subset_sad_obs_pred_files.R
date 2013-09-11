
setwd('~/maxent/geog')

datasets = c('bbs', 'cbc', 'nabc', 'fia')

for(i in seq_along(datasets)){
  out = read.csv(paste('./data/', datasets[i], '_out.csv', sep=''))
  obs_pred = read.csv(paste('./data/', datasets[i], '_obs_pred.csv', sep=''),
                      header=F)
  names(obs_pred) = c('SiteID', 'Year', 'Obs', 'Pred')
  uni_site_ids = unique(out$SiteID)
  obs_pred = obs_pred[obs_pred$SiteID %in% uni_site_ids, ]
  write.table(obs_pred, file=paste('./data/', datasets[i], '_obs_pred_sub.csv', sep=''),
              row.names=F, col.names=F, sep=',')
}

