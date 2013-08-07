"""Modeling geographic variation in macroecological patterns

Project code for attempts to model geographic variation in macroecological
patterns by using evironmental variables to model species richness and total
abundance and then using MaxEnt to model the macroecological patterns based
on those constraints.

At the moment all of the environmental data processing and the modeling
linking richness and abundance to environmental variables is being conducted
in R.

Currently this is just a quick test run focusing on how well this idea can
work using the Breeding Bird Survey data alone.

"""
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from pandas import DataFrame, read_csv
from math import log
import numpy as np
import matplotlib.pylab as plt

from mete import get_beta, get_mete_rad
from macroecotools import plot_color_by_pt_dens, obs_pred_rsquare
#import working_functions *


def import_data(datadir, dataset):
    envpred_data = read_csv(datadir + dataset + '_out.csv')
    sad_data = read_csv(datadir + dataset +'_obs_pred_sub.csv', names=['SiteID',
                                                                  'Year',
                                                                  'Obs',
                                                                  'Pred'])
    #Deal with the fact that in some cases there are different sites in
    #*_out.csv and *_obs_pred.csv
    #sad_data = sad_data.ix[np.in1d(sad_data['SiteID'].values, envpred_data['SiteID'].values)]
    return envpred_data, sad_data

def get_envpred_sads(envpred_data):
    envpred_sads = DataFrame(columns=['SiteID', 'EnvPred'])
    for index, site in envpred_data.iterrows():
        obs_S = site['S']
        envpred_S = 10 ** site['predlogS']
        envpred_N = 10 ** site['predlogN']
        beta = get_beta(envpred_S, envpred_N)
        #To produce a comparable number of species use obs_S; IS THIS RIGHT?
        site_sad, p = get_mete_rad(obs_S, envpred_N, beta=beta)
        site_ids = [site['SiteID'] for i in range(0, len(site_sad))]
        site_sad_with_id = DataFrame(np.column_stack([site_ids, site_sad]),
                                     columns=['SiteID', 'EnvPred'])
        envpred_sads = envpred_sads.append(site_sad_with_id, ignore_index=True)
    return envpred_sads

def plot_obs_pred(sad_data, envpred_sads, dest_file='./obs_pred.png'):
    plot_color_by_pt_dens(envpred_sads['EnvPred'], sad_data['Obs'], 3, loglog=1,
                          plot_obj=ax)
    plt.loglog([min(envpred_sads['EnvPred']), max(envpred_sads['EnvPred'])], 
               [min(envpred_sads['EnvPred']), max(envpred_sads['EnvPred'])], 'k-')
    plt.savefig(dest_file, dpi = 400)
    
datasets = ['bbs', 'cbc', 'fia', 'nabc']

for dataset in datasets:
    envpred_data, sad_data = import_data('./data/', dataset)
    envpred_sads = get_envpred_sads(envpred_data)
    envpred_sads.to_csv('./data/%s_envpred_sads.csv' % dataset)
    #envpred_sads = read_csv('./data/%s_envpred_sads.csv' % dataset)
    log_pred = [log(float(i)) for i in envpred_sads['EnvPred'].values]
    log_obs = [log(float(i)) for i in sad_data['Obs'].values]    
    print obs_pred_rsquare(np.array(log_pred), np.array(log_obs))
    plot_obs_pred(sad_data, envpred_sads)
    
    