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

from pandas import DataFrame, read_csv
import numpy as np
import matplotlib.pylab as plt

from mete import get_beta, get_mete_rad
from macroeco import plot_color_by_pt_dens

def import_data(datadir, dataset):
    envpred_data = read_csv(datadir + dataset + '_out.csv')
    sad_data = read_csv(datadir + dataset +'_obs_pred.csv', names=['SiteID',
                                                                  'Year',
                                                                  'Obs',
                                                                  'Pred'])
    #Deal with the fact that in some cases there are different sites in
    #*_out.csv and *_obs_pred.csv
    sad_data = sad_data.ix[np.in1d(sad_data['SiteID'].values,
                                   envpred_data['SiteID'].values)]
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
        site_sad_with_id = DataFrame(np.column_stack([site['SiteID'] *
                                                      np.ones(len(site_sad)),
                                                      site_sad]),
                                     columns=['SiteID', 'EnvPred'])
        envpred_sads = envpred_sads.append(site_sad_with_id, ignore_index=True)
    return envpred_sads

def plot_obs_pred(sad_data, envpred_sads):
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.loglog(sad_data['Pred'], sad_data['Obs'], 'bo')
    plt.loglog([min(sad_data['Pred']), max(sad_data['Pred'])], 
               [min(sad_data['Pred']), max(sad_data['Pred'])], 'k-')
    plt.subplot(1, 2, 2)
    plt.loglog(envpred_sads['EnvPred'], sad_data['Obs'], 'ro')
    plt.loglog([min(envpred_sads['EnvPred']), max(envpred_sads['EnvPred'])], 
               [min(envpred_sads['EnvPred']), max(envpred_sads['EnvPred'])], 'k-')

datasets = ['bbs', 'cbc', 'fia', 'gentry', 'mcdb', 'naba']
for dataset in datasets:
    envpred_data, sad_data = import_data('./data/', dataset)
    envpred_sads = get_envpred_sads(envpred_data)
    plot_obs_pred(sad_data, envpred_sads)