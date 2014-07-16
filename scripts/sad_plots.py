import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from pandas import read_csv
from math import log
import numpy as np
import matplotlib.pylab as plt
import sys

from macroecotools import plot_color_by_pt_dens, obs_pred_rsquare

def import_data(resultdir, rawdir, datatype, dataset):
    envpred_data = read_csv(resultdir + dataset + '_state_var_' + datatype + '_obs_pred.csv')
    sad_data = read_csv(rawdir + dataset +'_spab.csv')
    # in some cases there are different sites in *_state_var_obs_pred.csv then in *_spab.csv
    # Deal with the fact that ied.csv
    sad_data = sad_data.ix[np.in1d(sad_data['site_id'].values, envpred_data['site_id'].values)]
    return envpred_data, sad_data

def plot_obs_pred(sad_data, envpred_sads, dest_file='./obs_pred.png'):
    plot_color_by_pt_dens(envpred_sads['EnvPred'], sad_data['Obs'], 3, loglog=1)
    plt.loglog([min(envpred_sads['EnvPred']), max(envpred_sads['EnvPred'])], 
               [min(envpred_sads['EnvPred']), max(envpred_sads['EnvPred'])], 'k-')
    plt.savefig(dest_file, dpi = 400)

if len(sys.argv) > 1:
    datasets = [sys.argv[1]]
else:
    datasets = ['bbs_2012', 'bbs_2008_2012', 'cbc', 'fia', 'gentry', 'nabc']


for dataset in datasets:
    for datatype in ['fit', 'test']:
        envpred_sads = read_csv('./data/%s_envpred_sads.csv' % dataset)
        log_pred = [log(float(i)) for i in envpred_sads['EnvPred'].values]
        log_obs = [log(float(i)) for i in sad_data['Obs'].values]    
        print obs_pred_rsquare(np.array(log_pred), np.array(log_obs))
        plot_obs_pred(sad_data, envpred_sads, dest_file='./figs/%s_obs_pred.png' % dataset)        
    