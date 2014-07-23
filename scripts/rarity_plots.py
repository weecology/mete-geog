import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from pandas import read_csv
from math import log
import numpy as np
import matplotlib.pylab as plt
import sys

from macroecotools import plot_color_by_pt_dens, obs_pred_rsquare

def plot_obs_pred(obs_pred_data, dest_file='./obs_pred.png'):
    sites = list(set(obs_pred_data['site_id'].values))
    obs_tot = []
    pred_tot = []
    for i in sites:
        tmp = obs_pred_data[obs_pred_data['site_id'] == i]
        obs_tot.append(sum(tmp['obs'].values))
        pred_tot.append(sum(tmp['pred'].values))    
    plot_color_by_pt_dens(pred_tot, obs_tot, 3, loglog=1)
    plt.loglog([min(pred_tot), max(pred_tot)], 
               [min(pred_tot), max(pred_tot)], 'k-')
    plt.savefig(dest_file, dpi = 400)    

if len(sys.argv) > 1:
    datasets = [sys.argv[1]]
else:
    datasets = ['bbs_2012', 'bbs_2008_2012', 'cbc', 'fia', 'gentry', 'naba']

for dataset in datasets:
    for datatype in ['fit', 'test']:
        for predtype in ['sad']:
            obs_pred_data = read_csv('./results/' + dataset + '_' + predtype+ '_' + datatype + '_obs_pred_norm.csv')
            obs_pred_data = obs_pred_data[obs_pred_data['octave'] <= 4]
            sites = list(set(obs_pred_data['site_id'].values))
            obs_tot = []
            pred_tot = []
            for i in sites:
                tmp = obs_pred_data[obs_pred_data['site_id'] == i]
                obs_tot.append(sum(tmp['obs'].values))
                pred_tot.append(sum(tmp['pred'].values))
            print obs_pred_rsquare(np.array(pred_tot), np.array(obs_tot))
            fig_name = './figs/' + dataset + '_' + datatype +'_obs_pred_rarity.png'
            plot_obs_pred(obs_pred_data, dest_file=fig_name)

