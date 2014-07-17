import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from pandas import read_csv
from math import log
import numpy as np
import matplotlib.pylab as plt
import sys

from macroecotools import plot_color_by_pt_dens, obs_pred_rsquare

def plot_obs_pred(obs_pred_data, adj=0, dest_file='./obs_pred.png'):
    plot_color_by_pt_dens(obs_pred_data['pred'] + adj, obs_pred_data['obs'] + adj,
                          3, loglog=1)
    plt.loglog([min(obs_pred_data['pred'] + adj), max(obs_pred_data['pred'] + adj)], 
               [min(obs_pred_data['pred'] + adj), max(obs_pred_data['pred'] + adj)], 'k-')
    plt.savefig(dest_file, dpi = 400)

if len(sys.argv) > 1:
    datasets = [sys.argv[1]]
else:
    datasets = ['bbs_2012', 'bbs_2008_2012', 'cbc', 'gentry', 'naba']


for dataset in datasets:
    for datatype in ['fit', 'test']:
        for predtype in ['rad']:
            obs_pred_data = read_csv('./results/' + dataset + '_' + predtype+ '_' + datatype + '_obs_pred.csv')
            adj = 1
            log_pred = [log(float(i + adj)) for i in obs_pred_data['pred'].values]
            log_obs = [log(float(i + adj)) for i in obs_pred_data['obs'].values]                  
            print obs_pred_rsquare(np.array(log_pred), np.array(log_obs))
            fig_name = './figs/' + dataset + '_' + datatype +'_obs_pred_' + predtype + '.png'
            plot_obs_pred(obs_pred_data, adj=adj, dest_file=fig_name)

