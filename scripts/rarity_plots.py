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
    plot_color_by_pt_dens(obs_pred_data['pred'], obs_pred_data['obs'], 3, loglog=1)
    plt.loglog([min(obs_pred_data['pred']), max(obs_pred_data['pred'])], 
               [min(obs_pred_data['pred']), max(obs_pred_data['pred'])], 'k-')
    plt.savefig(dest_file, dpi = 400)    

if len(sys.argv) > 1:
    datasets = [sys.argv[1]]
else:
    datasets = ['bbs_2012', 'bbs_2008_2012', 'cbc', 'fia', 'gentry', 'naba']

for dataset in datasets:
    for datatype in ['fit', 'test']:
        for predtype in ['rare']:
            obs_pred_data = read_csv('./results/' + dataset + '_' + predtype+ '_' + datatype + '_obs_pred.csv')
            print obs_pred_rsquare(np.array(obs_pred_data['obs'].values), np.array(obs_pred_data['pred'].values))
            fig_name = './figs/' + dataset + '_' + datatype +'_obs_pred_rarity.png'
            plot_obs_pred(obs_pred_data, dest_file=fig_name)

