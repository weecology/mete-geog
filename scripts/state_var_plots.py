import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from pandas import read_csv
from math import log
import numpy as np
import matplotlib.pylab as plt
import sys

from macroecotools import plot_color_by_pt_dens, obs_pred_rsquare

def plot_obs_pred(obs_pred_data, yvar, dest_file='./obs_pred.png'):
    plot_color_by_pt_dens(10**obs_pred_data[yvar + 'pred'], 10**obs_pred_data[yvar],
                          3, loglog=1)
    plt.loglog([min(10**obs_pred_data['pred'] + adj), max(10**obs_pred_data['pred'] + adj)], 
               [min(10**obs_pred_data['pred'] + adj), max(10**obs_pred_data['pred'] + adj)], 'k-')
    plt.savefig(dest_file, dpi = 400)

if len(sys.argv) > 1:
    datasets = [sys.argv[1]]
else:
    datasets = ['bbs_2012', 'bbs_2008_2012', 'cbc', 'fia', 'gentry', 'naba']

yvars = ['lsr', 'lab']

for dataset in datasets:
    for datatype in ['fit', 'test']:
        obs_pred_data = read_csv('./results/' + dataset + '_state_var_' + datatype + '_obs_pred.csv')
        for yvar in yvars:
            fig_name = './figs/' + dataset + '_sate_var_' + datatype + '_obs_pred.png'
            plot_obs_pred(obs_pred_data, adj=adj, dest_file=fig_name)


