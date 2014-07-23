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
from math import log, ceil
import numpy as np
import sys

from mete import get_mete_rad, get_mete_sad, which
from macroecotools import plot_color_by_pt_dens, obs_pred_rsquare

def get_log_bins(abu, logbase=2):
    max_bin = log(max(abu), logbase)
    if int(max_bin) == max_bin: # if is an integer then we need one additional bin
        max_bin = int(max_bin + 1)
    else:
        max_bin = int(ceil(max_bin))
    bins = [logbase ** x for x in range(0, max_bin + 1)]
    return bins

def get_envpred(envpred_data, predtype=['sad', 'rad']):
    if predtype is 'sad':
        envpred = DataFrame(columns=['site_id', 'octave', 'env_pred'])
    if predtype is 'rad':
        envpred = DataFrame(columns=['site_id', 'rank', 'env_pred'])
    if pretype is 'rare':
        envpred = DataFrame(columns=['site_id', 'env_pred'])
    for index, site in envpred_data.iterrows():
        obs_S = site['S']
        envpred_S = 10 ** site['logSpred']
        envpred_N = 10 ** site['logNpred']
        if predtype is 'sad':        
            sad_bins = get_log_bins([envpred_N])
            octave = range(0, len(sad_bins) - 1)
            site_pred = get_mete_sad(envpred_S, envpred_N, bin_edges=sad_bins)
            site_ids = [site['site_id'] for i in range(0, len(site_pred))]
            site_pred_with_id = DataFrame(np.column_stack([site_ids, octave, site_pred]),
                                          columns=['site_id', 'octave', 'env_pred'])    
        if predtype is 'rad':
            # note using observed S here for time being            
            rank = range(1, int(obs_S + 1))
            site_pred, p = get_mete_rad(obs_S, envpred_N)
            site_ids = [site['site_id'] for i in range(0, len(site_pred))]
            site_pred_with_id = DataFrame(np.column_stack([site_ids, rank, site_pred]),
                                          columns=['site_id', 'rank', 'env_pred'])
        if predtype is 'rare':
            pred_rad = get_mete_rad(int(envpred_S), envpred_N)[0]
            site_pred = sum([i <= 10 for i in pred_rad])
            site_ids = [site['site_id'] for i in range(0, len(site_pred))]
            site_pred_with_id = DataFrame(np.column_stack([site_ids, site_pred]),
                                          columns=['site_id', 'env_pred']) 
        envpred = envpred.append(site_pred_with_id, ignore_index=True)
    return envpred

if len(sys.argv) > 1:
    datasets = [sys.argv[1]]
else:
    datasets = ['bbs_2012', 'bbs_2008_2012', 'cbc', 'fia', 'gentry', 'nabc']

for dataset in datasets:
    for datatype in ['fit', 'test']:
        for predtype in ['sad', 'rad', 'rare']:
            envpred_data = read_csv('./results/' + dataset + '_state_var_' + datatype + '_obs_pred.csv')
            # check that predicted S smaller than predicted N
            rows_to_keep1 = (envpred_data.logNpred - envpred_data.logSpred) > 0
            rows_to_keep2 = (envpred_data.logNpred - envpred_data.logS) > 0
            rows_to_keep = rows_to_keep1 & rows_to_keep2
            envpred_data = envpred_data[rows_to_keep]
            envpred = get_envpred(envpred_data, predtype)
            envpred.to_csv('./results/' + dataset + '_' + predtype + '_' + datatype + '_pred.csv')
            print dataset + ' ' + datatype + ' ' + predtype + ' ' + 'complete!'

