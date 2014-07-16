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

from mete import get_mete_sad
from macroecotools import plot_color_by_pt_dens, obs_pred_rsquare

def import_data(resultdir, rawdir, datatype, dataset):
    envpred_data = read_csv(resultdir + dataset + '_state_var_' + datatype + '_obs_pred.csv')
    return envpred_data, sad_data

def get_log_bins(abu, logbase=2):
    max_bin = log(max(abu), logbase)
    if int(max_bin) == max_bin: # if is an integer then we need one additional bin
        max_bin = int(max_bin + 1)
    else:
        max_bin = int(ceil(max_bin))
    bins = [logbase ** x for x in range(0, max_bin + 1)]
    return bins

def get_envpred_sads(envpred_data):
    envpred_sads = DataFrame(columns=['site_id', 'octave', 'env_pred'])
    for index, site in envpred_data.iterrows():
        obs_S = site['S']
        envpred_S = 10 ** site['logSpred']
        envpred_N = 10 ** site['logNpred']
        #To produce a comparable number of species use obs_S; IS THIS RIGHT?
        sad_bins = get_log_bins([envpred_N])
        octave = range(0, len(sad_bins) - 1)
        site_sad = get_mete_sad(envpred_S, envpred_N, bin_edges=sad_bins)
        site_ids = [site['site_id'] for i in range(0, len(site_sad))]
        site_sad_with_id = DataFrame(np.column_stack([site_ids, octave, site_sad]),
                                     columns=['site_id', 'octave', 'env_pred'])
        envpred_sads = envpred_sads.append(site_sad_with_id, ignore_index=True)
    return envpred_sads

if len(sys.argv) > 1:
    datasets = [sys.argv[1]]
else:
    datasets = ['bbs_2012', 'bbs_2008_2012', 'cbc', 'fia', 'gentry', 'nabc']

for dataset in datasets:
    for datatype in ['fit', 'test']:
        envpred_data = read_csv('./results/' + dataset + '_state_var_' + datatype + '_obs_pred.csv')
        envpred_sads = get_envpred_sads(envpred_data)
        envpred_sads.to_csv('./results/' + dataset + '_sad_' + datatype + '_pred.csv')
        print dataset + ' ' + datatype + ' complete!' 

