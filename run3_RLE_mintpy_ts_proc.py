#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import glob
import time
import os
import subprocess
import matplotlib.pyplot as plt
from src.RLE_utils import simple_SBAS_stats, mintpy_SBAS_stats

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='MintPy time-series processing after pycuampcor offset tracking')
    parser.add_argument("--inputDir", dest="inputDir",
            default='./offsets', type=str, help='Input directory of the offsets with geotiff files of range/azimuth offsets and SNR (default: offsets)')
    parser.add_argument("--burst_id", dest='burst_id',
                         required=True, type=str, help='burst ID to be processed')    
    parser.add_argument("--snr", dest="snr",
            default=10, type=int, help='SNR threshold to be used for SBAS approach (default: 10)')
    parser.add_argument("--tsmethod", dest="tsmethod",
            default='mintpy', type=str, help='method for time-series inversion: mintpy (default), sbas (simple SBAS method)')
    parser.add_argument("--nprocs", dest="nprocs",
                         default=2, type=int, help='Number of processes to run (default: 2)')
    return parser.parse_args(args=iargs)

def main(inps):

    offset_dir = f'{inps.inputDir}/offsets'
    snr_th = inps.snr
    tsmethod = inps.tsmethod
    burst_id = inps.burst_id
    nprocs = inps.nprocs

    # Create necessary folders
    os.makedirs(f'{inps.inputDir}/mintpy', exist_ok=True)
    os.makedirs(f'{inps.inputDir}/summary', exist_ok=True)
    savedir = f'{inps.inputDir}/mintpy'

    list_rg_tif = glob.glob(f'{offset_dir}/*.rg_off.tif')  #geotiff of range offsets

    refDates = []
    secDates = []
    for _ in list_rg_tif:
        _rgtif = _.split('/')[-1]
        refDate = _rgtif[0:8]
        refDates.append(refDate)
        secDate = _rgtif[9:17]
        secDates.append(secDate)

    _ = {'ref': refDates, 'sec': secDates}
    df = pd.DataFrame.from_dict(_)

    days = refDates + secDates
    days = list(np.unique(sorted(days)))

    num_pairs = df.shape[0]   #number of pairs
    n_days = len(days)      #number of unique days

    print('')
    print(f'Number of pairs for offset tracking {num_pairs}\n')

    start_time = time.time()
    #applying SBAS approach
    list_rgoff = df['ref'] + '_' + df['sec'] + '.rg_off.tif'
    list_azoff = df['ref'] + '_' + df['sec'] + '.az_off.tif'
    list_snr = df['ref'] + '_' + df['sec'] + '.snr.tif'

    if ( tsmethod == 'sbas'):
        rg_avg, rg_std, _ = simple_SBAS_stats(list_rgoff,list_snr,savedir,snr_th)
        az_avg, az_std, _ = simple_SBAS_stats(list_azoff,list_snr,savedir,snr_th)
    else:
        rg_avg, rg_std, az_avg, az_std = mintpy_SBAS_stats(list_rgoff,list_azoff,list_snr,savedir,snr_th,0.25,nprocs) 

    end_time = time.time()
    time_taken = (end_time - start_time)/60.
    print(f'{time_taken} min taken for SBAS processing')

    df = None
    _ = {'date':days, 'rg_avg':rg_avg, 'rg_std':rg_std, 'az_avg':az_avg, 'az_std':az_std}
    df = pd.DataFrame.from_dict(_)
    df.to_csv(f'{inps.inputDir}/summary/RLE_{burst_id.upper()}.csv', index=False)

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()

    print("=========================================================================")
    print("Running Step 3 of the RLE: SBAS analysis")
    print("=========================================================================")

    # Run workflow
    main(inps)
