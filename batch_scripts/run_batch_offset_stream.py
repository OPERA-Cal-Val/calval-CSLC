#!/usr/bin/env python3
import argparse
import pandas as pd
import datetime as dt
import itertools
import numpy as np
import os, time
import subprocess
import matplotlib.pyplot as plt

import sys
sys.path.append('../')
from src.RLE_utils import hdf_stream, convert_to_slcvrt, array2raster, simple_SBAS_stats, mintpy_SBAS_stats

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='Batch processing of offset tracking and Generating a final output')
    parser.add_argument("--s3path", dest='s3path',
                         required=True, type=str, help='aws S3 bucket location (e.g., s3://opera-provisional/...)')
    parser.add_argument("--burstID", dest='bid',
                         required=True, type=str, help='burst ID to be processed')
    parser.add_argument("--datefile", dest='dfile',
                         required=True, type=str, help='text file with 1-column dates (YYYYMMDD) to be processed')
    parser.add_argument("--slc_dir", dest="slc_dir",
            default='SLCDIR', type=str, help='slc directory (default: SLCDIR)')
    parser.add_argument("--out_dir", dest="out_dir",
            default='outputs', type=str, help='output directory for offset results (default: outputs)')
    parser.add_argument("--minTemp", dest="minTemp",
            default=5, type=int, help='minimum temporal baseline (days) (default: 5)')
    parser.add_argument("--maxTemp", dest="maxTemp",
            default=370, type=int, help='maximum temporal baseline (days) (default: 370)')
    parser.add_argument("--ww", dest="ww",
            default=64, type=int, help='window width for offset tracking of pycuampcor (default: 64)')
    parser.add_argument("--wh", dest="wh",
            default=64, type=int, help='window height for offset tracking of pycuampcor (default: 64)')
    parser.add_argument("--nwdc", dest="nwdc",
            default=20, type=int, help='number of windows processed in a chunk along lines (default: 20)')
    parser.add_argument("--nwac", dest="nwac",
            default=20, type=int, help='number of windows processed in a chunk along columns (default: 20')
    parser.add_argument("--snr", dest="snr",
            default=10, type=int, help='SNR threshold to be used for SBAS approach (default: 10)')
    parser.add_argument("--tsmethod", dest="tsmethod",
            default='mintpy', type=str, help='method for time-series inversion: mintpy (default), sbas (simple SBAS method)')
    parser.add_argument("--pngfile", dest='png',
            default='RLE_ts.png',type=str, help='PNG file name for time-series RLE (default: RLE_ts.png)')
    parser.add_argument("--csvfile", dest='csv',             
            default='RLE_ts.csv',type=str, help='CSV file name for time-series RLE (default: RLE_ts.csv)')
    return parser.parse_args(args=iargs)

def run(inps):

    data_dir = inps.s3path
    burst_id = inps.bid

    slc_dir = inps.slc_dir
    out_dir = inps.out_dir
    f = open(inps.dfile)
    datels = f.read().splitlines()
    
    minDelta = inps.minTemp    #min time interval
    maxDelta = inps.maxTemp    #max time interval

    ##parameters for pycuampcor
    windowSizeWidth = inps.ww     #window size (width) for pycuampcor
    windowSizeHeight = inps.wh    #window size (height) for pycuampcor

    #parameters for GPU parallel processing
    numberWindowDownInChunk = inps.nwdc     #The number of windows processed in a batch/chunk, along lines
    numberWindowAcrossInChunk = inps.nwac   #The number of windows processed in a batch/chunk, along columns

    snr_th = inps.snr

    tsmethod = inps.tsmethod

    #parameters for gpu processing
    num_gpu = subprocess.getoutput('nvidia-smi --list-gpus | wc -l')
    num_gpu = int(num_gpu)
    print(f'number of GPU: {num_gpu} \n')

    #RLE with multi-reference pairs
    date_pair = list(itertools.combinations(datels,2))

    refDates = []
    secDates = []
    deltas = []
    gpuIDs = []
    _gpuID = 0

    for refDate, secDate in date_pair:
        delta = dt.datetime.strptime(secDate, "%Y%m%d") - dt.datetime.strptime(refDate, "%Y%m%d")
        delta = int(delta.days)
    
        if (delta > minDelta) & (delta < maxDelta):
            refDates.append(refDate)
            secDates.append(secDate)
            deltas.append(delta)
            gpuIDs.append(_gpuID)
            _gpuID = _gpuID + 1
            _gpuID = _gpuID % num_gpu

    _ = {'ref': refDates, 'sec': secDates, 'deltaT':deltas, 'gpuID':gpuIDs}
    df = pd.DataFrame.from_dict(_)
    print('slc pairs to be processed:')
    print(df)

    days = refDates + secDates
    days = list(np.unique(sorted(days)))

    num_pairs = df.shape[0]   #number of pairs
    n_days = len(days)      #number of unique days

    print(f'\nminimum temporal baseline: {minDelta} days')
    print(f'maximum temporal baseline: {maxDelta} days')
    print(f'number of pairs for offset tracking {num_pairs}\n')

    st_time = time.time()

    #creating slc inputs from COMPASS hdf file
    for day in days:
        outSLC = slc_dir + '/' + day + '.slc'
        outSLCvrt = outSLC + '.vrt'

        if os.path.isfile(outSLC) and (os.path.isfile(outSLCvrt)):
            print(f'{day}.slc exist. \n')
        else:
            path_h5 = f'{data_dir}/{burst_id}/{day}/{burst_id}_{day}.h5'   #path to COMPASS CSLC h5 file in aws s3 bucket

            xcoor, ycoor, dx, dy, epsg, slc, date = hdf_stream(path_h5)   
            convert_to_slcvrt(xcoor, ycoor, dx, dy, epsg, slc, date, slc_dir)   #generating slc with vrt 

    max_processes = num_gpu
    processes = set()

    #main offset tracking with pycuampcor
    for refd, secd, deviceID in zip(df['ref'],df['sec'],df['gpuID']):

        rgoff_file = out_dir + '/' + refd + '_' + secd + '.rg_off.tif'
        azoff_file = out_dir + '/' + refd + '_' + secd + '.az_off.tif'
        snr_file = out_dir + '/' + refd + '_' + secd + '.snr.tif'

        if os.path.isfile(rgoff_file) and os.path.isfile(azoff_file) and os.path.isfile(snr_file):
            print(f'{rgoff_file}, {azoff_file}, {snr_file} already exist \n')    #when files exist, skipping offset tracking
        else:
            cmd = f'python offset_pycuampcor.py --slc_dir {slc_dir} --dateref {refd} --datesec {secd} --deviceID {deviceID} --out_dir {out_dir} --ww {windowSizeWidth} --wh {windowSizeHeight} --nwdc {numberWindowDownInChunk} --nwac {numberWindowAcrossInChunk}'
            print(cmd)
            processes.add(subprocess.Popen(cmd.split(' ')))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update([p for p in processes if p.poll() is not None])

    #Check if all the child processes were closed
    for p in processes:
        if p.poll() is None:
            p.wait()

    end_time = time.time()
    time_taken = (end_time - st_time)/60.
    print(f'{time_taken} min taken for all pycuampcor processing')

    st_time = time.time()

    #applying SBAS approach
    list_rgoff = df['ref'] + '_' + df['sec'] + '.rg_off.tif'
    list_azoff = df['ref'] + '_' + df['sec'] + '.az_off.tif'
    list_snr = df['ref'] + '_' + df['sec'] + '.snr.tif'

    if ( tsmethod == 'sbas'):
        rg_avg, rg_std, _ = simple_SBAS_stats(list_rgoff,list_snr,out_dir,snr_th)
        az_avg, az_std, _ = simple_SBAS_stats(list_azoff,list_snr,out_dir,snr_th)
    else:
        rg_avg, rg_std, az_avg, az_std = mintpy_SBAS_stats(list_rgoff,list_azoff,list_snr,out_dir,snr_th) 

    end_time = time.time()
    time_taken = (end_time - st_time)/60.
    print(f'{time_taken} min taken for SBAS processing')

    df = None
    _ = {'date':days, 'rg_avg':rg_avg, 'rg_std':rg_std, 'az_avg':az_avg, 'az_std':az_std}
    df = pd.DataFrame.from_dict(_)
    df.to_csv(inps.csv, index=False)
    df['date'] = pd.to_datetime(df['date'])
    df['date'].dt.strftime('%Y%m%d')

    fig, ax = plt.subplots(2,1,figsize=(15,10),sharex=True)

    ax[0].set_title('RLE in Ground Range (m)')
    ax[0].axhspan(-0.5,0.5,color='red', alpha=0.05)    #OPERA requirements in ground range
    ax[0].errorbar(df['date'],df['rg_avg'],df['rg_std'],marker='o',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0)
    ax[0].set_ylim(-5,5)
    ax[0].grid(axis='x',linestyle='--')

    ax[1].set_title('RLE in Azimuth (m)')
    ax[1].axhspan(-0.75,0.75,color='red', alpha=0.05)    #OPERA requirements in azimuth
    ax[1].errorbar(df['date'],df['az_avg'],df['az_std'],marker='o',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0)
    ax[1].set_xlabel('dates')
    ax[1].set_ylim(-5,5)
    ax[1].grid(axis='x',linestyle='--')
    fig.savefig(inps.png,dpi=300,bbox_inches='tight')

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()
    
    # Run workflow
    run(inps)
