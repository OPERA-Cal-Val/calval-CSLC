#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.pylab as pylab
from src.ALE_utils import enlos2rdr
from scipy.signal import detrend

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='plotting time-series offsets from pycuampcor')
    parser.add_argument("--savedir", dest='savedir', default='./RLE',             
            required=True,type=str, help='Path to the parent RLE directory. i.e. ./RLE')
    parser.add_argument("--burst_id", dest='burst_id',
                         required=True, type=str, help='burst ID to be processed')        
    parser.add_argument("--refDate", dest='refDate',
             default='20150601', type=str, help='Reference date of the stack')
    parser.add_argument("--detrend", dest='detrend',
             default=False, action="store_true",help='detrending from time-series offset (default: False)') 
    return parser.parse_args(args=iargs)

def if_pass(ts,requirement):
    bool_pass = (ts > -requirement) & (ts < requirement)
    pass_rate = np.count_nonzero(bool_pass)/len(bool_pass)
    return bool_pass, pass_rate, pass_rate>0.8

def main(inps):
    # Plotting Parameters
    params = {'legend.fontsize': 'x-large',
            'figure.figsize': (15, 5),
            'axes.labelsize': 'x-large',
            'axes.titlesize':'x-large',
            'xtick.labelsize':'x-large',
            'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)

    burst_id = inps.burst_id
    savedir = inps.savedir
    f = open(f'{savedir}/{burst_id.upper()}/cslc/enlos2rdr_{burst_id}.csv')
    los_east, los_north =f.read().split(' ')
    los_east = float(los_east); los_north = float(los_north)
    
    detrend_flag = inps.detrend

    # Read summary offset
    df = pd.read_csv(f'{savedir}/{burst_id.upper()}/summary/RLE_{burst_id.upper()}.csv') 
    df['date'] = pd.to_datetime(df['date'], format='%Y%m%d')
   
    grng, azi = enlos2rdr(df['rg_avg'],df['az_avg'], los_east, los_north) 
    grng_std, azi_std = enlos2rdr(df['rg_std'],df['az_std'], los_east, los_north) 

    if detrend_flag:
        grng = detrend(grng, type='linear')
        azi  = detrend(azi, type='linear') 
    
    df['grng_avg'] = grng
    df['azi_avg'] = azi
    df['grng_std'] = np.abs(grng_std)
    df['azi_std'] = np.abs(azi_std)

    # Get closest date to user-input reference date
    refDate_ = pd.to_datetime(inps.refDate, format='%Y%m%d')
    refDate = df.loc[(df['date']-refDate_).abs().idxmin(),'date']
 
    # Reference the stack to the reference date
    df['grng_avg'] = df['grng_avg'] - (df[df.date==refDate]['grng_avg'].values[0])
    df['azi_avg'] = df['azi_avg'] - (df[df.date==refDate]['azi_avg'].values[0])

    grng_avg = np.round(df['grng_avg'].mean(),3)
    grng_std = np.round(df['grng_avg'].std(),3)
    azi_avg = np.round(df['azi_avg'].mean(),3)
    azi_std = np.round(df['azi_avg'].std(),3)

    rg_bool_pass, rg_pass_rate, rg_pass_or_not = if_pass(df['grng_avg'],0.5)
    az_bool_pass, az_pass_rate, az_pass_or_not = if_pass(df['azi_avg'],0.75)

    # Plot time-series
    fig, ax = plt.subplots(2,1,figsize=(15,10),sharex=True)
    ax[0].set_ylabel(f'RLE in Ground Range (m)')
    ax[0].axhspan(-0.5,0.5,color='red', alpha=0.05,label='requirements')    #OPERA requirements in ground range
    ax[0].errorbar(df['date'][rg_bool_pass],df['grng_avg'][rg_bool_pass],df['grng_std'][rg_bool_pass],marker='o',color='b',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='passed offset')
    ax[0].errorbar(df['date'][~rg_bool_pass],df['grng_avg'][~rg_bool_pass],df['grng_std'][~rg_bool_pass],marker='o',color='r',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='failed offset')
    ax[0].set_ylim(-5,5)
    ax[0].grid(axis='x',linestyle='--')
    if rg_pass_or_not:
        ax[0].text(0.92,0.93,'Passed',color='w',size=18,weight='bold',transform = ax[0].transAxes,bbox=dict(facecolor='green',boxstyle='round',edgecolor='none'))
    else:
        ax[0].text(0.92,0.93,'Failed',color='w',size=18,weight='bold',transform = ax[0].transAxes,bbox=dict(facecolor='red',boxstyle='round',edgecolor='none'))
    ax[0].legend(loc = 'lower right',frameon=True)
    ax[0].text(0.02,0.05,f'RLE in Ground Range (m): {grng_avg}+/-{grng_std} \nReference Date: {refDate.date()}',transform = ax[0].transAxes)
    
    ax[1].set_ylabel(f'RLE in Azimuth (m)')
    ax[1].axhspan(-0.75,0.75,color='red', alpha=0.05,label='requirements')    #OPERA requirements in azimuth
    ax[1].errorbar(df['date'][az_bool_pass],df['azi_avg'][az_bool_pass],df['azi_std'][az_bool_pass],marker='o',color='b',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='passed offset')
    ax[1].errorbar(df['date'][~az_bool_pass],df['azi_avg'][~az_bool_pass],df['azi_std'][~az_bool_pass],marker='o',color='r',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='failed offset')
    ax[1].set_xlabel('Date')
    ax[1].set_ylim(-5,5)
    ax[1].grid(axis='x',linestyle='--')
    if az_pass_or_not:
        ax[1].text(0.92,0.93,'Passed',color='w',size=18,weight='bold',transform = ax[1].transAxes,bbox=dict(facecolor='green',boxstyle='round',edgecolor='none'))
    else:
        ax[1].text(0.92,0.93,'Failed',color='w',size=18,weight='bold',transform = ax[1].transAxes,bbox=dict(facecolor='red',boxstyle='round',edgecolor='none'))
    ax[1].legend(loc = 'lower right',frameon=True)
    ax[1].text(0.02,0.05,f'RLE in Azimuth (m): {azi_avg}+/-{azi_std} \nReference Date: {refDate.date()}',transform = ax[1].transAxes)
    plt.tight_layout()
    savefn = f'{savedir}/{burst_id.upper()}/summary/RLE_{burst_id.upper()}.png'
    fig.savefig(savefn,dpi=300,bbox_inches='tight')

    # Plot Azimuth vs Ground Range
    fig, ax = plt.subplots(1,1,figsize=(10,10),sharex=True)
    rle_req = patches.Rectangle((0.5, 0.75), -1, -1.5, linewidth=2, edgecolor='r', facecolor='none',label='requirements')  #OPERA ALE requirements
    ax.errorbar(df['grng_avg'][rg_bool_pass & az_bool_pass],df['azi_avg'][rg_bool_pass & az_bool_pass],xerr=df['grng_std'][rg_bool_pass & az_bool_pass],yerr=df['azi_std'][rg_bool_pass & az_bool_pass],marker='o',color='b',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, alpha=0.6, zorder=0,label='passed offset')
    ax.errorbar(df['grng_avg'][~rg_bool_pass | ~az_bool_pass],df['azi_avg'][~rg_bool_pass | ~az_bool_pass],xerr=df['grng_std'][~rg_bool_pass | ~az_bool_pass],yerr=df['azi_std'][~rg_bool_pass | ~az_bool_pass],marker='o',color='r',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, alpha=0.6, zorder=0,label='failed offset')
    if rg_pass_or_not and az_pass_or_not:
        ax.text(0.85,0.95,'Passed',color='w',size=18,weight='bold',transform = ax.transAxes,bbox=dict(facecolor='green',boxstyle='round',edgecolor='none'))
    else:
        ax.text(0.85,0.95,'Failed',color='w',size=18,weight='bold',transform = ax.transAxes,bbox=dict(facecolor='red',boxstyle='round',edgecolor='none'))
    ax.grid(True)
    ax.set_xlim(-2,2)
    ax.set_ylim(-2,2)
    ax.axhline(0, color='black')
    ax.axvline(0, color='black')
    ax.set_title(f'$\Delta gr$: {grng_avg} +/- {grng_std} m, \
    $\Delta az$: {azi_avg}, +/- {azi_std} m')
    ax.set_xlabel('Ground Range Offset (m)')
    ax.set_ylabel('Azimuth Offset (m)')
    nl = '\n'
    fig.suptitle(f"{burst_id.upper()} 2015-2023 {nl} (Average RLE for each date)", fontsize=16)
    ax.add_patch(rle_req)
    ax.legend(loc = 'lower right',frameon=True)
    ax.text(0.02,0.05,f'Reference Date: {refDate.date()}',transform = ax.transAxes)
    fig.savefig(f'{savedir}/{burst_id.upper()}/summary/RLE_summary_{burst_id.upper()}_grazi_mean.png',dpi=300,bbox_inches='tight')

    # Save the RLE Summary Rating as CSV
    summary_df = pd.DataFrame()
    summary_df['burst_id'] = None
    summary_df['rating'] = None
    summary_df['ref_date'] = None
    summary_df['GRng_Mean'] = None
    summary_df['GRng_Stdev'] = None
    summary_df['Azi_Mean'] = None
    summary_df['Azi_Stdev'] = None

    burst_ids = []; rle_rating = []; ref_date = []; GRng_mean = []; Az_mean = []; GRng_stdev = []; Az_stdev = []

    if rg_pass_or_not and az_pass_or_not:
        rating = 'PASS'
    else:
        rating = 'FAIL'

    burst_ids.append(burst_id.upper())
    rle_rating.append(rating)
    ref_date.append(dt.datetime.strftime(refDate,'%Y%m%d'))
    GRng_mean.append(grng_avg)
    GRng_stdev.append(grng_std)
    Az_mean.append(azi_avg)
    Az_stdev.append(azi_std)

    summary_df['burst_id'] = burst_ids
    summary_df['rating'] = rle_rating
    summary_df['ref_date'] = ref_date
    summary_df['GRng_Mean'] = GRng_mean
    summary_df['GRng_Stdev'] = GRng_stdev
    summary_df['Azi_Mean'] = Az_mean
    summary_df['Azi_Stdev'] = Az_stdev

    summary_df.to_csv(f'{savedir}/RLE_summary.csv', index=False, mode='a', header=False)

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()
    
    # Run workflow
    main(inps)