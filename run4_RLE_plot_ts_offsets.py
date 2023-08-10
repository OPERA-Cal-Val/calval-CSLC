#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from src.ALE_utils import en2rdr
from scipy.signal import detrend

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='plotting time-series offsets from pycuampcor')
    parser.add_argument("--csvfile", dest='csv',             
            required=True,type=str, help='CSV file name for time-series RLE')
    parser.add_argument("--pngfile", dest='png',
            required=True,type=str, help='PNG file name for time-series RLE')
    parser.add_argument("--en2rdr", dest='en2rdr',
             required=True,type=str, help='csv file of incidence and azimuth angle (unit: deg)') 
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

    f = open(inps.en2rdr)
    inc_angle, az_angle =f.read().split(' ')
    inc_angle = float(inc_angle); az_angle = float(az_angle)

    df = pd.read_csv(inps.csv) 
    df['date'] = pd.to_datetime(df['date'], format='%Y%m%d')

    grng, azi = en2rdr(df['rg_avg'],df['az_avg'], az_angle, inc_angle) 
    grng_std, azi_std = en2rdr(df['rg_std'],df['az_std'], az_angle, inc_angle) 

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

    # Evaluate if the point is passing the requirement or not
    rg_bool_pass, rg_pass_rate, rg_pass_or_not = if_pass(df['grng_avg'],0.5)
    az_bool_pass, az_pass_rate, az_pass_or_not = if_pass(df['azi_avg'],0.75)

    # Plot figure
    fig, ax = plt.subplots(2,1,figsize=(15,10))
    ax[0].set_ylabel(f'RLE in Ground Range (m) \n Reference Date: {refDate.date()}')
    ax[0].axhspan(-0.5,0.5,color='red', alpha=0.05,label='requirements')    #OPERA requirements in ground range
    ax[0].errorbar(df['date'][rg_bool_pass],df['grng_avg'][rg_bool_pass],df['grng_std'][rg_bool_pass],marker='o',color='b',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='passed offset')
    ax[0].errorbar(df['date'][~rg_bool_pass],df['grng_avg'][~rg_bool_pass],df['grng_std'][~rg_bool_pass],marker='o',color='r',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='failed offset')
    ax[0].set_ylim(-5,5)
    # ax[0].set_xlabel('Date')
    ax[0].grid(axis='x',linestyle='--')
    if rg_pass_or_not:
        ax[0].text(0.94,0.90,'Pass',color='w',size=15,weight='bold',transform = ax[0].transAxes,bbox=dict(facecolor='blue',boxstyle='round',edgecolor='none'))
    else:
        ax[0].text(0.94,0.90,'Fail',color='w',size=15,weight='bold',transform = ax[0].transAxes,bbox=dict(facecolor='red',boxstyle='round',edgecolor='none'))
    ax[0].legend(loc = 'lower right',frameon=True)
    
    ax[1].set_ylabel(f'RLE in Azimuth (m) \n Reference Date: {refDate.date()}')
    ax[1].axhspan(-0.75,0.75,color='red', alpha=0.05,label='requirements')    #OPERA requirements in azimuth
    ax[1].errorbar(df['date'][az_bool_pass],df['azi_avg'][az_bool_pass],df['azi_std'][az_bool_pass],marker='o',color='b',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='passed offset')
    ax[1].errorbar(df['date'][~az_bool_pass],df['azi_avg'][~az_bool_pass],df['azi_std'][~az_bool_pass],marker='o',color='r',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='failed offset')
    ax[1].set_xlabel('Date')
    ax[1].set_ylim(-5,5)
    ax[1].grid(axis='x',linestyle='--')
    if az_pass_or_not:
        ax[1].text(0.94,0.90,'Pass',color='w',size=15,weight='bold',transform = ax[1].transAxes,bbox=dict(facecolor='blue',boxstyle='round',edgecolor='none'))
    else:
        ax[1].text(0.94,0.90,'Fail',color='w',size=15,weight='bold',transform = ax[1].transAxes,bbox=dict(facecolor='red',boxstyle='round',edgecolor='none'))
    ax[1].legend(loc = 'lower right',frameon=True)
    plt.tight_layout()
    fig.savefig(inps.png,dpi=300,bbox_inches='tight')

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()
    
    # Run workflow
    main(inps)
