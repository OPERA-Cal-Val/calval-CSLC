#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='plotting time-series offsets from pycuampcor')
    parser.add_argument("--csvfile", dest='csv',             
            required=True,type=str, help='CSV file name for time-series RLE')
    parser.add_argument("--pngfile", dest='png',
            required=True,type=str, help='PNG file name for time-series RLE')
    return parser.parse_args(args=iargs)

def if_pass(ts,requirement):
    bool_pass = (ts > -requirement) & (ts < requirement)
    pass_rate = np.count_nonzero(bool_pass)/len(bool_pass)
    return bool_pass, pass_rate, pass_rate>0.8

def run(inps):

    df = pd.read_csv(inps.csv) 
    df['date'] = pd.to_datetime(df['date'], format='%Y%m%d')

    rg_bool_pass, rg_pass_rate, rg_pass_or_not = if_pass(df['rg_avg'],0.5)
    az_bool_pass, az_pass_rate, az_pass_or_not = if_pass(df['az_avg'],0.75)

    fig, ax = plt.subplots(2,1,figsize=(15,10),sharex=True)
    ax[0].set_title('RLE in Ground Range (m)')
    ax[0].axhspan(-0.5,0.5,color='red', alpha=0.05,label='requirements')    #OPERA requirements in ground range
    ax[0].errorbar(df['date'][rg_bool_pass],df['rg_avg'][rg_bool_pass],df['rg_std'][rg_bool_pass],marker='o',color='b',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='passed offset')
    ax[0].errorbar(df['date'][~rg_bool_pass],df['rg_avg'][~rg_bool_pass],df['rg_std'][~rg_bool_pass],marker='o',color='r',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='failed offset')
    ax[0].set_ylim(-5,5)
    ax[0].grid(axis='x',linestyle='--')
    if rg_pass_or_not:
        ax[0].text(0.94,0.90,'Pass',color='w',size=15,weight='bold',transform = ax[0].transAxes,bbox=dict(facecolor='blue',boxstyle='round',edgecolor='none'))
    else:
        ax[0].text(0.94,0.90,'Fail',color='w',size=15,weight='bold',transform = ax[0].transAxes,bbox=dict(facecolor='red',boxstyle='round',edgecolor='none'))
    ax[0].legend(loc = 'lower right',frameon=True)
    
    ax[1].set_title('RLE in Azimuth (m)')
    ax[1].axhspan(-0.75,0.75,color='red', alpha=0.05,label='requirements')    #OPERA requirements in azimuth
    ax[1].errorbar(df['date'][az_bool_pass],df['az_avg'][az_bool_pass],df['az_std'][az_bool_pass],marker='o',color='b',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='passed offset')
    ax[1].errorbar(df['date'][~az_bool_pass],df['az_avg'][~az_bool_pass],df['az_std'][~az_bool_pass],marker='o',color='r',linestyle=' ',ecolor='lightgray', elinewidth=3, capsize=0, zorder=0,label='failed offset')
    ax[1].set_xlabel('dates')
    ax[1].set_ylim(-5,5)
    ax[1].grid(axis='x',linestyle='--')
    if az_pass_or_not:
        ax[1].text(0.94,0.90,'Pass',color='w',size=15,weight='bold',transform = ax[1].transAxes,bbox=dict(facecolor='blue',boxstyle='round',edgecolor='none'))
    else:
        ax[1].text(0.94,0.90,'Fail',color='w',size=15,weight='bold',transform = ax[1].transAxes,bbox=dict(facecolor='red',boxstyle='round',edgecolor='none'))
    ax[1].legend(loc = 'lower right',frameon=True)
    fig.savefig(inps.png,dpi=300,bbox_inches='tight')

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()
    
    # Run workflow
    run(inps)
