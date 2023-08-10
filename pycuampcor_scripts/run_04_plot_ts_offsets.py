#!/usr/bin/env python3
import argparse
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='plotting time-series offsets from pycuampcor')
    parser.add_argument("--csvfile", dest='csv',             
            type=str, help='CSV file name for time-series RLE')
    parser.add_argument("--pngfile", dest='png',
            type=str, help='PNG file name for time-series RLE')
    
    return parser.parse_args(args=iargs)

def run(inps):

    df = pd.read_csv(inps.csv) 
    df['date'] = pd.to_datetime(df['date'], format='%Y%m%d')

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
