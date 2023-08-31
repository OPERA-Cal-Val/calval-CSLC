#!/usr/bin/env python3

import argparse
import os
import warnings
from pathlib import Path
import datetime as dt
import requests
import time

import geopandas as gpd
import pandas as pd
import papermill as pm
import concurrent.futures
import timeit
warnings.filterwarnings('ignore')

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='Absolute Location Error estimation to validate OPERA CSLCs')
    parser.add_argument("--savedir", dest='savedir',
                         required=True, type=str, help='Save directory')
    parser.add_argument("--burst_ids", dest="burst_ids",
                         required=True, nargs='+', help="List of burst_ids to process, ['t064_135523_iw2', 't071_151224_iw2'] ")
    parser.add_argument("--startDate", dest="startDate",
                         default='20140101', type=str, help="Start date of ALE evaluation")
    parser.add_argument("--endDate", dest="endDate",
                         default=dt.datetime.today().strftime('%Y%m%d'), type=str, help="Start date of ALE evaluation")
    parser.add_argument("--snr", dest="snr",
                         default=15, type=int, help='Signal-to-noise ratio threshold to remove corner reflectors outliers (default:15)')    
    parser.add_argument("--ovsFactor", dest="ovsFactor",
                         default=128, type=int, help='Oversampling factor size to locate the peak amplitude (default: 128)')
    parser.add_argument("--nprocs", dest="nprocs",
                         default=2, type=int, help='Number of processes to run (default: 2)')
    parser.add_argument("--validation_bursts", dest="validation_bursts",
                        default=Path('validation_data/validation_bursts.csv'), type=str, help='Validation burst table (default: validation_data/validation_bursts.csv)')
    parser.add_argument("--validation_csv", dest="validation_csv",
                        default=Path('validation_data/validation_table.csv'), type=str, help='Validation table (default: validation_data/validation_table.csv')
    return parser.parse_args(args=iargs)

def run_papermill(p):
    # Set Parameters
    cslc_date = p[0]
    burst_id = p[1]
    cslc_url = p[2]
    cslc_static_url = p[3]
    cr_network = p[4]
    snr_threshold = p[5]
    if cr_network == 'Rosamond':
        solidtide='True'
    else:
        solidtide='False'
    ovsFactor = p[6]
    
    save_dir = f'{p[-1]}/{cr_network}/{burst_id.upper()}'

    # Run the ALE for each date via papermill
    pm.execute_notebook('./util_notebooks/ALE_template.ipynb',
                f'{save_dir}/ipynbs/ALE_{burst_id.upper()}_{cslc_date}.ipynb',
                parameters={'cslc_url': cslc_url,
                            'cslc_static_url': cslc_static_url,
                            'save_dir': save_dir,
                            'burst_id': burst_id,
                            'date': cslc_date, 
                            'snr_threshold': snr_threshold, 
                            'solidtide': solidtide,
                            'cr_network': cr_network,
                            'ovsFactor': ovsFactor},
                kernel_name='calval_CSLC')
    
    return (f'Finished processing AO ({cr_network}) burst ({burst_id}), for date ({cslc_date})')

def download_crdata(p):
    # Set Parameters
    date = p[0]
    burst_id = p[1]
    cslc_url = p[2]
    cr_network = p[4]
    save_dir = f'{p[-1]}/{cr_network}/{burst_id.upper()}'
    sensing_time = dt.datetime.strptime(cslc_url.split('/')[-1].split('_')[-3], '%Y%m%dT%H%M%SZ').strftime('%Y-%m-%d+%H\u0021%M')

    print(f'Downloading crdata_{burst_id.upper()}_{date}.csv')
    # Get date and download
    if cr_network=='Rosamond':
        # Download corner reflector data from UAVSAR/NISAR based on the date of the CSLC product
        res = requests.get(f'https://uavsar.jpl.nasa.gov/cgi-bin/corner-reflectors.pl?date={sensing_time}&project=uavsar')
    elif cr_network=='Oklahoma':
        # Download corner reflector data from UAVSAR/NISAR based on the date of the CSLC product
        res = requests.get(f'https://uavsar.jpl.nasa.gov/cgi-bin/corner-reflectors.pl?date={sensing_time}&project=nisar')
    elif cr_network=='Alaska':
        # Download corner reflector data from NISAR based on the date of the CSLC product
        res = requests.get(f'https://uavsar.jpl.nasa.gov/cgi-bin/corner-reflectors.pl?date={sensing_time}&project=alaska')
    else:
        raise SystemExit(f'No corner reflector data found for {burst_id}_{date}. Terminating process.')
    
    # Write the crdata to file
    with open(f'{save_dir}/crdata/crdata_{burst_id.upper()}_{date}.csv', 'wb') as crfile:
        crfile.write(res.content)
        crfile.flush()
    
    return

def main(inps):
    # Specify valid burst(s)
    # Default is to loop through all
    sample_bursts = inps.burst_ids
    savedir = inps.savedir
    snr_threshold = inps.snr
    ovsFactor = inps.ovsFactor
    nprocs = inps.nprocs
    startDate = inps.startDate
    endDate = inps.endDate

    # read list of bursts used for validation
    validation_bursts = Path(inps.validation_bursts) #Path('validation_data/validation_bursts.csv')
    if validation_bursts.is_file():
        burstId_df = pd.read_csv(validation_bursts)
    else:
        raise Exception(f'Expected burst record {validation_bursts.absolute()} '
                        'not found. Check working directory.')

    # only pass records matching specied AOI(s)
    if sample_bursts == []:
        sample_bursts =  burstId_df['burst_id'].unique().tolist()
    else:
        burstId_df = burstId_df[burstId_df['burst_id'].isin(sample_bursts)]

    # access table of all S3 links
    validation_csv = Path(inps.validation_csv) #Path('validation_data/validation_table.csv')
    df_ = pd.read_csv(validation_csv)
    df = df_.drop_duplicates(subset=['burst_id', 'date'])

    validation_bursts_df = gpd.GeoDataFrame(
        df.loc[:, [c for c in df.columns if c != "geometry"]],
        geometry=gpd.GeoSeries.from_wkt(df["geometry"])
        )

    # Loop over all valid bursts for specified AOIs
    for burst_index, burst_row in burstId_df.iterrows():
        burstId = burst_row['burst_id']
        burst_cr_network = burst_row['cr_network']
        # Start runtime evaluation
        start = timeit.default_timer()

        # Create folders
        os.makedirs(f'{savedir}/{burst_cr_network}/{burstId.upper()}/pngs', exist_ok=True)
        os.makedirs(f'{savedir}/{burst_cr_network}/{burstId.upper()}/ipynbs', exist_ok=True)
        os.makedirs(f'{savedir}/{burst_cr_network}/{burstId.upper()}/summary', exist_ok=True)
        os.makedirs(f'{savedir}/{burst_cr_network}/{burstId.upper()}/crdata', exist_ok=True)

        # Check if file exist
        outcsv1 = f'ALE_{burst_cr_network}_{burstId}_allDates.csv'
        outcsv2 = f'ALE_{burst_cr_network}_allCRs.csv'
        if os.path.isfile(outcsv1):
           os.remove(outcsv1)
        if os.path.isfile(outcsv2):
           os.remove(outcsv2)

        params = []
        for val_index, val_row in validation_bursts_df.iterrows():
            if (val_row['burst_id'] == burstId) and (dt.datetime.strptime(str(val_row['date']),'%Y%m%d') >= dt.datetime.strptime(startDate,'%Y%m%d')) \
                and (dt.datetime.strptime(str(val_row['date']),'%Y%m%d') <= dt.datetime.strptime(endDate,'%Y%m%d')):
                # Set parameters
                params.append([val_row['date'],val_row['burst_id'],val_row['cslc_url'],val_row['cslc_static_url'],burst_cr_network,snr_threshold,ovsFactor,savedir])
        
        print(f'Number of CPUs your computer have: {os.cpu_count()}')
        print(f'Using {nprocs} CPUs for this processing.')

        # Download CR data first
        with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
            executor.map(download_crdata,params)

        # Run papermill
        time.sleep(5)
        with concurrent.futures.ProcessPoolExecutor(max_workers=nprocs) as executor:
            for result in executor.map(run_papermill,params):
                print(result)

        # End runtime evaluation
        stop = timeit.default_timer()
        print(f'Finished run for {burst_cr_network}')
        print(f'Time: ', (stop - start)/60, 'min.')


if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()

    print("==================================================================")
    print("Running the Absolute Location Estimation via Point-Target Analysis")
    print("==================================================================")
    
    # Run the main function
    main(inps)