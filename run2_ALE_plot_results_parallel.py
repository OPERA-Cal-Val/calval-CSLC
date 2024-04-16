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
    parser.add_argument("--nprocs", dest="nprocs",
                         default=2, type=int, help='Number of processes to run (default: 2)')
    parser.add_argument("--ver", dest="prod_version",
                    default='v0.2', type=str, help='Product version to validate (default: v0.2)')
    return parser.parse_args(args=iargs)

def run_papermill(p):
    # Set Parameters
    burst_id = p[0]
    cr_network = p[1]
    save_dir = f'{p[-2]}/{cr_network}/{burst_id.upper()}'
    prod_version = p[-1]
    print(save_dir)
    print(burst_id)
    # Run the ALE for each date via papermill
    pm.execute_notebook('util_notebooks/plot_ALE_results_template.ipynb',
                f'{save_dir}/ipynbs/ALE_results_{cr_network}_{burst_id.upper()}.ipynb',
                parameters={'save_dir': save_dir,
                            'burst_id': burst_id,
                            'cr_network': cr_network},
                kernel_name='calval_CSLC')
    
    return (f'Finished plotting results for ({cr_network}) burst ({burst_id})')

def main(inps):
    # Specify valid burst(s)
    # Default is to loop through all
    sample_bursts = inps.burst_ids
    savedir = inps.savedir
    nprocs = inps.nprocs
    prod_version = inps.prod_version

    # read list of bursts used for validation
    validation_bursts = Path(f'validation_data/validation_bursts_priority_{prod_version}.csv')
    if validation_bursts.is_file():
        burstId_df = pd.read_csv(validation_bursts)
    else:
        raise Exception(f'Expected burst record {validation_bursts.absolute()} '
                        'not found. Check working directory.')

    # Start runtime evaluation
    start = timeit.default_timer()
    
    # Loop over all burst_ids
    params = []
    for burst_id in sample_bursts:
        print(burst_id)
        burst_cr_network = burstId_df[burstId_df['burst_id']==burst_id]['cr_network'].values[0]

        # Set parameters
        params.append([burst_id,burst_cr_network,savedir,prod_version])

    print(params)     
    print(f'Number of CPUs your computer have: {os.cpu_count()}')
    print(f'Using {nprocs} CPUs for this processing.')
    # Run papermill
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

    print("Running ALE now")
    
    # Run the main function
    main(inps)