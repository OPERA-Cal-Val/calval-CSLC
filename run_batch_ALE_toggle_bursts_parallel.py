import os
import warnings
from pathlib import Path

import geopandas as gpd
import pandas as pd
import papermill as pm
import concurrent.futures
import timeit
warnings.filterwarnings('ignore')


def run_papermill(p):
    # Set Parameters
    cslc_date = p[0]
    burst_id = p[1]
    cslc_url = p[2]
    cslc_static_url = p[3]
    cr_network = p[4]
    snr_threshold = p[5]
    solidtide = p[6]
    ovsFactor = p[7]
    save_dir = f'{p[-1]}/{cr_network}/{burst_id}'

    # Run the ALE for each date via papermill
    pm.execute_notebook('ALE_template_gamma.ipynb',
                f'{save_dir}/ipynbs/ALE_{burst_id}_{cslc_date}.ipynb',
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

def main():
    # Specify valid burst(s)
    # Default is to loop through all
    sample_bursts = ['t064_135523_iw2', 't071_151224_iw2']
    savedir = '/Users/bato/work/OPERA/CSLC/CalVal/TESTPARALLEL'

    # read list of bursts used for validation
    validation_bursts = Path('validation_data/validation_bursts.csv')
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
    validation_csv = Path('validation_data/validation_table.csv')
    df = pd.read_csv(validation_csv)
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

        # Burst parameters
        snr_threshold = 15
        solidtide = 'True'
        ovsFactor = 128

        # Create folders
        os.makedirs(f'{savedir}/{burst_cr_network}/{burstId}/pngs', exist_ok=True)
        os.makedirs(f'{savedir}/{burst_cr_network}/{burstId}/ipynbs', exist_ok=True)
        os.makedirs(f'{savedir}/{burst_cr_network}/{burstId}/summary', exist_ok=True)
        os.makedirs(f'{savedir}/{burst_cr_network}/{burstId}/crdata', exist_ok=True)

        # Check if file exist
        outcsv1 = f'ALE_{burst_cr_network}_{burstId}_allDates.csv'
        outcsv2 = f'ALE_{burst_cr_network}_allCRs.csv'
        if os.path.isfile(outcsv1):
           os.remove(outcsv1)
        if os.path.isfile(outcsv2):
           os.remove(outcsv2)

        params = []
        for val_index, val_row in validation_bursts_df.iterrows():
            if val_row['burst_id'] == burstId:
                # Set parameters
                params.append([val_row['date'],val_row['burst_id'],val_row['cslc_url'],val_row['cslc_static_url'],burst_cr_network,snr_threshold,solidtide,ovsFactor,savedir])
        
        # # For debugging
        # print(list(map(testparallel,params)))
        
        print(f'Number of CPUs your computer has: {os.cpu_count()}')
        cpuworkers = int(os.cpu_count()/4)
        print(f'Using {cpuworkers} CPUs for this processing.')
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpuworkers) as executor:
            for result in executor.map(run_papermill,params):
                print(result)

        # End runtime evaluation
        stop = timeit.default_timer()
        print(f'Finished run for {burst_cr_network}')
        print(f'Time: ', (stop - start)/60, 'min.')


if __name__ == '__main__':
    main()