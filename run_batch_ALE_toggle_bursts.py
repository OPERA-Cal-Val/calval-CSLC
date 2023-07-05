import os
import warnings
from pathlib import Path

import geopandas as gpd
import pandas as pd
import papermill as pm

import timeit

warnings.filterwarnings('ignore')

# Specify valid burst(s)
# Default is to loop through all
sample_bursts = ['t064_135523_iw2', 't071_151224_iw2']
savedir = '/Users/bato/work/OPERA/CSLC/CalVal/'
ovsFactor = 128

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

    for val_index, val_row in validation_bursts_df.iterrows():
        if val_row['burst_id'] == burstId:
            # Set parameters
            cslc_url = val_row['cslc_url']
            cslc_static_url = val_row['cslc_static_url']
            burst_id = val_row['burst_id']
            cr_network = burst_cr_network
            save_dir = f'{savedir}/{burst_cr_network}/{burst_id}'
            snr_threshold = 15
            solidtide = 'True'
            cslc_date = val_row['date']
            ovsFactor = ovsFactor

            # Create folders
            os.makedirs(f'{save_dir}/pngs', exist_ok=True)
            os.makedirs(f'{save_dir}/ipynbs', exist_ok=True)
            os.makedirs(f'{save_dir}/summary', exist_ok=True)
            os.makedirs(f'{save_dir}/crdata', exist_ok=True)

            # Check if file exist
            outcsv1 = f'ALE_{cr_network}_{burst_id}_allDates.csv'
            outcsv2 = f'ALE_{cr_network}_allCRs.csv'
            if os.path.isfile(outcsv1):
               os.remove(outcsv1)
            if os.path.isfile(outcsv2):
               os.remove(outcsv2)

            # Run the ALE for each date
            print(f'Processing AO ({cr_network}) burst ({burst_id}), for date ({cslc_date})')
            pm.execute_notebook('ALE_Stream_gamma.ipynb',
                        f'{save_dir}/ipynbs/ALE_COMPASS_{burst_id}_{cslc_date}.ipynb',
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
    
    # End runtime evaluation
    stop = timeit.default_timer()
    print(f'Finished run for {burst_cr_network}')
    print(f'Time: ', (stop - start)/60, 'min.')
