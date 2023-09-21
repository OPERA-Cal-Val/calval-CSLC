###!/usr/bin/env python3

import argparse
import os
import warnings
from pathlib import Path
import datetime as dt
import time
import numpy as np

import fsspec

import geopandas as gpd
import pandas as pd
import concurrent.futures
import timeit
warnings.filterwarnings('ignore')
from src.RLE_utils import stream_cslc, convert_to_slcvrt, stream_static_layers

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='creating geotiff files for relative geolocation error estimation via streaming')
    parser.add_argument("--savedir", dest='savedir',
                         required=True, type=str, help='Save directory')
    parser.add_argument("--burst_ids", dest="burst_ids",
                         required=True, nargs='+', help="List of burst_ids to process, ['t064_135523_iw2', 't071_151224_iw2'] ")
    parser.add_argument("--startDate", dest="startDate",
                         default='20140101', type=str, help="Start date of RLE evaluation (default: 20140101)")
    parser.add_argument("--endDate", dest="endDate",
                         default=dt.datetime.today().strftime('%Y%m%d'), type=str, help="End date of RLE evaluation (default: today)")
    parser.add_argument("--nprocs", dest="nprocs",
                         default=2, type=int, help='Number of processes to run (default: 2)')
    parser.add_argument("--validation_bursts", dest="validation_bursts",
                        default=Path('validation_data/validation_bursts.csv'), type=str, help='Validation burst table (default: validation_data/validation_bursts.csv)')
    parser.add_argument("--validation_csv", dest="validation_csv",
                        default=Path('validation_data/validation_table.csv'), type=str, help='Validation table (default: validation_data/validation_table.csv')
    return parser.parse_args(args=iargs)

def cslc2tiff(p):
    date = p[0]
    burst_id = p[1]
    cslc_url = p[2]
    save_dir = f'{p[-1]}/{burst_id.upper()}/cslc'  

    print(f'Product: {cslc_url}')
    cslc,xcoor,ycoor,dx,dy,epsg,sensing_start,sensing_stop,dims,bounding_polygon,orbit_direction,center_lon,center_lat = stream_cslc(cslc_url)
    convert_to_slcvrt(xcoor, ycoor, dx, dy, epsg, cslc, date, save_dir)   #generating slc with vrt
    
    return f'OPERA CSLC with burst_id ({burst_id}) for date ({date}) successfully stored in ({save_dir})'

def main(inps):
    # Specify valid burst(s)
    # Default is to loop through all
    sample_bursts = inps.burst_ids
    savedir = inps.savedir
    nprocs = inps.nprocs
    startDate = inps.startDate
    endDate = inps.endDate
    # valBursts = inps.valBursts
    # valTable = inps.valTable

    # read list of bursts used for validation
    validation_bursts = Path(inps.validation_bursts) #Path(valBursts)
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
    validation_csv = Path(inps.validation_csv) #Path(valTable)
    df_ = pd.read_csv(validation_csv)
    df = df_.drop_duplicates(subset=['burst_id', 'date'])

    validation_bursts_df = gpd.GeoDataFrame(
        df.loc[:, [c for c in df.columns if c != "geometry"]],
        geometry=gpd.GeoSeries.from_wkt(df["geometry"])
        )

    # Loop over all valid bursts for specified AOIs
    for burst_index, burst_row in burstId_df.iterrows():
        burstId = burst_row['burst_id']

        # Start runtime evaluation
        start = timeit.default_timer()

        # Create folders
        os.makedirs(f'{savedir}/{burstId.upper()}/cslc',exist_ok=True)
        
        params = []
        for val_index, val_row in validation_bursts_df.iterrows():
            if (val_row['burst_id'] == burstId) and (dt.datetime.strptime(str(val_row['date']),'%Y%m%d') >= dt.datetime.strptime(startDate,'%Y%m%d')) \
                and (dt.datetime.strptime(str(val_row['date']),'%Y%m%d') <= dt.datetime.strptime(endDate,'%Y%m%d')):

                enlos2rdr = f'{savedir}/{burstId.upper()}/cslc/enlos2rdr_{burstId}.csv'  #enlos2rdr file (los_east, los_north) for converting EN to RDR

                path_enlos2rdr = Path(enlos2rdr) 
                if path_enlos2rdr.is_file():
                    pass
                else:
                    #reading static layer
                    los_east, los_north = stream_static_layers(val_row['cslc_static_url'])
                    los_east = np.nanmean(los_east)
                    los_north = np.nanmean(los_north)
                    with open(enlos2rdr,'w') as f:
                        f.write(f'{los_east} {los_north}')

                # Set parameters
                params.append([val_row['date'],val_row['burst_id'],val_row['cslc_url'],savedir])

        print(f'Number of CPUs your computer have: {os.cpu_count()}')
        print(f'Using {nprocs} CPUs for this processing.')

        # Run cslc2tiff
        with concurrent.futures.ProcessPoolExecutor(max_workers=nprocs) as executor:
            for result in executor.map(cslc2tiff,params):
                print(result)

        # End runtime evaluation
        stop = timeit.default_timer()
        print(f'Finished run for {burstId}')
        print(f'Time: ', (stop - start)/60, 'min.')

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()

    print("=========================================================================")
    print("Running Step 1 of the RLE: Streaming CSLC and saving to PyCuAmpcor format")
    print("=========================================================================")
    
    # Run the main function
    main(inps)

