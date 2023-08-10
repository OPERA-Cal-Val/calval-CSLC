#!/usr/bin/env python3
import argparse
import datetime as dt
import os, time
import numpy as np
from multiprocessing import Process

import sys
sys.path.append('../')
from src.RLE_utils import hdf_stream, hdf_read, convert_to_slcvrt

def fetch_cslc(day):
    outSLC = slc_dir + '/' + day + '.slc'
    outSLCvrt = outSLC + '.vrt'

    if os.path.isfile(outSLC) and (os.path.isfile(outSLCvrt)):
        print(f'{day}.slc exist. \n')
    else:
        path_h5 = f'{data_dir}/{burst_id}/{day}/{burst_id}_{day}.h5'   #path to COMPASS CSLC h5 file
        cslc_source = data_dir[0:2]
        if (cslc_source == 's3'):
            xcoor, ycoor, dx, dy, epsg, slc, date = hdf_stream(path_h5)   
            convert_to_slcvrt(xcoor, ycoor, dx, dy, epsg, slc, date, slc_dir)   #generating slc with vrt
        else:
            xcoor, ycoor, dx, dy, epsg, slc, date = hdf_read(path_h5)
            convert_to_slcvrt(xcoor, ycoor, dx, dy, epsg, slc, date, slc_dir)

def multi_stream(days,n_multi):

    n_dates = len(days)    #number of all dates
    n_group = int(np.ceil(n_dates/n_multi))   #number of groups for multiprocessing

    st_time = time.time()
    print('Streaming started')

    for i in range(n_group):
         _i_min = i*n_multi
         _i_max = min(_i_min+n_multi,n_dates)

         _days = days[_i_min:_i_max] 

         print(f'{_days} to be downloaded')

         processes = [Process(target=fetch_cslc, args=(day,)) for day in _days]

         for process in processes:
             process.start()
         for process in processes:  
             process.join()

         processes = None

    end_time = time.time()
    time_taken = (end_time - st_time)/60.
    print(f'{time_taken} min taken for multi-streaming')

def check_missing(days):
    #checking if all files are streamed correctly
    missing_dates = []

    for _ in days:
        outSLC = slc_dir + '/' + _ + '.slc'
        outSLCvrt = outSLC + '.vrt'

        if os.path.isfile(outSLC) and os.path.isfile(outSLCvrt):
            pass
        else:
            missing_dates.append(_)

    n_miss = len(missing_dates)    
 
    if not missing_dates:
        pass
    else:
        print(f'missing dates: {missing_dates}')
        print('Warning: missing dates are streamed...')
        multi_stream(missing_dates,n_multi)
        print('finished...') 

    return n_miss

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='Reading CSLC hdf files in AWS S3 bucket or local directory and exporting tiff files before pycuampcor')
    parser.add_argument("--inputLoc", dest='inputs',
                         required=True, type=str, help='e.g., aws S3 bucket location (e.g., s3://opera-provisional/...)')
    parser.add_argument("--burstID", dest='bid',
                         required=True, type=str, help='burst ID to be processed')
    parser.add_argument("--datefile", dest='dfile',
                         required=True, type=str, help='text file with 1-column dates (YYYYMMDD) to be processed')
    parser.add_argument("--slc_dir", dest="slc_dir",
            default='SLCDIR', type=str, help='slc directory (default: SLCDIR)')
    parser.add_argument("--n_multi", dest="n_multi",
            default=4, type=int, help='number of multi streaming/downloading (default: 4)')

    return parser.parse_args(args=iargs)

def run(inps):

    global data_dir, burst_id, slc_dir, n_multi
    data_dir = inps.inputs
    burst_id = inps.bid

    slc_dir = inps.slc_dir
    f = open(inps.dfile)
    datels = f.read().splitlines()

    n_multi = inps.n_multi

    print(f'input location: {data_dir}')
    print(f'burst ID: {burst_id}')
    print(f'output SLC directory {slc_dir}')
    print(f'number of multi processing {n_multi}\n')
  
    multi_stream(datels,n_multi)  #streaming via multiprocessing

    #downloading missed dates
    n_iter = 100
    
    for ii in range(n_iter):
        n_miss = check_missing(datels)
  
        if (n_miss == 0):
            break

    if (n_miss == 0):
        print('all files are successfully streamed')
    else:
        print('still missing dates to be downloaded')

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()
    
    # Run workflow
    run(inps)
