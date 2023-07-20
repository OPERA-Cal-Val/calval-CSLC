#!/usr/bin/env python3
import argparse
import datetime as dt
import os, time
from src.RLE_utils import hdf_stream, hdf_read, convert_to_slcvrt

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

    return parser.parse_args(args=iargs)

def run(inps):

    data_dir = inps.inputs
    burst_id = inps.bid

    slc_dir = inps.slc_dir
    f = open(inps.dfile)
    datels = f.read().splitlines()

    st_time = time.time()

    #creating slc inputs from COMPASS hdf file
    for day in datels:
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
                convert_to_slcvrt(xcoor, ycoor, dx, dy, epsg, slc, date, slc_dir)   #generating slc with vrt 

    end_time = time.time()
    time_taken = (end_time - st_time)/60.
    print(f'{time_taken} min taken for preparing CSLC before pycuampcor')

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()
    
    # Run workflow
    run(inps)
