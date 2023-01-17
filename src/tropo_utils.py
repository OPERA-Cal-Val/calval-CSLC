import numpy as np
import datetime as dt
import os

import re
import math

import pyaps3 as pa
from mintpy.utils import ptime, readfile, writefile, utils as ut

WEATHER_MODEL_HOURS = {
    'ERA5'   : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
    'ERAINT' : [0, 6, 12, 18],
    'MERRA'  : [0, 6, 12, 18],
}

def get_grib_filenames(date_list, hour, model, grib_dir, snwe=None):
    """Get default grib file names based on input info.
    Parameters: date_list  - list of str, date in YYYYMMDD format
                hour       - str, hour in 2-digit with zero padding
                model      - str, global atmospheric model name
                grib_dir   - str, local directory to save grib files
                snwe       - tuple of 4 int, for ERA5 only.
    Returns:    grib_files - list of str, local grib file path
    """
    # area extent
    area = snwe2str(snwe)

    grib_files = []
    
    if model == 'ERA5':
        if area:
            grib_file = 'ERA5{}_{}_{}.grb'.format(area, date_list, hour)
        else:
            grib_file = 'ERA5_{}_{}.grb'.format(date_list, hour)

    elif model == 'ERAINT': grib_file = 'ERA-Int_{}_{}.grb'.format(date_list, hour)
    elif model == 'MERRA' : grib_file = 'merra-{}-{}.nc4'.format(date_list, hour)
    elif model == 'NARR'  : grib_file = 'narr-a_221_{}_{}00_000.grb'.format(date_list, hour)
    elif model == 'ERA'   : grib_file = 'ERA_{}_{}.grb'.format(date_list, hour)
    elif model == 'MERRA1': grib_file = 'merra-{}-{}.hdf'.format(date_list, hour)
    grib_files.append(os.path.join(grib_dir, grib_file))
    return grib_files

def snwe2str(snwe):
    """Get area extent in string"""
    if not snwe:
        return None
    s, n, w, e = snwe

    area = ''
    area += '_S{}'.format(abs(s)) if s < 0 else '_N{}'.format(abs(s))
    area += '_S{}'.format(abs(n)) if n < 0 else '_N{}'.format(abs(n))
    area += '_W{}'.format(abs(w)) if w < 0 else '_E{}'.format(abs(w))
    area += '_W{}'.format(abs(e)) if e < 0 else '_E{}'.format(abs(e))

    return area

def dload_grib_files(grib_files, tropo_model='ERA5', snwe=None):
    """Download weather re-analysis grib files using PyAPS
    Parameters: grib_files : list of string of grib files
    Returns:    grib_files : list of string
    """
    print('\n------------------------------------------------------------------------------')
    print('downloading weather model data using PyAPS ...')

    # Get date list to download (skip already downloaded files)
    grib_files_exist = check_exist_grib_file(grib_files, print_msg=True)
    grib_files2dload = sorted(list(set(grib_files) - set(grib_files_exist)))
    date_list2dload = [str(re.findall('\d{8}', os.path.basename(i))[0]) for i in grib_files2dload]
    print('number of grib files to download: %d' % len(date_list2dload))
    print('------------------------------------------------------------------------------\n')

    # Download grib file using PyAPS
    if len(date_list2dload) > 0:
        hour = re.findall('\d{8}[-_]\d{2}', os.path.basename(grib_files2dload[0]))[0].replace('-', '_').split('_')[1]
        grib_dir = os.path.dirname(grib_files2dload[0])

        # Check for non-empty account info in PyAPS config file
        check_pyaps_account_config(tropo_model)

        # try 3 times to download, then use whatever downloaded to calculate delay
        i = 0
        while i < 3:
            i += 1
            try:
                if tropo_model in ['ERA5', 'ERAINT']:
                    pa.ECMWFdload(date_list2dload, hour, grib_dir,
                                  model=tropo_model,
                                  snwe=snwe,
                                  flist=grib_files2dload)

                elif tropo_model == 'MERRA':
                    pa.MERRAdload(date_list2dload, hour, grib_dir)

                elif tropo_model == 'NARR':
                    pa.NARRdload(date_list2dload, hour, grib_dir)
            except:
                if i < 3:
                    print('WARNING: the {} attampt to download failed, retry it.\n'.format(i))
                else:
                    print('\n\n'+'*'*50)
                    print('WARNING: downloading failed for 3 times, stop trying and continue.')
                    print('*'*50+'\n\n')
                pass

    # check potentially corrupted files
    grib_files = check_exist_grib_file(grib_files, print_msg=False)
    return grib_files

def check_pyaps_account_config(tropo_model):
    """Check for input in PyAPS config file. If they are default values or are empty, then raise error.
    Parameters: tropo_model - str, tropo model being used to calculate tropospheric delay
    Returns:    None
    """
    # Convert MintPy tropo model name to data archive center name
    # NARR model included for completeness but no key required
    MODEL2ARCHIVE_NAME = {
        'ERA5' : 'CDS',
        'ERAI' : 'ECMWF',
        'MERRA': 'MERRA',
        'NARR' : 'NARR',
    }
    SECTION_OPTS = {
        'CDS'  : ['key'],
        'ECMWF': ['email', 'key'],
        'MERRA': ['user', 'password'],
    }

    # Default values in cfg file
    default_values = [
        'the-email-address-used-as-login@ecmwf-website.org',
        'the-user-name-used-as-login@earthdata.nasa.gov',
        'the-password-used-as-login@earthdata.nasa.gov',
        'the-email-adress-used-as-login@ucar-website.org',
        'your-uid:your-api-key',
    ]

    # account file for pyaps3 < and >= 0.3.0
    cfg_file = os.path.join(os.path.dirname(pa.__file__), 'model.cfg')
    rc_file = os.path.expanduser('~/.cdsapirc')

    # for ERA5: ~/.cdsapirc
    if tropo_model == 'ERA5' and os.path.isfile(rc_file):
        pass

    # check account info for the following models
    elif tropo_model in ['ERA5', 'ERAI', 'MERRA']:
        section = MODEL2ARCHIVE_NAME[tropo_model]

        # Read model.cfg file
        cfg_file = os.path.join(os.path.dirname(pa.__file__), 'model.cfg')
        cfg = ConfigParser()
        cfg.read(cfg_file)

        # check all required option values
        for opt in SECTION_OPTS[section]:
            val = cfg.get(section, opt)
            if not val or val in default_values:
                raise ValueError('PYAPS: No account info found for {} in {} section in file: {}'.format(tropo_model, section, cfg_file))

    return

def check_exist_grib_file(gfile_list, print_msg=True):
    """Check input list of grib files, and return the existing ones with right size."""
    gfile_exist = ut.get_file_list(gfile_list)
    if gfile_exist:
        file_sizes = [os.path.getsize(i) for i in gfile_exist] # if os.path.getsize(i) > 10e6]
        if file_sizes:
            comm_size = ut.most_common([i for i in file_sizes])
            if print_msg:
                print('common file size: {} bytes'.format(comm_size))
                print('number of grib files existed    : {}'.format(len(gfile_exist)))

            gfile_corrupt = []
            for gfile in gfile_exist:
                if os.path.getsize(gfile) < comm_size * 0.9:
                    gfile_corrupt.append(gfile)
        else:
            gfile_corrupt = gfile_exist

        if gfile_corrupt:
            if print_msg:
                print('------------------------------------------------------------------------------')
                print('corrupted grib files detected! Delete them and re-download...')
                print('number of grib files corrupted  : {}'.format(len(gfile_corrupt)))

            for gfile in gfile_corrupt:
                print('remove {}'.format(gfile))
                os.remove(gfile)
                gfile_exist.remove(gfile)

            if print_msg:
                print('------------------------------------------------------------------------------')
    return gfile_exist

def get_delay(grib_file, tropo_model, delay_type, dem, inc, lat, lon, mask=None, verbose=False):
    """Get delay matrix using PyAPS for one acquisition
    Parameters: grib_file       - str, grib file path
                tropo_model     - str, GAM model
                delay_type      - str, dry/wet/comb
                dem/inc/lat/lon - 2D np.ndarray in float32 for DEM, incidence angle, latitude/longitude
                verbose         - bool, verbose message
    Returns:    pha             - 2D np.ndarray in float32, single path tropospheric delay
                                  temporally absolute, spatially referenced to ref_y/x
    """
    if verbose:
        print('GRIB FILE: {}'.format(grib_file))

    # initiate pyaps object
    aps_obj = pa.PyAPS(grib_file,
                       grib=tropo_model,
                       Del=delay_type,
                       dem=dem,
                       inc=inc,
                       lat=lat,
                       lon=lon,
                       mask=mask,
                       verb=verbose)

    # estimate delay
    pha = np.zeros((aps_obj.ny, aps_obj.nx), dtype=np.float32)
    aps_obj.getdelay(pha)

    # reverse the sign for consistency between different phase correction steps/methods
    pha *= -1
    return pha

