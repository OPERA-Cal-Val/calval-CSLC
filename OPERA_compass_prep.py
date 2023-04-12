#!/usr/bin/env python

'''Wrapper to prep OPERA beta products for ALE/RLE analysis'''

import os
import glob
from pathlib import Path
import argparse


def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description= \
            'Wrapper to prep OPERA beta products for ALE/RLE analysis')
    parser.add_argument('-i', '--input-dir', dest='input_dir', type=str,
        required=True, help='Path to OPERA beta products')
    parser.add_argument('-p', '--path-id', dest='path_id', type=str,
        required=True,
        help='Specify a single valid path ID')
    parser.add_argument('-b', '--burst-id', dest='burst_id', type=str,
        default=None,
        help='List of burst IDs, with each entry separated by commas. '
             'If None, all available burst IDs are prepped. '
             'e.g. t094_200136_iw3 for path ID 94, id 20016, and subswith 3 '
             '(default: None)')
    parser.add_argument('-w', '--working-dir', dest='work_dir',
        default='stack',
        help='Directory to store intermediate and final results '
             '(default: stack)')

    return parser.parse_args(args=iargs)

def run_prep(inps):
    '''Prep OPERA beta products'''

    # get all OPERA product dirs
    product_dirs = glob.glob(inps.input_dir + '/*')

    # filter by path ID
    if len(inps.path_id) > 3:
        raise Exception(f'Specified path id {inps.path_id} invalid. '
                        f'Check input')
    if len(inps.path_id) < 3:
        # prepend 0s if necessary
        num_0s = 3 - len(inps.path_id)
        inps.path_id = num_0s * '0' + inps.path_id
    inps.path_id = 'T' + inps.path_id
    product_dirs = [i for i in product_dirs if inps.path_id in Path(i).name]

    # filter by specified burst ID(s)
    if inps.burst_id is None:
        inps.burst_id = []
        for i in product_dirs:
            p_name = Path(i).name.split('_')
            burst_id_str = p_name[4].lower().replace('-','_')
            if burst_id_str not in inps.burst_id:
                inps.burst_id.append(burst_id_str)
    else:
        inps.burst_id = list(set(inps.burst_id.split(',')))
        filt_product_dirs = []
        for i in product_dirs:
            p_name = Path(i).name.split('_')
            burst_id_str = p_name[4].lower().replace('-','_')
            if burst_id_str in inps.burst_id:
                filt_product_dirs.append(i)
        product_dirs = filt_product_dirs

    if product_dirs == []:
        raise Exception('No scenes available with specified inputs')

    # create work directory(ies) and populate with symlinks
    for i in inps.burst_id:
        burst_path = Path(inps.work_dir, i)
        # create path
        Path(burst_path).mkdir(parents=True, exist_ok=True)

        # find products with matching burst IDs
        burst_id_str = i.upper().replace('_','-')
        filt_product_dirs = [j for j in product_dirs if burst_id_str in j]
        for j in filt_product_dirs:
            p_name = Path(j).name.split('_')
            date_name = p_name[6][:8]
            prod_path = Path(burst_path, date_name)
            # create path
            Path(prod_path).mkdir(parents=True, exist_ok=True)
            # create symlink with OPERA beta h5 file
            prod_rename = i + '_' + date_name + '.h5'
            prod_rename = Path(prod_path, prod_rename)
            h5_file = glob.glob(j + '/*.h5')[0]
            Path(prod_rename).symlink_to(h5_file)
        

if __name__ == "__main__":
    '''Run ALE/RLE prep workflow from command line'''
    # load arguments from command line
    inps = createParser()

    # Run workflow
    run_prep(inps)
