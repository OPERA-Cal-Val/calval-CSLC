#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import datetime as dt
from src.ALE_utils import en2rdr

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='summarizing the outcomes of pycuampcor')
    parser.add_argument("--txtfile", dest='txt',
            required=True,type=str, help='two column text file (burst_id ref_date)')
    parser.add_argument("--dirCSV", dest='dirCSV',
            required=True,type=str, help='directory of RLE ts results')
    parser.add_argument("--dir_en2rdr", dest='dir_en2rdr',
            required=True,type=str, help='directory of en2rdr files (incidence and azimuth angle (unit: deg))')
    parser.add_argument("--summaryCSV", dest='summaryCSV',
            default='summary_RLE.csv',type=str, help='summary csv file (summary_RLE.csv)')                        
    return parser.parse_args(args=iargs)

def if_pass(ts,requirement):
    bool_pass = (ts > -requirement) & (ts < requirement)
    pass_rate = np.count_nonzero(bool_pass)/len(bool_pass)
    return bool_pass, pass_rate, pass_rate>0.8

def run(inps):

    txtfile = inps.txt 
    df_txt = pd.read_csv(txtfile, sep='\s+', names=['burstId','refDate'], header=None)

    dirCSV = inps.dirCSV
    dir_en2rdr = inps.dir_en2rdr
    summaryCSV = inps.summaryCSV

    summary_df = pd.DataFrame()
    # add placeholder columns
    summary_df['BurstID'] = None
    summary_df['RLE'] = None
    summary_df['ReferenceDate'] = None
    summary_df['dGroundRange'] = None
    summary_df['dAzimuth'] = None

    _bid = []; _rle = []; _refD = []; _dGrng = []; _dAz = []

    for index, row in df_txt.iterrows():

        _id = row.burstId
        _id = _id.lower()
        _ref = row.refDate

        _csv = f'{dirCSV}/RLE_ts_{_id.upper()}.csv'
        _en2rdr = f'{dir_en2rdr}/en2rdr_{_id}.csv'

        f = open(_en2rdr)
        inc_angle, az_angle =f.read().split(' ')
        inc_angle = float(inc_angle); az_angle = float(az_angle)

        _df = pd.read_csv(_csv)
        _df['date'] = pd.to_datetime(_df['date'], format='%Y%m%d') 

        grng, azi = en2rdr(_df['rg_avg'],_df['az_avg'], az_angle, inc_angle) 
        grng_std, azi_std = en2rdr(_df['rg_std'],_df['az_std'], az_angle, inc_angle)

        _df['grng_avg'] = grng
        _df['azi_avg'] = azi
        _df['grng_std'] = np.abs(grng_std)
        _df['azi_std'] = np.abs(azi_std)
        # Get closest date to user-input reference date
        refDate_ = pd.to_datetime(_ref, format='%Y%m%d')
        refDate = _df.loc[(_df['date']-refDate_).abs().idxmin(),'date']

        # Reference the stack to the reference date
        _df['grng_avg'] = _df['grng_avg'] - (_df[_df.date==refDate]['grng_avg'].values[0])
        _df['azi_avg'] = _df['azi_avg'] - (_df[_df.date==refDate]['azi_avg'].values[0])

        _grng_avg = np.round(_df['grng_avg'].mean(),3)
        _grng_std = np.round(_df['grng_avg'].std(),3)
        _azi_avg = np.round(_df['azi_avg'].mean(),3)
        _azi_std = np.round(_df['azi_avg'].std(),3)

        rg_bool_pass, rg_pass_rate, rg_pass_or_not = if_pass(_df['grng_avg'],0.5)
        az_bool_pass, az_pass_rate, az_pass_or_not = if_pass(_df['azi_avg'],0.75)

        if rg_pass_or_not and az_pass_or_not:
            _dec = 'PASS'
        else:
            _dec = 'FAIL'

        _bid.append(_id.upper())
        _rle.append(_dec)
        _refD.append(dt.datetime.strftime(refDate,'%Y%m%d'))
        _dGrng.append(f'{_grng_avg} +/- {_grng_std}')
        _dAz.append(f'{_azi_avg} +/- {_azi_std}')

    summary_df['BurstID'] = _bid
    summary_df['RLE'] = _rle
    summary_df['ReferenceDate'] = _refD
    summary_df['dGroundRange'] = _dGrng
    summary_df['dAzimuth'] =_dAz

    summary_df.to_csv(summaryCSV, index=False)

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()
    
    # Run workflow
    run(inps)
