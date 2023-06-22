from osgeo import gdal,osr
import h5py
import os
import datetime as dt
import numpy as np
from mintpy.utils import readfile
from mintpy.cli import ifgram_inversion, load_data
import pandas as pd
import fsspec

def hdf_read(input_hdf):

    grid_path = f'data'
    metadata_path = f'metadata'
    burstmetadata_path = f'metadata/processing_information/input_burst_metadata'
    id_path = f'identification'
    pol = 'VV'

    with h5py.File(input_hdf,'r') as h5:
        xcoor = h5[f'{grid_path}/x_coordinates'][:]
        ycoor = h5[f'{grid_path}/y_coordinates'][:]
        dx = h5[f'{grid_path}/x_spacing'][()].astype(int)
        dy = h5[f'{grid_path}/y_spacing'][()].astype(int)
        epsg = h5[f'{grid_path}/projection'][()].astype(int)
        slc = h5[f'{grid_path}/{pol}'][:]
        sensing_start = h5[f'{burstmetadata_path}/sensing_start'][()].astype(str)
        sensing_stop = h5[f'{burstmetadata_path}/sensing_stop'][()].astype(str)
        dims = h5[f'{burstmetadata_path}/shape'][:]
        bounding_polygon =h5[f'{id_path}/bounding_polygon'][()].astype(str) 
        orbit_direction = h5[f'{id_path}/orbit_pass_direction'][()].astype(str)
        center_lon, center_lat = h5[f'{burstmetadata_path}/center'] 
        date = dt.datetime.strptime(sensing_start.astype(str),'%Y-%m-%d %H:%M:%S.%f').strftime('%Y%m%d')

        #deramping and adding flattening phase back
        azimuth_carrier_phase = h5[f'{grid_path}/azimuth_carrier_phase'][:]
        flatten_phase = h5[f'{grid_path}/flattening_phase'][:]        
        ramp = np.exp(1j*azimuth_carrier_phase)
        flat_phase = np.exp(1j*flatten_phase)
        slc = slc*np.conj(ramp)*np.conj(flat_phase)
        
    return xcoor, ycoor, dx, dy, epsg, slc, date

def hdf_stream(path_h5):

    s3f = fsspec.open(path_h5, mode='rb', anon=True, default_fill_cache=False)
    print(f'streaming: {path_h5}')

    grid_path = f'data'
    metadata_path = f'metadata'
    burstmetadata_path = f'metadata/processing_information/input_burst_metadata'
    id_path = f'identification'
    pol = 'VV'

    with h5py.File(s3f.open(),'r') as h5:
        xcoor = h5[f'{grid_path}/x_coordinates'][:]
        ycoor = h5[f'{grid_path}/y_coordinates'][:]
        dx = h5[f'{grid_path}/x_spacing'][()].astype(int)
        dy = h5[f'{grid_path}/y_spacing'][()].astype(int)
        epsg = h5[f'{grid_path}/projection'][()].astype(int)
        slc = h5[f'{grid_path}/{pol}'][:]
        sensing_start = h5[f'{burstmetadata_path}/sensing_start'][()].astype(str)
        sensing_stop = h5[f'{burstmetadata_path}/sensing_stop'][()].astype(str)
        dims = h5[f'{burstmetadata_path}/shape'][:]
        bounding_polygon =h5[f'{id_path}/bounding_polygon'][()].astype(str)
        orbit_direction = h5[f'{id_path}/orbit_pass_direction'][()].astype(str)
        center_lon, center_lat = h5[f'{burstmetadata_path}/center']
        date = dt.datetime.strptime(sensing_start.astype(str),'%Y-%m-%d %H:%M:%S.%f').strftime('%Y%m%d')

        #deramping and adding flattening phase back
        azimuth_carrier_phase = h5[f'{grid_path}/azimuth_carrier_phase'][:]
        flatten_phase = h5[f'{grid_path}/flattening_phase'][:]
        ramp = np.exp(1j*azimuth_carrier_phase)
        flat_phase = np.exp(1j*flatten_phase)
        slc = slc*np.conj(ramp)*np.conj(flat_phase)

    return xcoor, ycoor, dx, dy, epsg, slc, date 
    
def convert_to_slcvrt(xcoor, ycoor, dx, dy, epsg, slc, date, outdir):

     os.makedirs(outdir,exist_ok=True)

     height, width = slc.shape

     slc_file = outdir + '/' + date+'.slc'
     slc_vrt = slc_file+'.vrt'

     outtype = '<f'  #little endian (float)
     dtype = gdal.GDT_CFloat32
     drvout = gdal.GetDriverByName('ENVI')
     raster_out = drvout.Create(slc_file, width,height, 1, dtype)
     raster_out.SetGeoTransform([xcoor[0],dx,0.0,ycoor[0],0.0,dy])

     srs = osr.SpatialReference()
     srs.ImportFromEPSG(int(epsg))
     raster_out.SetProjection(srs.ExportToWkt())

     band_out = raster_out.GetRasterBand(1)
     band_out.WriteArray(slc)
     band_out.FlushCache()
     del band_out

     command = 'gdal_translate ' + slc_file + ' ' + slc_vrt + f' > {outdir}/tmp.LOG'
     os.system(command)

def array2raster(outrasterfile,OriginX, OriginY, pixelWidth,pixelHeight,epsg,array):
    #generating geotiff file from 2D array

    cols = array.shape[1]
    rows = array.shape[0]
    originX = OriginX
    originY = OriginY

    driver = gdal.GetDriverByName('ENVI')
    outRaster = driver.Create(outrasterfile, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def simple_SBAS_stats(offlist,snrlist,out_dir,snr_thr):
    #offlist: offset filelist
    #snrlist: snr filelist
    #out_dir: output directory
    #snr_thr: snr threshold

    num_pairs = offlist.shape[0]

    refd = []
    secd = []

    for _ in offlist:
        refd.append(_[0:8])
        secd.append(_[9:17])

    days = refd + secd
    days = (np.unique(sorted(days))).tolist()
    n_days = len(days)      #number of unique days

    #building a design matrix
    D = np.zeros((num_pairs,n_days))   #initialization

    for ii in range(num_pairs):
        ref_index = days.index(refd[ii])
        sec_index = days.index(secd[ii])
        D[ii,ref_index] = -1
        D[ii,sec_index] = 1

    invD = np.linalg.pinv(D)  #inverse of a Design matrix

    #opening first tiff file for obtaining parameters
    _ = out_dir + '/' + offlist[0] 
    ds = gdal.Open(_, gdal.GA_ReadOnly)
    _ = ds.GetRasterBand(1).ReadAsArray()
    row, col = _.shape
    transform = ds.GetGeoTransform()
    minX = transform[0]
    maxY = transform[3]
    spacingX = transform[1]
    spacingY = transform[5]
    maxX = minX + (col-1)*spacingX
    minY = maxY + (row-1)*spacingY
    proj = ds.GetProjection()
    _ = None

    off_3d = np.zeros((row,col,num_pairs))

    for ii, (offF,snrF) in enumerate(zip(offlist, snrlist)):
        offFile = out_dir + '/' + offF
        snrFile = out_dir + '/' + snrF

        ds = gdal.Open(snrFile, gdal.GA_ReadOnly)
        snr = ds.GetRasterBand(1).ReadAsArray()

        ds = gdal.Open(offFile, gdal.GA_ReadOnly)
        off = ds.GetRasterBand(1).ReadAsArray()
        off[snr<snr_thr] = np.nan
        off_3d[:,:,ii] = off

    #SBAS inversion for time-series estimates
    ts_off = np.einsum('ijk,lk->ijl', off_3d, invD)
    norm_res_off = np.sqrt(np.sum((off_3d - np.einsum('ijk,lk->ijl', ts_off, D))**2,axis = 2))  #L2 norm residual

    #removing pixels with a large L2 norm residual
    normResThr = np.nanmin(norm_res_off) + (np.nanmax(norm_res_off) 
                                              - np.nanmin(norm_res_off))*0.75
                                              #threshold of L2 norm residual #np.nanquantile(norm_res_off, 0.75)   
    indRes = (norm_res_off>normResThr)
    off_3d[indRes,:] = np.nan

    first_ts_off = ts_off[:,:,0]
    ts_off_all = dict()

    for ii, ID in enumerate(days):
        dat = ts_off[:,:,ii] - first_ts_off     #first data becomes zero
        ts_off_all[ID] = dat    #time-series range offset

    _avg = []
    _std = []

    for day in days:
        _avg.append(np.nanmean(ts_off_all[day]))
        _std.append(np.nanstd(ts_off_all[day]))

    return _avg, _std, days

def mintpy_SBAS_stats(rgofflist,azofflist,snrlist,out_dir,snr_thr,q=0.25):
    #rgofflist: list of range offset files
    #azofflist: list of azimuth offset files
    #snrlist: list of snr files
    #out_dir: output directory
    #snr_thr: snr threshold
    #q: quantile threshold for excluding outliers

    refd = []
    secd = []

    for _ in rgofflist:
        refd.append(_[0:8])
        secd.append(_[9:17])

    metafile = 'metadata.txt'

    #creating metadata file
    with open(metafile, 'w') as f:
        for tif_, ref_, sec_ in zip(rgofflist, refd, secd):
            f.write(f"{tif_} {ref_} {sec_}\n")
    
    with open(metafile, 'a') as f:
        for tif_, ref_, sec_ in zip(azofflist, refd, secd):
            f.write(f"{tif_} {ref_} {sec_}\n")

    with open(metafile, 'a') as f:
        for tif_, ref_, sec_ in zip(snrlist, refd, secd):
            f.write(f"{tif_} {ref_} {sec_}\n")

    #input configure file for Mintpy
    cfgfile = 'smallbaselineApp.cfg'

    script = f'''
    ##------------------------ smallbaselineApp.cfg ------------------------##
    ########## computing resource configuration
    mintpy.compute.maxMemory = 10 #[float > 0.0], auto for 4, max memory to allocate in GB
    mintpy.compute.cluster   = local #[local / slurm / pbs / lsf / none], auto for none, cluster type
    mintpy.compute.numWorker = 10 #[int > 1 / all], auto for 4 (local) or 40 (non-local), num of workers
    mintpy.compute.config    = auto #[none / slurm / pbs / lsf ], auto for none (same as cluster), config name
    ########## 1. load_data
    mintpy.load.processor      = cosicorr  #[isce, aria, hyp3, gmtsar, snap, gamma, roipac], auto for isce
    mintpy.load.autoPath       = auto  #[yes / no], auto for no, use pre-defined auto path
    mintpy.load.updateMode     = no  #[yes / no], auto for yes, skip re-loading if HDF5 files are complete
    mintpy.load.compression    = auto  #[gzip / lzf / no], auto for no.

    mintpy.load.metaFile       = ./{metafile}

    ##---------interferogram datasets:
    mintpy.load.unwFile        = auto  #[path pattern of unwrapped interferogram files]
    mintpy.load.corFile        = auto  #[path pattern of spatial coherence       files]
    mintpy.load.connCompFile   = auto  #[path pattern of connected components    files], optional but recommended
    mintpy.load.intFile        = auto  #[path pattern of wrapped interferogram   files], optional
    mintpy.load.ionoFile       = auto  #[path pattern of ionospheric delay       files], optional
    mintpy.load.magFile        = auto  #[path pattern of interferogram magnitude files], optional
    ##---------offset datasets (optional):
    mintpy.load.azOffFile      = ./{out_dir}/*az_off.tif  #[path pattern of azimuth offset file], optional
    mintpy.load.rgOffFile      = ./{out_dir}/*rg_off.tif  #[path pattern of range   offset file], optional
    mintpy.load.offSnrFile     = ./{out_dir}/*snr.tif  #[path pattern of offset signal-to-noise ratio file], optional

    ##---------geometry datasets:
    mintpy.load.demFile        = auto  #[path of DEM file]
    mintpy.load.lookupYFile    = auto  #[path of latitude /row   /y coordinate file], not required for geocoded data
    mintpy.load.lookupXFile    = auto  #[path of longitude/column/x coordinate file], not required for geocoded data
    mintpy.load.incAngleFile   = auto  #[path of incidence angle file], optional but recommended
    mintpy.load.azAngleFile    = auto  #[path of azimuth   angle file], optional
    mintpy.load.shadowMaskFile = auto  #[path of shadow mask file], optional but recommended
    mintpy.load.waterMaskFile  = auto  #[path of water  mask file], optional but recommended
    mintpy.load.bperpFile      = auto  #[path pattern of 2D perpendicular baseline file], optional
    ##---------multilook (optional):
    ## multilook while loading data with nearest interpolation, to reduce dataset size
    mintpy.load.ystep          = auto    #[int >= 1], auto for 1 - no multilooking
    mintpy.load.xstep          = auto    #[int >= 1], auto for 1 - no multilooking
    ##---------subset (optional):
    ## if both yx and lalo are specified, use lalo option unless a) no lookup file AND b) dataset is in radar coord
    mintpy.subset.yx           = auto    #[y0:y1,x0:x1 / no], auto for no
    mintpy.subset.lalo         = auto    #[S:N,W:E / no], auto for no
    '''

    with open(cfgfile,"w") as f:
        f.writelines(script)

    #prep data for MintPy
    load_data.main('-t smallbaselineApp.cfg'.split())

    #time-series inversion with MintPy
    tsRgFile = 'timeseriesRg.h5'  #time series hn5 file in range
    tsAzFile = 'timeseriesAz.h5'  #time series h5 file in azimuth

    cmd = f'inputs/offsetStack.h5 -i rangeOffset -w no --min-norm-phase --md offsetSNR --mt {snr_thr} -o {tsRgFile} residualInvRg.h5 numInvOffsetRg.h5'
    ifgram_inversion.main(cmd.split())

    cmd = f'inputs/offsetStack.h5 -i azimuthOffset -w no --min-norm-phase --md offsetSNR --mt {snr_thr} -o {tsAzFile} residualInvAz.h5 numInvOffsetAz.h5'
    ifgram_inversion.main(cmd.split())

    ts_rg_data, meta = readfile.read(tsRgFile)
    ts_az_data, meta = readfile.read(tsAzFile)
    
    tslist= readfile.get_slice_list(tsAzFile)
    days = [ ii[11:] for ii in tslist ]

    rg_avg = []; rg_std = []
    az_avg = []; az_std = []

    for ii, day in enumerate(days):
        
        _ = ts_rg_data[ii,:]
        if (ii != 0): _[_==0.] = np.nan    
        _thr1 = np.nanquantile(_,q)  #removing outliers
        _thr2 = np.nanquantile(_,1-q)
        _[_<_thr1] = np.nan
        _[_>_thr2] = np.nan

        rg_avg.append(np.nanmean(_))
        rg_std.append(np.nanstd(_))
        
        _ = ts_az_data[ii,:]
        if (ii != 0): _[_==0.] = np.nan 
        _thr1 = np.nanquantile(_,q)
        _thr2 = np.nanquantile(_,1-q)
        _[_<_thr1] = np.nan
        _[_>_thr2] = np.nan
        
        az_avg.append(np.nanmean(_))
        az_std.append(np.nanstd(_)) 

    return rg_avg, rg_std, az_avg, az_std
