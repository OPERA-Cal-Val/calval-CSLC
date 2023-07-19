#!/usr/bin/env python3
import argparse
import numpy as np
from osgeo import gdal,osr
import os
import time
from PyCuAmpcor import PyCuAmpcor
import h5py
from src.RLE_utils import array2raster

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='running pycuampcor on a pair of slcs to obtain offset')
    parser.add_argument("--slc_dir", dest='slc_dir',
                         required=True, type=str, help='slc directory')
    parser.add_argument("--dateref", dest="dateref",
                         required=True, type=str, help='reference date (YYYYMMDD)')
    parser.add_argument("--datesec", dest="datesec",
                         required=True, type=str, help='secondary date (YYYYMMDD)')    
    parser.add_argument("--out_dir", dest="out_dir",
                         required=True, type=str, help='output directory')
    parser.add_argument("--deviceID", dest="deviceID",
                         default=0, type=int, help='device ID (default: 0)')
    parser.add_argument("--ww", dest="ww",
                         default=64, type=int, help='window size width (default: 64)')
    parser.add_argument("--wh", dest="wh",
                         default=64, type=int, help='window size height (default: 64)')
    parser.add_argument("--nwdc", dest="nwdc",
                         default=20, type=int, help='number of window down in chunk for GPU processing (default: 20)')
    parser.add_argument("--nwac", dest="nwac",
                         default=20, type=int, help='number of window accross in chunk for GPU processing (default: 20)')
    return parser.parse_args(args=iargs)

def run(inps):
    #input parameters
    slc_dir = inps.slc_dir
    date_ref = inps.dateref
    date_sec = inps.datesec
    out_dir = inps.out_dir
    deviceID = inps.deviceID
    windowSizeWidth = inps.ww
    windowSizeHeight = inps.wh
    numberWindowDownInChunk = inps.nwdc
    numberWindowAcrossInChunk = inps.nwac
    
    os.makedirs(out_dir,exist_ok=True)

    ref_file = slc_dir +'/' + date_ref+'.slc'
    refSLCvrt = ref_file+'.vrt'
    sec_file = slc_dir +'/' + date_sec+'.slc'
    secSLCvrt = sec_file+'.vrt'
    
    ds = gdal.Open(ref_file)
    width, height = ds.RasterXSize, ds.RasterYSize
    geotransform = ds.GetGeoTransform()
    x_start = geotransform[0]
    y_start = geotransform[3]
    dx = geotransform[1]
    dy = geotransform[5]
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    epsg = proj.GetAttrValue('AUTHORITY',1)

    #pycuampcor parameters
    halfSearchRangeDown = 20
    halfSearchRangeAcross = 20

    skipSampleDown = 16
    skipSampleAcross = 16

    referenceStartPixelDownStatic = halfSearchRangeDown
    referenceStartPixelAcrossStatic = halfSearchRangeAcross

    numberWindowAcross = int(np.floor((width 
                                      - referenceStartPixelAcrossStatic - windowSizeWidth 
                                      - halfSearchRangeAcross)/ skipSampleAcross/10))*10

    numberWindowDown = int(np.floor((height
                                     - referenceStartPixelDownStatic - windowSizeHeight 
                                     - halfSearchRangeDown)/skipSampleDown/10))*10

    assert referenceStartPixelAcrossStatic+windowSizeWidth+halfSearchRangeAcross+numberWindowAcross*skipSampleAcross < width
    assert referenceStartPixelDownStatic+windowSizeHeight+halfSearchRangeDown+numberWindowDown*skipSampleDown < height

    ref_start_center_across_pix = referenceStartPixelAcrossStatic+(windowSizeWidth-1)/2
    ref_start_center_down_pix = referenceStartPixelDownStatic+(windowSizeHeight-1)/2
    ref_start_center_across_coor = x_start+ref_start_center_across_pix*dx
    ref_start_center_down_coor = y_start+ref_start_center_down_pix*dy
    offset_xcoor = np.arange(0,numberWindowAcross)*skipSampleAcross*dx+ref_start_center_across_coor
    offset_ycoor = np.arange(0,numberWindowDown)*skipSampleDown*dy+ref_start_center_down_coor

    print("Run the Pycuampcor...", end=" ")
    st_time = time.time()

    objOffset = PyCuAmpcor() # create the processor
    objOffset.algorithm = 0 # cross-correlation method 0=freq 1=time

    objOffset.corrStatWindowSize = 7
    objOffset.corrSurfaceOverSamplingMethod = 0
    objOffset.rawDataOversamplingFactor = 2

    objOffset.nStreams = 2  # cudaStreams; multiple streams to overlap data transfer with gpu calculations
    objOffset.referenceImageName = refSLCvrt
    objOffset.referenceImageHeight = height # RasterYSize
    objOffset.referenceImageWidth = width # RasterXSize
    objOffset.secondaryImageName = secSLCvrt
    objOffset.secondaryImageHeight = height
    objOffset.secondaryImageWidth = width
    objOffset.windowSizeWidth = windowSizeWidth # template window size
    objOffset.windowSizeHeight = windowSizeHeight
    objOffset.halfSearchRangeDown = halfSearchRangeDown # search range
    objOffset.halfSearchRangeAcross = halfSearchRangeAcross
    objOffset.derampMethod = 0   #0=mag for TOPS, 1=deramping, else=skip deramping

    objOffset.skipSampleDown = skipSampleDown # strides between windows
    objOffset.skipSampleAcross = skipSampleAcross
    objOffset.numberWindowDownInChunk = numberWindowDownInChunk 
    objOffset.numberWindowAcrossInChunk = numberWindowAcrossInChunk
    objOffset.corrSurfaceOverSamplingFactor = 32 # oversampling factor for correlation surface
    objOffset.corrSurfaceZoomInWindow = 16  # area in correlation surface to be oversampled

    objOffset.useMmap = 1 # default using memory map as buffer, if having troubles, set to 0
    objOffset.mmapSize = 32 # mmap or buffer size used for transferring data from file to gpu, in GB

    objOffset.numberWindowDown = numberWindowDown # number of windows to be processed
    objOffset.numberWindowAcross = numberWindowAcross

    covFile = out_dir+f'/{date_ref}_{date_sec}.cov'
    snrFile = out_dir+f'/{date_ref}_{date_sec}.snr'    #BIP
    offsetFile = out_dir+f'/{date_ref}_{date_sec}.off'    #BIP
    grossoffFile = out_dir+f'/{date_ref}_{date_sec}.grossoff'

    objOffset.covImageName = covFile
    objOffset.snrImageName = snrFile
    objOffset.offsetImageName = offsetFile
    objOffset.grossOffsetImageName = grossoffFile
    objOffset.deviceID = deviceID

    objOffset.setupParams()
    objOffset.referenceStartPixelDownStatic = halfSearchRangeDown # starting pixel offset
    objOffset.referenceStartPixelAcrossStatic = halfSearchRangeDown
    objOffset.setConstantGrossOffset(0, 0) # gross offset between reference and secondary images
    objOffset.checkPixelInImageRange() # check whether there is something wrong with

    objOffset.runAmpcor()

    end_time = time.time()
    time_taken = end_time - st_time

    print(f'Done, {time_taken}s taken for {date_ref} - {date_sec} pair \n')

    snrF = open(snrFile,"rb")
    snr = np.fromfile(snrF, dtype='<f4', count=numberWindowAcross*numberWindowDown)
    snr = np.reshape(snr,(-1,numberWindowAcross))

    offsetF = open(offsetFile,"rb")
    offset = np.fromfile(offsetF, dtype='<f4', count=2*numberWindowAcross*numberWindowDown)
    offset = np.reshape(offset,(-1,numberWindowAcross*2))

    rg_offset = offset[:,1::2]
    az_offset = offset[:,::2]

    rg_offset[np.isnan(snr)] = np.nan
    az_offset[np.isnan(snr)] = np.nan

    offset_dx = offset_xcoor[1]-offset_xcoor[0]
    offset_dy = offset_ycoor[1]-offset_ycoor[0]
    extent = (offset_xcoor[0]-offset_dx/2,offset_xcoor[-1]+offset_dx/2, offset_ycoor[-1]+offset_dy/2, offset_ycoor[0]-offset_dy/2)

    rgoffsetFile = out_dir+f'/{date_ref}_{date_sec}.rg_off.tif'
    azoffsetFile = out_dir+f'/{date_ref}_{date_sec}.az_off.tif'
    snrFile = out_dir+f'/{date_ref}_{date_sec}.snr.tif'

    if os.path.exists(rgoffsetFile): os.remove(rgoffsetFile)
    if os.path.exists(azoffsetFile): os.remove(azoffsetFile)
    if os.path.exists(snrFile): os.remove(snrFile)
    
    #writing offset and snr files in UTM coordinates
    array2raster(rgoffsetFile,extent[0],extent[3],offset_dx,offset_dy,int(epsg),rg_offset*dx)
    array2raster(azoffsetFile,extent[0],extent[3],offset_dx,offset_dy,int(epsg),az_offset*np.abs(dy))
    array2raster(snrFile,extent[0],extent[3],offset_dx,offset_dy,int(epsg),snr)

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()

    # Run workflow
    run(inps)
