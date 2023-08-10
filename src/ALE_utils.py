import numpy as np
import math
import scipy
import isce3

'''
Collection utility functions to find the corner reflectors based on intensity peak
'''

def oversample_slc(slc,sampling=1,y=None,x=None):
    if y is None:
        y = np.arange(slc.shape[0])
    if x is None:
        x = np.arange(slc.shape[1])

    [rows, cols] = np.shape(slc)
    
    slcovs = isce3.signal.point_target_info.oversample(slc,sampling)

    y_orign_step = y[1]-y[0]
    y_ovs_step = y_orign_step/sampling
    x_orign_step = x[1]-x[0]
    x_ovs_step = x_orign_step/sampling

    y = np.arange(y[0],y[-1]+y_orign_step,y_ovs_step)
    x = np.arange(x[0],x[-1]+x_orign_step,x_ovs_step)

    return slcovs,y,x

def findCR(data,y,x,x_bound=[-np.inf,np.inf],y_bound=[-np.inf,np.inf],method="sinc"):
    '''
    Find the location of CR with fitting
    '''
    max_ind = np.argmax(data)
    max_data = data[max_ind]
    
    def _sinc2D(x,x0,y0,a,b,c):
        return c*np.sinc(a*(x[0]-x0))*np.sinc(b*(x[1]-y0))
    
    def _para2D(x,x0,y0,a,b,c,d):
        return a*(x[0]-x0)**2+b*(x[1]-y0)**2+c*(x[0]-x0)*(x[1]-y0)+d

    if method == "sinc":
        # using sinc function for fitting 
        xdata = np.vstack((x,y))
        p0 = [x[max_ind],y[max_ind],0.7,0.7,max_data]
        bounds = ([x_bound[0],y_bound[0],0,0,0],[x_bound[1],y_bound[1],1,1,np.inf])
        popt = scipy.optimize.curve_fit(_sinc2D,xdata,data,p0=p0,bounds=bounds)[0]
        xloc = popt[0]; yloc = popt[1]
    elif method == "para":
        #using paraboloid function for fitting
        xdata = np.vstack((x,y))
        p0 = [x[max_ind],y[max_ind],-1,-1,1,1]
        bounds = ([x_bound[0],y_bound[0],-np.inf,-np.inf,-np.inf,0],[x_bound[1],y_bound[1],0,0,np.inf,np.inf])
        popt = scipy.optimize.curve_fit(_para2D,xdata,data,p0=p0,bounds=bounds)[0]
        xloc = popt[0]; yloc = popt[1]

    return yloc,xloc

def interpolate_correction_layers(xcoor, ycoor, data, method):
    '''
    Interpolate the correction layers
        xcoor: UTM coordinates of the CSLC data along x-axis, shape: (N,)
        ycoor: UTM coordinates of the CSLC data along y-axis, shape: (N,)
        data: Correction layer values
        method: Interpolation method ('nearest', 'linear', 'cubic')
    '''
    
    Xcslc, Ycslc = np.meshgrid(xcoor, ycoor)
    xx = np.linspace(xcoor.min(),xcoor.max(),num=data.shape[1])
    yy = np.linspace(ycoor.min(),ycoor.max(),num=data.shape[0])
    Xdata, Ydata = np.meshgrid(xx,yy)

    # Interpolate
    points = list(zip(Xdata.ravel(), Ydata.ravel()))
    values = data.ravel()

    data_resampl = scipy.interpolate.griddata(points, values, (Xcslc, Ycslc), method=method)

    return np.flipud(data_resampl)

def en2rdr(E, N, az_angle, inc_angle):
    rng = E * np.sin(np.deg2rad(inc_angle)) * np.cos(np.deg2rad(az_angle - 90)) * -1 + N * np.sin(np.deg2rad(inc_angle)) * np.sin(np.deg2rad(az_angle - 90)) 
    grng = rng / np.sin((np.deg2rad(inc_angle)))
    azi = E * np.sin(np.deg2rad(az_angle - 90)) * -1 + N * np.cos(np.deg2rad(az_angle - 90))

    return grng, azi

def get_snr_peak(img: np.ndarray, cutoff_percentile: float=3.0):
    '''
    Estimate the signal-to-noise ration (SNR) of the peak
    in the input image patch
    Parameter
    ---------
    img: numpy.ndarray
        SLC image patch to calculate the SNR
    cutout: float
        Cutout ratio of high and low part of the signal to cutoff
    Returns
    -------
    snr_peak_db: float
        SNR of the peak in decibel (db)
    '''

    power_arr = img.real ** 2 + img.imag ** 2

    # build up the mask array
    thres_low = np.nanpercentile(power_arr, cutoff_percentile)
    thres_high = np.nanpercentile(power_arr, 100 - cutoff_percentile)
    mask_threshold = np.logical_and(power_arr < thres_low,
                                    power_arr > thres_high)
    mask_invalid_pixel = np.logical_and(power_arr <= 0.0,
                                        np.isnan(power_arr))
    ma_power_arr = np.ma.masked_array(power_arr,
                                      mask=np.logical_and(mask_threshold,
                                                          mask_invalid_pixel))

    peak_power = np.nanmax(power_arr)
    mean_background_power = np.mean(ma_power_arr)

    snr_peak_db = np.log10(peak_power / mean_background_power) * 10.0

    return snr_peak_db
