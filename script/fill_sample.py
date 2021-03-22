import numpy as np
from astropy.io import fits
from scipy import interpolate
import matplotlib.pyplot as plt


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]


def ffill(arr):
    mask = np.isnan(arr)
    idx = np.where(~mask, np.arange(mask.shape[1]), 0)
    np.maximum.accumulate(idx, axis=1, out=idx)
    out = arr[np.arange(idx.shape[0])[:, None], idx]
    return out

# My modification to do a backward-fill


def bfill(arr):
    mask = np.isnan(arr)
    idx = np.where(~mask, np.arange(mask.shape[1]), mask.shape[1] - 1)
    idx = np.minimum.accumulate(idx[:, ::-1], axis=1)[:, ::-1]
    out = arr[np.arange(idx.shape[0])[:, None], idx]
    return out


def intfill(arr):
    grid_x, grid_y = np.mgrid[0:arr.shape[0]:1, 0:arr.shape[1]:1]
    out = interpolate.griddata((range(arr.shape[0]), range(
        arr.shape[1])), arr, (grid_x, grid_y), method='nearest')
    return out


def sample(arr, samplerate=100):
    newarr = arr[:int(arr.shape[0] / samplerate),
                 :int(arr.shape[1] / samplerate)].copy()
    for i in range(newarr.shape[0] - 1):
        for j in range(newarr.shape[1] - 1):
            i_start = i * samplerate
            i_end = i * samplerate + samplerate
            j_start = j * samplerate
            j_end = j * samplerate + samplerate
            newarr[i, j] = np.max(arr[i_start:i_end, j_start:j_end])
    return newarr


def unsample(arr, newarr, samplerate=100):
    newnewarr = arr.copy()
    nanlist = np.argwhere(np.isnan(newnewarr))
    for nancoor in nanlist:
        newnewarr[nancoor[0], nancoor[1]] = newarr[(int(nancoor[0]/samplerate)),(int(nancoor[1]/samplerate))]
    return newnewarr

def interfill(arr):
    #from https://stackoverflow.com/questions/37662180/interpolate-missing-values-2d-python
    x = np.arange(0, arr.shape[1])
    y = np.arange(0, arr.shape[0])
    #mask invalid values
    arr = np.ma.masked_invalid(arr)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~arr.mask]
    y1 = yy[~arr.mask]
    newarr = arr[~arr.mask]
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),
                              (xx, yy),
                                 method='cubic')
    return GD1


def fillnan(rms_file, out_file):
    hdu = fits.open(rms_file)[0]
    data = hdu.data
    newdata = interfill(sample(data))
    newnewdata = unsample(data, newdata)
    # plt.imshow(newnewdata)
    # plt.show()
    # print(newdata)
    fits.writeto(out_file, data=newnewdata,
                 header=hdu.header, overwrite=True)


# fillnan('rg_006_rms_100.fits', 'rg_006_rms_100_int.fits')
# fillnan('rg_010_rms_100.fits', 'rg_010_rms_100_int.fits')
# fillnan('rg_096_rms_100.fits', 'rg_096_rms_100_int.fits')
# fillnan('abcd_006_lpf_rms_100.fits', 'abcd_006_lpf_rms_100_int.fits')
# fillnan('rg_096_rms_100_1e-3.fits', "rg_096_rms_100_1e-3_cubic.fits")
fillnan("dustonly_rms_100_1e-3.fits","dustonly_rms_100_cubic.fits")
