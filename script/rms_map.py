from astropy.io import fits
import numpy as np
from scipy import ndimage
from sklearn.preprocessing import Imputer


def rms_map(inputfits, outputfits, masklimit):

    def rmsfunc(x):
        return np.std(x)

    im = fits.open(inputfits)[0].data.squeeze()
    im_header = fits.open(inputfits)[0].header
    im[im > masklimit] = np.nan
    print(inputfits)

    res = ndimage.generic_filter(im, rmsfunc, size=(100, 100))
    fits.writeto(outputfits, data=res,
                 header=im_header, overwrite=True)


filelist = ['../data/rg_022.4.fits']
limitlist = [4e-3]
for i in range(len(filelist)):
    rms_map(filelist[i],
            filelist[i].replace('.fits', '_rms_100_5e-4.fits'),
            limitlist[i])
