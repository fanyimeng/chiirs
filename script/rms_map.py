from astropy.io import fits
import numpy as np
from scipy import ndimage
from sklearn.preprocessing import Imputer


def rms_map(inputfits, outputfits, masklimit):
    # imp = Imputer(strategy="mean")

    def rmsfunc(x):
        return np.std(x)

    im = fits.open(inputfits)[0].data.squeeze()
    im_header = fits.open(inputfits)[0].header
    im[im > masklimit] = np.nan
    print(inputfits)

    res = ndimage.generic_filter(im, rmsfunc, size=(100, 100))
    # res = imp.fit_transform(res.transpose()).transpose()
    fits.writeto(outputfits, data=res,
                 header=im_header, overwrite=True)


rootdir = '/Users/meng/Dropbox/astro/densecore_data/'
# filelist = ['rg_006.fits',
#             'rg_010.fits',
#             'rg_022.fits',
#             'rg_033.fits',
#             'rg_046.fits',
#             'rg_096.fits',
#             'rg_274.fits',
#             'rg_343.fits']
filelist = ['abcd_006_lpf.fits']
limitlist = [5e-4,
             7e-4,
             5e-4,
             5e-4,
             5e-4,
             5e-4,
             4e-2,
             6e-2]
for i in range(len(filelist)):
    rms_map(rootdir + filelist[i],
            filelist[i].replace('.fits', '_rms_100.fits'),
            limitlist[i])
