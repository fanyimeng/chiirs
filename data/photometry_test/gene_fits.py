from astropy.io import fits
import numpy as np

def circPixList(center, rad):
    pixlist = []
    x1 = int(center[0]) - int(rad) - 1
    x2 = int(center[0]) + int(rad) + 1
    y1 = int(center[1]) - int(rad) - 1
    y2 = int(center[1]) + int(rad) + 1
    for i in range(x1, x2):
        for j in range(y1, y2):
            if (i-center[0])**2 + (j-center[1])**2 <= rad**2:
                pixlist.append([i, j])
    if rad < 1:
        pixlist = [center]
    return pixlist


hdu = fits.open(
    '/Users/meng/Dropbox/astro/densecore_data/ccsgrb2/formattedData/rg_096.fits')[0]
data = hdu.data
header = hdu.header

data = np.zeros_like(data[:2000, :2000]).copy()
header['NAXIS1'] = 2000
header['NAXIS2'] = 2000
header['CDELT1'] = -1.388888888889E-06
header['CRPIX1'] = 1.000000000000E+03
header['CDELT2'] = 1.388888888889E-06
header['CRPIX2'] = 1.000000000000E+03
header['BMAJ'] = 1.388888888889E-06
header['BMIN'] = 1.388888888889E-06

radlist = np.arange(0, 2.0, 0.01)
for rad in radlist:
    outdata = data.copy()
    rad_inpix = rad / 3600 / header['CDELT2']
    print(rad_inpix)
    for pix in circPixList([1000, 1000], rad_inpix):
        outdata[pix[0], pix[1]] = 1
    fits.writeto('photomodel_%.2f.fits' % rad, outdata, header, overwrite=True)
    print('photomodel_%.1f.fits' % rad)
