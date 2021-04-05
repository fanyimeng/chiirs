from astropy.io import fits


hdu = fits.open('SGRB2_1.3CM.fits')[0]
header = hdu.header
header['BMAJ'] = 0.27/3600
header['BMAJ'] = 0.23/3600
header['BPA'] = 70

fits.writeto('022.4.fits', data = hdu.data, header=header)