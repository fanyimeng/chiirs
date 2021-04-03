import numpy as np
import os

radlist = np.arange(0, 2.0, 0.01)
for rad in radlist:
    file = 'photomodel_%.2f.fits' % rad
    imsmooth(imagename=file,
             major='0.65arcsec',
             minor='0.65arcsec',
             pa='0deg',
             targetres=True,
             outfile='sm_' + file.replace('fits', 'im'), overwrite=True)
    exportfits(imagename='sm_' + file.replace('fits', 'im'),
               fitsimage='sm_' + file,
               dropdeg=True, overwrite=True)
    os.system('rm -rf ' + 'sm_' + file.replace('fits', 'im'))
