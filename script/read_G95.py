import pandas as pd
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates
import numpy as np
from astropy.coordinates import FK5
from astropy.time import Time


df = pd.read_csv("gaume1995_format.tab", header=0, delim_whitespace=True)
# df = pd.read_csv('gaume1995.csv')

header = fits.open('../data/abcd_006_lpf_new.fits')[0].header
mywcs = wcs.WCS(header).celestial
df['G95_X'] = np.ones_like(df['G95_RA'])
df['G95_Y'] = np.ones_like(df['G95_DEC'])
df['G95_A'] = df['G95_A']/0.05
df['G95_B'] = df['G95_B']/0.05
for i in range(df.shape[0]):
    cen_coor = coordinates.SkyCoord(str(df['G95_RA'][i]),
                                    str(df['G95_DEC'][i]),
                                    unit=(u.h, u.deg),
                                    frame='fk4',
                                    equinox='B1950')
    cen_coor = cen_coor.transform_to(FK5(equinox='J2000'))
    # print(cen_coor)
    (df['G95_X'][i], df['G95_Y'][i]) =  mywcs.wcs_world2pix([[cen_coor.ra.deg,
                                                             cen_coor.dec.deg]], 0)[0].astype(int)
    # print('%s %s' % (df['G95_RA'][i], df['G95_DEC'][i])
    print('%.0f %.0f' % (df['G95_X'][i], df['G95_Y'][i]))


df.to_csv('gaume1995.csv', index=False)
