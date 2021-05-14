from astropy.table import Table
import numpy as np
import os
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import coordinates
from astropy import wcs
from astropy.io import fits


def getRADEC(header, cen_coor):
    w = wcs.WCS(header).celestial
    x = cen_coor[0]
    y = cen_coor[1]
    ra, dec = w.all_pix2world(x, y, 0)
    c = coordinates.SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    return (c.ra.hms[:]), (c.dec.dms[0], np.abs(c.dec.dms[1]), np.abs(c.dec.dms[2]))


def vErrFormatter(value1, err):
    dist = -np.int(np.log10(np.abs(err))) + 1
    formatstr = '$%.*f \\pm %.*f$'
    if dist > 0:
        line = formatstr % (dist, value1, dist, err)
    else:
        line = formatstr % (0, value1, 0, err)
    return line


# def vErrFormatter2(value1, err1, err2):
#     err1 = np.abs(err1)
#     err2 = np.abs(err2)
#     dist = -np.int(np.log10(np.abs(np.nanmax([err1, err2])))) + 1
#     # dist = np.int(np.nanmax([-np.int(np.log10(np.abs(err2))) + 1, dist]))
#     formatstr = '$%.*f^{+%.*f}_{-%.*f}$'
#     if dist > 0:
#         line = formatstr % (dist, value1, dist, err1, dist, err2)
#     else:
#         line = formatstr % (0, value1, 0, err1, 0, err2)
#     return line

def vErrFormatter2(value1, err1, err2):
    formatstr = '$%.*f^{+%.*f}_{-%.*f}$'
    dist = 2
    if dist > 0:
        line = formatstr % (dist, value1, dist, err1, dist, err2)
    return line


df = pd.read_csv('006.csv')


f = open('../obspara.tex', "w")
header = fits.open('../data/img_006.fits')[0].header
for i in range(df.shape[0]):
    row = df.iloc[i, :]
    num_str = '%i' % (row['006_idx'])

    ra_str = '%07.4f' % (getRADEC(header,
                                 [row['x_img'], row['y_img']])[0][-1])
    dec_str = '%02i:%05.2f' % (getRADEC(header,
                                      [row['x_img'], row['y_img']])[1][1],
                             getRADEC(header,
                                      [row['x_img'], row['y_img']])[1][2])

    shape_str = '%.2f' % (row['r_img'] * 0.05)
    s006_str = vErrFormatter(row['flx_006'], row['err_006'])
    s022_str = '...'
    s096_str = '...'
    separator = '&'
    outstr = separator.join([num_str,
                             ra_str,
                             dec_str,
                             shape_str,
                             s006_str,
                             s022_str,
                             s096_str
                             ])
    outstr = outstr + '\\\\ \n'
    print(outstr)
    f.write(outstr)
f.close()
os.system('pdflatex ../chiirs.tex')

# print(row['ALPHAPEAK_J2000'])


'''
X_IMAGE
Y_IMAGE
FLUX_RADIUS
FLUX_ISO
FLUXERR_ISO
FLAGS
A_IMAGE
B_IMAGE
THETA_IMAGE
XPEAK_IMAGE
YPEAK_IMAGE
ALPHAPEAK_J2000
DELTAPEAK_J2000
FLUX_PEAK
FLUXERR_PEAK
FLX_PEAK3SIGMA_006
FLX_PEAK3SIGMA_010
ERR_PEAK3SIGMA_006
ERR_PEAK3SIGMA_010
FLX_PEAK1PIX_006
FLX_PEAK1PIX_010
ERR_PEAK1PIX_006
ERR_PEAK1PIX_010
EM_6_065
EM_6_065_UP
EM_6_065_DOWN
EM_6
EM_6_UP
EM_6_DOWN
EM_610
EM_610_UP
EM_610_DOWN
'''

# for c in t:
#     num = num + 1
#     nstr = '%3i & ' % (num)
#     corelongid = '%1.0f%04.0f-%04.0f & ' % (c['dec_m'] - 20,
#                                             c['dec_s'] * 1e2,
#                                             c['ra_s'] * 1e2)
#     ra = '%02.0f:%02.0f:%07.4f & ' % (c['ra_h'],
#                                       c['ra_m'],
#                                       c['ra_s'])
#     dec = '%+02.0f:%02.0f:%05.2f & ' % (c['dec_d'],
#                                         c['dec_m'],
#                                         c['dec_s'])

#     coretype = c['band'].replace('096', 'II').replace('both', 'I')

#     flux096 = vErrFormatter(c['096f'], c['096p'], c['096r'])

#     agname = '&' + c['ag_name'].replace('core_', '').split('_')[0]
#     # if c['ag_rapix'] > 0:
#     #     shiftx = c(['ag_rapix'] - c['ximg'])*0.05
#     #     shifty = c(['ag_rapix'] - c['ximg'])*0.05
#     if c['ag_peakflux'] > -1e8:
#         agflux = '&%.2f' % c['ag_peakflux']
#     else:
#         agflux = '&--'
#     corecluster = '&' + c['cluster']
#     coreasso = '&' + c['asso']
#     outflow = '&' + outflowdic[c['outflow']]
#     line = ''.join([nstr, ra, dec, coretype, flux096,
#                     agname, agflux, corecluster, outflow, coreasso, ' \\\\ \n'])
#     print(line)
#     f.write(line)
# f.close()

# os.system('pdflatex ../chiirs.tex')
