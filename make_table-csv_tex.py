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
    err = abs(err)
    dist = -np.int(np.log10(np.abs(err))) + 1
    if ('%.*f' % (dist, err))[-1] == '0' and ('%.*f' % (dist, err))[-2] != '.':
        dist = dist - 1
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


df = pd.read_csv('script/006.csv')
df022 = pd.read_csv('script/022.csv')
df096 = pd.read_csv('script/096.csv')

header = fits.open('./data/img_006.fits')[0].header


f = open('./obspara.tex', "w")
for i in range(df.shape[0]):
    row = df.iloc[i, :]
    num_str = '%i' % (row['006_idx'])

    ra_str = '%07.4f' % (getRADEC(header,
                                  [row['x_img'], row['y_img']])[0][-1])
    dec_str = '%02i:%05.2f' % (getRADEC(header,
                                        [row['x_img'], row['y_img']])[1][1],
                               getRADEC(header,
                                        [row['x_img'], row['y_img']])[1][2])

    r006_str = '%.2f' % (row['r_img'] * 0.05)
    s006_str = vErrFormatter(row['flx_006'],
                             np.nanmax([row['err_006'], 0.1]))

    if row['w022'] < 0.5:
        s022_str = '...'
        r022_str = '...'
    else:
        row022 = df022.iloc[int(row['022_idx']), :]
        s022_str = vErrFormatter(row022['flx_022'],
                                 np.nanmax([row022['err_022'], 0.7]))
        r022_str = '%.2f' % (row022['r_img'] * 0.05)

    if row['w096'] < 0.5:
        s096_str = '...'
        r096_str = '...'
    else:
        # print(row['096_idx'])
        row096 = df096.iloc[int(row['096_idx']), :]
        s096_str = vErrFormatter(row096['flx_096'],
                                 np.nanmax([row096['err_096'], 0.2]))
        r096_str = '%.2f' % (row096['r_img'] * 0.05)
    separator = '&'
    outstr = separator.join([num_str,
                             ra_str,
                             dec_str,
                             r006_str,
                             s006_str,
                             r022_str,
                             s022_str,
                             r096_str,
                             s096_str
                             ])
    outstr = outstr + '\\\\ \n'
    print(outstr)
    f.write(outstr)
f.close()


f = open('./derivedpara.tex', "w")
'''
    \\# &
    $r_{\rm calc}$ & 
    ${\rm EM}_{\rm calc}$ &
    $\\log_{10)(\\dot{N}_{\rm Ly}/{\rm s^{-1}})$ &
    $n_{\rm e^{-}}$ &
    $n_{\rm H_2}$ &
    $r_i$ &
    $t$\\


        &
    $\times 10^-2$ pc & 
    $\times 10^7 {\rm pc cm^{-6}}$ &
    &
    $\times 10^4 {\rm cm^{-3}}$&
    $\times 10^8 {\rm cm^{-3}}$&
    pc &
    $\times 10^4 {\rm yrs}$\\

'''
for i in range(df.shape[0]):
    row = df.iloc[i, :]
    num_str = '%i' % (row['006_idx'])
    r_calc = '%.2f' % (row['R'] / 25. / 1e-3)
    em_calc = '%.2f' % (row['EM'] / 1e7)
    nly = '%.2f' % (row['NLY'])
    n_e = '%.2f' % (np.sqrt(row['EM'] / (row['R'] / 25.)) / 1e4)
    n_d = '%.2f' % (row['gasn'] / 1e6)
    r_i = '%.3f' % (row['r_i'] / 1e-3)
    t = '%.1f' % (row['t'] / 1e4)
    if row['w022'] < 0.5:
        # n_e = '%.2f' % (np.sqrt(row['EM'] / (row['R'] / 25.)) / 1e4)
        r_calc = '$^*$' + r_calc
        em_calc = '$^*$' + em_calc
        nly = '$^*$' + nly
        n_e = '$^*$' + n_e
        # n_d = '$^*$' + n_d
        r_i = '$^*$' + r_i
        t = '$^*$' + t
    if row['dustf_096'] < 0:
        n_d = '...'
        r_i = '$^{\\dagger}$' + r_i 
        t = '$^{\\dagger}$' + t 

    separator = '&'
    outstr = separator.join([num_str,
                             r_calc,
                             em_calc,
                             nly,
                             n_e,
                             n_d,
                             r_i,
                             t
                             ])
    outstr = outstr + '\\\\ \n'
    print(outstr)
    f.write(outstr)
f.close()


os.system('pdflatex chiirs.tex')

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
