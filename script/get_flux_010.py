import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy.io import fits


hdu006 = fits.open('../data/abcd_006_lpf_new.fits')[0]
hdu010 = fits.open('../data/abcd_010_lpf_new.fits')[0]
data006 = hdu006.data
data010 = hdu010.data
data006[data006 < 0] = 0
data010[data010 < 0] = 0
header006 = hdu006.header
header010 = hdu010.header
flxf006 = 4e3 * np.log(2.0) * header006['CDELT2'] * header006['CDELT1'] / (
    3.14159 * header006['BMAJ'] * header006['BMIN'])
flxf010 = 4e3 * np.log(2.0) * header010['CDELT2'] * header010['CDELT1'] / (
    3.14159 * header010['BMAJ'] * header010['BMIN'])
print(flxf006, flxf010)
rmshdu006 = fits.open('../data/006_sex_test_check_rms.fits')[0]
rmshdu010 = fits.open('../data/010_sex_test_check_rms.fits')[0]
rms006 = rmshdu006.data
rms010 = rmshdu010.data

x1d = np.arange(0, data006.shape[0])
y1d = np.arange(0, data006.shape[1])
x2d = x1d[np.newaxis, :]
y2d = y1d[:, np.newaxis]

cenf_010 = []
cenf_006 = []
cene_010 = []
cene_006 = []
r = 0.62 * 1 / 0.05

df = pd.read_csv('./cores006.csv')

for i in tqdm(range(len(df['cenx']))):
    # mask = (x2d - int(df['cenx'][i]))**2 + (y2d - int(df['ceny'][i]))**2 < r**2

    mask = ((((x2d - int(df['cenx'][i])) * np.cos(np.radians(df['ellt'][i])) -
              (y2d - int(df['ceny'][i])) * np.sin(np.radians(df['ellt'][i]))) /
             df['ella'][i])**2 +
            (((x2d - int(df['cenx'][i])) * np.sin(np.radians(df['ellt'][i])) -
              (y2d - int(df['ceny'][i])) * np.cos(np.radians(df['ellt'][i]))) /
             df['ellb'][i])**2
            ) < 1
    '''
            xct = xc * cos_angle - yc * sin_angle
        yct = xc * sin_angle + yc * cos_angle
    '''
    # cenf_010.append(np.nanmax(data010[mask]))
    cenf_010.append(data010[data006[mask].argmax()])
    cenf_006.append(np.nanmax(data006[mask]))
    # print(data010[mask].argmax(), rms010[mask].shape)
    cene_010.append(rms010[data010[mask].argmax()])
    cene_006.append(rms006[data006[mask].argmax()])
cenf_010 = np.array(cenf_010)  # * abs(flxf010)
cenf_006 = np.array(cenf_006)  # * abs(flxf006)
cene_010 = np.array(cene_010)  # * abs(flxf010)
cene_006 = np.array(cene_006)  # * abs(flxf006)
df['ELLPEAK_FLX_010'] = cenf_010
df['ELLPEAK_FLX_006'] = cenf_006
df['ELLPEAK_ERR_010'] = cene_010
df['ELLPEAK_ERR_006'] = cene_006
df.to_csv('./cores006.csv')
