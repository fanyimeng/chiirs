import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
import pandas as pd
from tqdm import tqdm
from astropy.io import fits


def ell_collide(ell1, ell2):
    def ell_rad(ell_para, angle):
        '''
        input ell_para[a, b, theta]
        output r from a spot (costheta*r, sintheta*r) on the ell.
        '''
        a = ell_para[2]
        b = ell_para[3]
        cos_angle = np.cos(np.radians(angle))
        sin_angle = np.sin(np.radians(angle))
        r = b * a / np.sqrt(cos_angle**2 * b ** 2 + sin_angle**2 * a**2)
        x = ell_para[0] + np.cos(np.radians(angle + ell_para[4])) * r
        y = ell_para[1] + np.sin(np.radians(angle + ell_para[4])) * r
        return (x, y)

    def ell_in(ell_para, x, y):
        g_ell_center = ell_para[:2]
        g_ell_width = ell_para[2]
        g_ell_height = ell_para[3]
        angle = ell_para[4]

        cos_angle = np.cos(np.radians(180. - angle))
        sin_angle = np.sin(np.radians(180. - angle))

        xc = x - g_ell_center[0]
        yc = y - g_ell_center[1]

        xct = xc * cos_angle - yc * sin_angle
        yct = xc * sin_angle + yc * cos_angle

        rad_cc = (xct**2 / (g_ell_width)**2) + \
            (yct**2 / (g_ell_height)**2)
        return rad_cc

    # collide_degree = 0
    ang_list = np.arange(0, 360, 1)
    x, y = ell_rad(ell2, ang_list)
    ell_in_list = ell_in(ell1, x, y)
    # print(ell_in_list)
    collide_degree = (ell_in_list < 1).sum()
    # print(collide_degree)
    # for i in range(36):
    #     x, y = ell_rad(ell2, i * 10.)
    #     if ell_in(ell1, x, y) < 1:
    #         collide_degree = collide_degree + 1
    return collide_degree


def collidelist(coresa, coresb, outputcsv='cores_006_coll.csv'):
    collist = []
    for i in tqdm(range(len(coresa['cenx'].to_numpy()))):
        colldeg = 0
        coll_obj = -1
        for j in range(len(coresb['cenx'].to_numpy())):
            newcolldeg = ell_collide(coresa.iloc[i].to_numpy()[:5],
                                     coresb.iloc[j].to_numpy()[:5])
            if newcolldeg > colldeg:
                # print(newcolldeg)
                coll_obj = j
        collist.append(coll_obj)
    collist = np.array(collist)
    coresa['collide'] = collist
    coresa.to_csv(outputcsv)
    print((collist >= 0).sum(), '/', collist.shape)


hdu096 = fits.open('../data/rg_096.fits')[0]
hdu006 = fits.open('../data/abcd_006_lpf_new.fits')[0]
data096 = hdu096.data
data006 = hdu006.data
data096[data096 < 0] = 0
data006[data006 < 0] = 0
header096 = hdu096.header
header006 = hdu006.header
flxf096 = 4e3 * np.log(2.0) * header096['CDELT2'] * header096['CDELT1'] / (
    3.14159 * header096['BMAJ'] * header096['BMIN'])
flxf006 = 4e3 * np.log(2.0) * header006['CDELT2'] * header006['CDELT1'] / (
    3.14159 * header006['BMAJ'] * header006['BMIN'])

rmshdu096 = fits.open('../data/rg_096_rms_100_1e-3_cubic.fits')[0]
rmshdu006 = fits.open('../data/abcd_006_lpf_new.fits')[0]
rms096 = rmshdu096.data
rms006 = rmshdu006.data

x1d = np.arange(0, data096.shape[0])
y1d = np.arange(0, data096.shape[1])
x2d = x1d[np.newaxis, :]
y2d = y1d[:, np.newaxis]

cenf_006 = []
cenf_096 = []
cene_006 = []
cene_096 = []
r = 0.62 * 1 / 0.05

df = pd.read_csv('./cores_096_coll.csv')

for i in tqdm(range(len(df['cenx']))):
    mask = (x2d - int(df['cenx'][i]))**2 + (y2d - int(df['ceny'][i]))**2 < r**2
    cenf_006.append(np.nansum(data006[mask]))
    cenf_096.append(np.nansum(data096[mask]))
    cene_006.append(np.nansum(rms006[mask]))
    cene_096.append(np.nansum(rms096[mask]))
cenf_006 = np.array(cenf_006) * abs(flxf006)
cenf_096 = np.array(cenf_096) * abs(flxf096)
cene_006 = np.array(cene_006) * abs(flxf006)
cene_096 = np.array(cene_096) * abs(flxf096)
df['cenf_006'] = cenf_006
df['cenf_096'] = cenf_096
df['cene_006'] = cene_006
df['cene_096'] = cene_096
df.to_csv('./cores_096_flxcompare.csv')
