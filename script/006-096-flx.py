from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy import coordinates
from scipy import stats
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from astropy.table import QTable, Table, Column
from matplotlib.ticker import *
import seaborn as sns

plt.rcParams['figure.dpi'] = 75.
plt.rcParams['savefig.dpi'] = 75.
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.clf()
fig = plt.figure(figsize=(6.0, 6.0))
ax = plt.subplot(1, 1, 1)

datatable = pd.read_csv('cores_096_flxcompare.csv')

yarr = datatable['cenf_096'].copy()
xarr= datatable['cenf_006'].copy()
ratio = yarr/xarr
xlimarr = xarr.copy()
ylimarr = yarr.copy()
xlimarr[xarr<0.1] = 0.1
ylimarr[yarr<0.1] = 0.1
xarr[xarr<0.1] = 0.1
yarr[yarr<0.1] = 0.1

ax.set_xscale("log", nonposx='clip')
ax.xaxis.set_major_formatter(LogFormatter())
ax.xaxis.set_minor_formatter(LogFormatter())
ax.set_yscale("log", nonposy='clip')
ax.yaxis.set_major_formatter(LogFormatter())
ax.yaxis.set_minor_formatter(LogFormatter())
# for axis in [ax.xaxis, ax.yaxis]:
#     axis.set_major_formatter(ScalarFormatter())


# ax.scatter(datatable['ag_peakflux'], datatable['096f'])


print([np.nanmin(xarr) * 0.1, np.nanmax(xarr) * 10])
print([np.nanmin(yarr) * 0.1, np.nanmax(yarr) * 10])
ax.set_xlim([np.nanmin(xlimarr) * 0.1, np.nanmax(xlimarr) * 10])
ax.set_ylim([np.nanmin(ylimarr) * 0.1, np.nanmax(ylimarr) * 10])


def plotscatter(ax, x, y, subsize, color, datatable):
    zx1 = x  - subsize
    zx2 = x  + subsize
    zy1 = y  - subsize
    zy2 = y  + subsize
    yarr = datatable['cenf_096'].copy()
    xarr= datatable['cenf_006'].copy()
    xarr[xarr<0.1] = 0.1
    yarr[yarr<0.1] = 0.1
    yarr[datatable['cenx']<zx1] = np.nan
    xarr[datatable['cenx']<zx1] = np.nan
    yarr[datatable['cenx']>zx2] = np.nan
    xarr[datatable['cenx']>zx2] = np.nan
    yarr[datatable['ceny']<zy1] = np.nan
    xarr[datatable['ceny']<zy1] = np.nan
    yarr[datatable['ceny']>zy2] = np.nan
    xarr[datatable['ceny']>zy2] = np.nan
    ax.errorbar(xarr, yarr,
                yerr=0, 
                xerr=0,
                ls='none',
                marker='o', markersize='0.5', lw=0.1,
                ecolor=color,
                color=color,
                capsize=4,
                zorder = 10000)
    return 0 


plotscatter(ax, 4780, 5780, 350, (0.9,0.3,0.3,0.5), datatable)
plotscatter(ax, 4700, 4860, 350, (0.3,0.9,0.3,0.5), datatable)
plotscatter(ax, 4700, 3900, 350, (0.3,0.3,0.9,0.5), datatable)
plotscatter(ax, 4380, 2050, 750, (0.5,0.5,0.5,0.5), datatable)

ax.set_aspect(np.log(ax.get_xlim()[
              1] / ax.get_xlim()[0]) / np.log(ax.get_ylim()[1] / ax.get_ylim()[0]))

ax.set_ylabel('$S_\\mathrm{cen}$(3 mm) [mJy]')
ax.set_xlabel('$S_\\mathrm{cen}$(6 cm) [mJy]')

ax.tick_params(axis='both', which='major', direction='in', size=10,
               bottom=True, top=True, left=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', size=6,
               bottom=True, top=True, left=True, right=True)

xx = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], 1)
yy = xx*(96/6.)**-0.1
ax.plot(xx, yy, '-', lw=0.6, color='gray')
yy = xx*(96/6.)**0.3
ax.plot(xx, yy, '-', lw=0.6, color='gray')
yy = xx*(96/6.)**2
ax.plot(xx, yy, '-', lw=0.6, color='gray')

ax.text(100, 1, 'N',   color = (0.9,0.3,0.3,0.5))
ax.text(100, 1/2, 'M', color = (0.3,0.9,0.3,0.5))
ax.text(100, 1/4, 'S', color = (0.3,0.3,0.9,0.5))
ax.text(100, 1/8, 'DS',color = (0.5,0.5,0.5,0.5))



'''
xx = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], 1)

te = 1e4


emcolor = (1, 0.3, 0.3, 0.9)
em = 1e6
xx = np.arange(ax.get_xlim()[0], 40, 1)
yy = (xx / 1e3 / 2 * 50) ** 2 * em / (5.36e3 * 6**0.1 * te**0.35)
print(np.max(yy), np.min(yy))
ax.plot(xx, yy, '-', lw=0.6, color=emcolor)
ax.text(40, np.max(yy), '$10^6$ pc cm$^{-6}$',
        rotation=40,
        color=emcolor,
        fontsize=10)
xx = np.arange(80, ax.get_xlim()[1], 1)
yy = (xx / 1e3 / 2 * 50) ** 2 * em / (5.36e3 * 6**0.1 * te**0.35)
ax.plot(xx, yy, '-', lw=0.6, color=emcolor)

em = 1e7
xx = np.arange(ax.get_xlim()[0], 40, 1)
yy = (xx / 1e3 / 2 * 50) ** 2 * em / (5.36e3 * 6**0.1 * te**0.35)
print(np.max(yy), np.min(yy))
ax.plot(xx, yy, '-', lw=0.6, color=emcolor)
ax.text(40, np.max(yy), '$10^7$ pc cm$^{-6}$',
        rotation=40,
        color=emcolor,
        fontsize=10)
xx = np.arange(80, ax.get_xlim()[1], 1)
yy = (xx / 1e3 / 2 * 50) ** 2 * em / (5.36e3 * 6**0.1 * te**0.35)
ax.plot(xx, yy, '-', lw=0.6, color=emcolor)


em = 1e8
xx = np.arange(ax.get_xlim()[0], 26, 1)
yy = (xx / 1e3 / 2 * 50) ** 2 * em / (5.36e3 * 6**0.1 * te**0.35)
print(np.max(yy), np.min(yy))
ax.plot(xx, yy, '-', lw=0.6, color=emcolor)
ax.text(26, np.max(yy), '$10^8$ pc cm$^{-6}$',
        rotation=40,
        color=emcolor,
        fontsize=10)
xx = np.arange(50, ax.get_xlim()[1], 1)
yy = (xx / 1e3 / 2 * 50) ** 2 * em / (5.36e3 * 6**0.1 * te**0.35)
ax.plot(xx, yy, '-', lw=0.6, color=emcolor)


ncolor = (0.5, 0.5, 1, 0.8)

n = 1e4
xx = np.arange(ax.get_xlim()[0], 40, 1)
yy = (xx / 1e3 / 2 * 50) ** 2 / (5.36e3 * 6 **
                                  0.1 * te**0.35) * n**2 * np.sqrt(xx / 1e3)
print(np.max(yy), np.min(yy))
ax.plot(xx, yy, '--', lw=0.6, color=ncolor)
ax.text(40, np.max(yy), '$10^4$ cm$^{-3}$',
        rotation=50,
        color=ncolor,
        fontsize=8)
xx = np.arange(58, ax.get_xlim()[1], 1)
yy = (xx / 1e3 / 2 * 50) ** 2 / (5.36e3 * 6 **
                                  0.1 * te**0.35) * n**2 * np.sqrt(xx / 1e3)
ax.plot(xx, yy, '--', lw=0.6, color=ncolor)

n = 1e3
xx = np.arange(ax.get_xlim()[0], 40, 1)
yy = (xx / 1e3 / 2 * 50) ** 2 / (5.36e3 * 6 **
                                  0.1 * te**0.35) * n**2 * np.sqrt(xx / 1e3)
print(np.max(yy), np.min(yy))
ax.plot(xx, yy, '--', lw=0.6, color=ncolor)
ax.text(40, np.max(yy), '$10^3$ cm$^{-3}$',
        rotation=50,
        color=ncolor,
        fontsize=8)
xx = np.arange(58, ax.get_xlim()[1], 1)
yy = (xx / 1e3 / 2 * 50) ** 2 / (5.36e3 * 6 **
                                  0.1 * te**0.35) * n**2 * np.sqrt(xx / 1e3)
ax.plot(xx, yy, '--', lw=0.6, color=ncolor)


taucolor = (0.5, 0.5, 0.5, 0.8)

tau = 0.01
xx = np.arange(5, ax.get_xlim()[1], 1)
yy = tau * 6.**2 * (xx / 1e3 / 2 * 50) ** 2 * 1e4
print(np.max(yy), np.min(yy))
ax.plot(xx, yy, '--', lw=0.6, color=taucolor)
ax.text(5, np.min(yy), '$\\tau = 10^{-2}$',
        rotation=40,
        ha = 'right',
        va = 'top',
        color=taucolor,
        fontsize=8)
xx = np.arange(ax.get_xlim()[0], 3.5, 1)
yy = tau * 6.**2 * (xx / 1e3 / 2 * 50) ** 2 * 1e4
ax.plot(xx, yy, '--', lw=0.6, color=taucolor)

tau = 0.001
xx = np.arange(5, ax.get_xlim()[1], 1)
yy = tau * 6.**2 * (xx / 1e3 / 2 * 50) ** 2 * 1e4
print(np.max(yy), np.min(yy))
ax.plot(xx, yy, '--', lw=0.6, color=taucolor)
ax.text(5, np.min(yy), '$\\tau = 10^{-3}$',
        rotation=40,
        ha = 'right',
        va = 'top',
        color=taucolor,
        fontsize=8)
xx = np.arange(ax.get_xlim()[0], 3.5, 1)
yy = tau * 6.**2 * (xx / 1e3 / 2 * 50) ** 2 * 1e4
ax.plot(xx, yy, '--', lw=0.6, color=taucolor)

'''



# tau = 0.01
# yy3 = tau * 6.**2 * (xx / 1e3 / 2 * 50) ** 2 * 1e4
# print(np.max(yy3), np.min(yy3))
# ax.plot(xx, yy3, '--', lw=0.6, color=(0.5, 0.5, 0.5, 0.8))

# tau = 0.001
# yy3 = tau * 6.**2 * (xx / 1e3 / 2 * 50) ** 2 * 1e4
# print(np.max(yy3), np.min(yy3))
# ax.plot(xx, yy3, '--', lw=0.6, color=(0.5, 0.5, 0.5, 0.8))


# ax.axvline(x = np.sqrt(0.62 * 0.28)/2/50*2*1e3)
# xx = np.arange(0.04, 2e3, 1e2)
# yy = xx * 2.0
# ax.plot(xx, yy, 'r--', lw=1)
# xx = np.arange(0.04, 2e3 * 1.41, 1e2)
# yy = xx * 1.0
# ax.plot(xx, yy, 'r--', lw=1)
# xx = np.arange(0.04, 2e3 * 2, 1e2)
# yy = xx * 0.5
# ax.plot(xx, yy, 'r--', lw=1)

# xx = np.arange(8e3, 1e5, 1e2)
# yy = xx * 2.0
# ax.plot(xx, yy, 'r--', lw=1)
# xx = np.arange(8e3 * 1.41, 1e5, 1e2)
# yy = xx * 1.0
# ax.plot(xx, yy, 'r--', lw=1)
# xx = np.arange(8e3 * 2, 1e5, 1e2)
# yy = xx * 0.5
# ax.plot(xx, yy, 'r--', lw=1)


# textposx = 4000
# ax.text(textposx, textposx * 2,
#         '$S_{\\rm 3\\sigma} = 2S_{\\rm p}$',
#         size=8, rotation=45.,
#         ha="center", va="center")
# ax.text(textposx * 1.41, textposx * 1.41,
#         '$S_{\\rm 3\\sigma} = S_{\\rm p}$',
#         size=8, rotation=45.,
#         ha="center", va="center")
# ax.text(textposx * 2, textposx,
#         '$S_{\\rm 3\\sigma} = 0.5S_{\\rm p}$',
#         size=8, rotation=45.,
#         ha="center", va="center")


# ax.text(0.4, 2e4,
# '96 GHz',
# size=10, rotation=0.,
# ha="left", va="top")

plt.savefig('../plot/006-096-flx_2.pdf')
plt.clf()
