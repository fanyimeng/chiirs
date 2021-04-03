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

datatable = pd.read_csv('cores006.csv')
a = datatable['A_IMAGE'] * 0.05
b = datatable['B_IMAGE'] * 0.05
ax = plt.subplot(1, 1, 1)


ax.set_xscale("log", nonposx='clip')
ax.xaxis.set_major_formatter(LogFormatter())
ax.xaxis.set_minor_formatter(LogFormatter())
ax.set_yscale("log", nonposy='clip')
ax.yaxis.set_major_formatter(LogFormatter())
ax.yaxis.set_minor_formatter(LogFormatter())


# ax.scatter(datatable['ag_peakflux'], datatable['096f'])


xarr = (np.sqrt(a * b) / 50. * 2.) * 1e3
yarr = datatable['FLUX_ISO'].copy()

ax.set_xlim([np.nanmin(xarr) * 0.5, np.nanmax(xarr) * 4])
ax.set_ylim([np.nanmin(yarr) * 0.5, np.nanmax(yarr) * 4])


cmap = plt.cm.get_cmap('jet', 5)


im = ax.scatter(xarr, yarr,
                c=np.log(datatable['FLX_PEAK1PIX_010'] /
                         datatable['FLX_PEAK1PIX_006']) / np.log(10. / 6.),
                cmap=cmap,
                vmin = -1,
                vmax = 1,
                alpha=0.5,
                zorder=10000)

mask = np.log(datatable['FLX_PEAK1PIX_010'] /
                         datatable['FLX_PEAK1PIX_006']) / np.log(10. / 6.)>-0.1
print("Thermal", np.nansum(np.ones_like(xarr)[mask]))

mask = np.log(datatable['FLX_PEAK1PIX_010'] /
                         datatable['FLX_PEAK1PIX_006']) / np.log(10. / 6.)>-1e20
print("Total", np.nansum(np.ones_like(xarr)[mask]))




# yarr[datatable['em'] < 1e7] = np.nan
# xarr[datatable['em'] < 1e7] = np.nan
# ax.errorbar(xarr, yarr,
#             yerr=datatable['flx_err'] * 1, xerr=(a - b) / 2 * 1e0,
#             ls='none',
#             marker='o', markersize='0.5', lw=0.1,
#             ecolor=(0.3, 0.3, 0.8, 0.8),
#             color=(0.3, 0.3, 0.8, 0.8),
#             capsize=2,
#             zorder=10000)
# ax.errorbar(xarr, yarr,
#                 xerr=3*datatable['cene_006'], yerr=3*datatable['cene_010'],
#                 fmt='none',
#                 color=(0.5, 0.5, 0.5, 0.5))

# xarr = (np.sqrt(a * b) / 50. * 2.) * 1e3
# yarr = datatable['flx_iso'].copy()
# yarr[datatable['em'] > 1e7] = np.nan
# xarr[datatable['em'] > 1e7] = np.nan
# ax.errorbar(xarr, yarr,
#             yerr=datatable['flx_err'] * 1, xerr=(a - b) / 2 * 1e0,
#             ls='none',
#             marker='none', markersize='0', lw=0.1,
#             ecolor=(0.5, 0.5, 0.5, 0.4),
#             color=(0.5, 0.5, 0.5, 0.4),
#             capsize=2,
#             zorder=10000)
# ax.scatter(xarr, yarr,
#            c=yarr / xarr,
#            cmap=plt.cm.jet,
#            alpha=0.5,
#            zorder=10000)

ax.set_aspect(np.log(ax.get_xlim()[
              1] / ax.get_xlim()[0]) / np.log(ax.get_ylim()[1] / ax.get_ylim()[0]))

ax.set_xlabel('$R$ [$10^{-3}$pc]')
ax.set_ylabel('$S_\\mathrm{6cm}$ [mJy]')

ax.tick_params(axis='both', which='major', direction='in', size=10,
               bottom=True, top=True, left=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', size=6,
               bottom=True, top=True, left=True, right=True)

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
        ha='right',
        va='top',
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
        ha='right',
        va='top',
        color=taucolor,
        fontsize=8)
xx = np.arange(ax.get_xlim()[0], 3.5, 1)
yy = tau * 6.**2 * (xx / 1e3 / 2 * 50) ** 2 * 1e4
ax.plot(xx, yy, '--', lw=0.6, color=taucolor)


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

ax.axvline(x = (0.62*0.28)**0.5/50*2*1e3)

plt.subplots_adjust(hspace=0.05)
pcbar = ax.get_position().get_points().flatten()
cbar_ax = fig.add_axes(
    [pcbar[2] + 0.01 * (pcbar[2] - pcbar[0]),
     pcbar[1],
     0.02,
     pcbar[3] - pcbar[1]])
cb = plt.colorbar(im, cax=cbar_ax, orientation='vertical')
cbar_ax.yaxis.set_ticks_position('right')
cbar_ax.yaxis.set_label_position('right')
cbar_ax.tick_params(axis='both', which='both', direction='in')
cb.set_label("$\\alpha$")


plt.savefig('../plot/r-flx2.pdf')
plt.clf()
