import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import *

plt.rcParams['figure.dpi'] = 75.
plt.rcParams['savefig.dpi'] = 75.
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.clf()
fig = plt.figure(figsize=(6.0, 6.0))
ax = plt.subplot(1, 1, 1)

datatable = pd.read_csv('cores006.csv')

yarr = datatable['FLX_PEAK_010'].copy()
xarr = datatable['FLX_PEAK_006'].copy()
ratio = yarr / xarr
xlimarr = xarr.copy()
ylimarr = yarr.copy()
xlimarr[xarr < 0.0001] = 0.0001
ylimarr[yarr < 0.0001] = 0.0001
xarr[xarr < 0.0001] = 0.0001
yarr[yarr < 0.0001] = 0.0001

ax.set_xscale("log", nonposx='clip')
ax.xaxis.set_major_formatter(LogFormatter())
ax.xaxis.set_minor_formatter(LogFormatter())
ax.set_yscale("log", nonposy='clip')
ax.yaxis.set_major_formatter(LogFormatter())
ax.yaxis.set_minor_formatter(LogFormatter())

print([np.nanmin(xarr) * 0.5, np.nanmax(xarr) * 2])
print([np.nanmin(yarr) * 0.5, np.nanmax(yarr) * 2])
ax.set_xlim([np.nanmin(xlimarr) * 0.5, np.nanmax(xlimarr) * 10])
ax.set_ylim([np.nanmin(ylimarr) * 0.5, np.nanmax(ylimarr) * 10])


def plotscatter(ax, x, y, subsize, color, datatable):
    zx1 = x - subsize
    zx2 = x + subsize
    zy1 = y - subsize
    zy2 = y + subsize
    yarr = datatable['FLX_PEAK_010'].copy()
    xarr = datatable['FLX_PEAK_006'].copy()
    xarr[xarr < 0.0001] = 0.0001
    yarr[yarr < 0.0001] = 0.0001
    # yarr[datatable['XPEAK_IMAGE'] < zx1] = np.nan
    # xarr[datatable['XPEAK_IMAGE'] < zx1] = np.nan
    # yarr[datatable['XPEAK_IMAGE'] > zx2] = np.nan
    # xarr[datatable['XPEAK_IMAGE'] > zx2] = np.nan
    # yarr[datatable['YPEAK_IMAGE'] < zy1] = np.nan
    # xarr[datatable['YPEAK_IMAGE'] < zy1] = np.nan
    # yarr[datatable['YPEAK_IMAGE'] > zy2] = np.nan
    # xarr[datatable['YPEAK_IMAGE'] > zy2] = np.nan
    ax.errorbar(xarr, yarr,
                yerr=datatable['ERR_PEAK_010'],
                xerr=datatable['ERR_PEAK_006'],
                ls='none',
                marker='o', markersize='2', lw=0.1,
                ecolor=color,
                color=color,
                capsize=0,
                zorder=10000)
    return 0


# plotscatter(ax, 4780, 5780, 350, (0.9, 0.3, 0.3, 0.5), datatable)
# plotscatter(ax, 4700, 4860, 350, (0.3, 0.9, 0.3, 0.5), datatable)
# plotscatter(ax, 4700, 3900, 350, (0.3, 0.3, 0.9, 0.5), datatable)
plotscatter(ax, 4380, 2050, 750000, (0.5, 0.5, 0.5, 0.5), datatable)

ax.set_aspect(np.log(ax.get_xlim()[
              1] / ax.get_xlim()[0]) / np.log(ax.get_ylim()[1] / ax.get_ylim()[0]))

ax.set_ylabel('$S_\\mathrm{cen}$(3 cm) [mJy]')
ax.set_xlabel('$S_\\mathrm{cen}$(6 cm) [mJy]')

ax.tick_params(axis='both', which='major', direction='in', size=10,
               bottom=True, top=True, left=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', size=6,
               bottom=True, top=True, left=True, right=True)

xx = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], 1)
yy = xx * (10. / 6.)**-0.5
ax.plot(xx, yy, '-', lw=0.6, color='gray')
yy = xx * (10. / 6.)**-0.1
ax.plot(xx, yy, '-', lw=0.6, color='gray')
yy = xx * (10. / 6.)**2
ax.plot(xx, yy, '-', lw=0.6, color='gray')

ax.text(100, 1, 'N', color=(0.9, 0.3, 0.3, 0.5))
ax.text(100, 1 / 2, 'M', color=(0.3, 0.9, 0.3, 0.5))
ax.text(100, 1 / 4, 'S', color=(0.3, 0.3, 0.9, 0.5))
ax.text(100, 1 / 8, 'DS', color=(0.5, 0.5, 0.5, 0.5))

plt.savefig('../plot/006-010-flx.pdf')
plt.clf()
