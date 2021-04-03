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

f_factor = 4e3 * np.log(2.0) * 0.05 * 0.05 / (
    3.14159 * 0.62 * 0.28)


ax.set_xscale("log", nonposx='clip')
ax.xaxis.set_major_formatter(LogFormatter())
ax.xaxis.set_minor_formatter(LogFormatter())
ax.set_yscale("log", nonposy='clip')
ax.yaxis.set_major_formatter(LogFormatter())
ax.yaxis.set_minor_formatter(LogFormatter())


print(datatable['XPEAK_IMAGE'][np.argmin(datatable['FLUX_ISO'])])
print(datatable['YPEAK_IMAGE'][np.argmin(datatable['FLUX_ISO'])])
print(datatable['FLX_PEAK1PIX_006'][np.argmin(datatable['FLUX_ISO'])])
print(datatable['FLX_PEAK3SIGMA_006'][np.argmin(datatable['FLUX_ISO'])])
print(datatable['A_IMAGE'][np.argmin(datatable['FLUX_ISO'])])
print(datatable['B_IMAGE'][np.argmin(datatable['FLUX_ISO'])])
print(datatable['THETA_IMAGE'][np.argmin(datatable['FLUX_ISO'])])
print(np.min(datatable['FLUX_ISO']), )

def plotscatter(ax, x, y, subsize, color, datatable):
    yarr = datatable['FLUX_ISO'].copy()
    xarr = datatable['FLX_PEAK3SIGMA_006'].copy()
    ax.set_xlim([np.nanmin(xarr) * 0.5, np.nanmax(xarr) * 10])
    ax.set_ylim([np.nanmin(yarr) * 0.5, np.nanmax(yarr) * 10])
    ax.errorbar(xarr, yarr,
                yerr=0 ,
                xerr=0,
                ls='none',
                marker='o', markersize='2', lw=0.1,
                ecolor=color,
                color=color,
                capsize=0,
                zorder=10000)
    return 0


plotscatter(ax, 4380, 2050, 750000, (0.5, 0.5, 0.5, 0.5), datatable)

ax.set_aspect(np.log(ax.get_xlim()[
              1] / ax.get_xlim()[0]) / np.log(ax.get_ylim()[1] / ax.get_ylim()[0]))

ax.set_ylabel('$S_\\mathrm{cen}$(6 cm ISO) [mJy]')
ax.set_xlabel('$S_\\mathrm{cen}$(6 cm 3SIGMA) [mJy]')

ax.tick_params(axis='both', which='major', direction='in', size=10,
               bottom=True, top=True, left=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', size=6,
               bottom=True, top=True, left=True, right=True)

xx = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], 1)
yy = xx
ax.plot(xx, yy, '-', lw=0.6, color='gray')
# yy = xx * (10. / 6.)**-0.1
# ax.plot(xx, yy, '-', lw=0.6, color='gray')
# yy = xx * (10. / 6.)**2
# ax.plot(xx, yy, '-', lw=0.6, color='gray')

ax.text(100, 1, 'N', color=(0.9, 0.3, 0.3, 0.5))
ax.text(100, 1 / 2, 'M', color=(0.3, 0.9, 0.3, 0.5))
ax.text(100, 1 / 4, 'S', color=(0.3, 0.3, 0.9, 0.5))
ax.text(100, 1 / 8, 'DS', color=(0.5, 0.5, 0.5, 0.5))

plt.savefig('../plot/2beam-iso-flx.pdf')
plt.clf()
