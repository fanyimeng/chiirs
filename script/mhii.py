from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy import coordinates
from scipy import stats
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from astropy.table import QTable, Table, Column
from matplotlib.ticker import *
import seaborn as sns
import pandas as pd


def readtable(tablename):
    import pandas as pd
    df = pd.read_csv(tablename,
                     header=None, delim_whitespace=True)
    x = (df[0].to_numpy())
    y = (df[1].to_numpy())
    a = (df[2].to_numpy())
    b = (df[3].to_numpy())
    theta = (df[4].to_numpy())
    r = (df[5].to_numpy())
    flx_iso = (df[6].to_numpy())
    flx_err = (df[7].to_numpy())
    data = {}
    data['x'] = x
    data['y'] = y
    data['a'] = a
    data['b'] = b
    data['theta'] = theta
    data['r'] = r
    data['flx_iso'] = flx_iso
    data['flx_err'] = flx_err
    return data


def coredist(x1, y1, x2, y2):
    if (len(x1) != len(y1)) or (len(x2) != len(y2)):
        return 0
    else:
        dist = []
        for i in range(len(x1)):
            xtmp = x1[i] * np.ones_like(x2)
            xtmp = np.array(xtmp) - np.array(x2)
            xtmp = np.power(xtmp, 2)
            ytmp = y1[i] * np.ones_like(y2)
            ytmp = np.array(ytmp) - np.array(y2)
            ytmp = np.power(ytmp, 2)
            dist.append(np.nanmin(np.power(xtmp + ytmp, 0.5)))
    return(np.array(dist) * 0.05)


def mhii(s, nu, te, d, a, b):
    mhii = 3.15e-12 * s**0.5 * nu**0.5 * \
        te**0.175 * d**2.5 * (a * b)**(1.5 / 2)
    return(mhii)


def kdeplot(data, kdecolor, kdeorder, kdeax, text):
    xx = np.arange(np.nanmin(data),
                   np.nanmax(data),
                   0.001)
    xx = np.log10(xx)
    kde = stats.gaussian_kde(np.log10(data))
    facecolor = (*kdecolor[:3], 0.2)
    ax.plot(np.power(10, xx), kde(xx), color=kdecolor,
            zorder=kdeorder)
    ax.fill_between(np.power(10, xx), kde(xx), 0,
                    facecolor=facecolor,  # The fill color
                    alpha=0.2,
                    zorder=kdeorder)
    ax.plot(data, np.full_like(data, -0.03 * kdeorder - 0.02), '|',
            markeredgewidth=1, color=kdecolor)
    # vlinecolor = (*kdecolor[:3], 0)
    # ax.vlines(x=np.nanmean(data), ymax=2, ymin=1.2, lw=1,
    #           linestyle='--', color=vlinecolor)
    # ax.vlines(x=np.nanmean(data), ymax=1.0, ymin=0, lw=1,
    #           linestyle='--', color=vlinecolor)
    ax.text(np.nanmean(data), 0.57,
            '$\\leftarrow$%.1f $M_{\\odot}$' % (np.nanmean(data)),
            size=10, rotation=90.,
            ha="center", va="bottom", color = kdecolor)
    ax.text(0.004, 0.9 - 0.05 * kdeorder, text, size=10, rotation=0.,
            ha="left", va="top",
            bbox=dict(boxstyle="square",
                      ec=kdecolor,
                      fc=facecolor))

    return [np.power(10, xx), kde(xx)]


plt.rcParams['figure.dpi'] = 75.
plt.rcParams['savefig.dpi'] = 75.
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.clf()
fig = plt.figure(figsize=(4, 4.))


ax = plt.subplot(1, 1, 1)

ax.set_xscale("log", nonposx='clip')
ax.xaxis.set_major_formatter(LogFormatter())
ax.xaxis.set_minor_formatter(LogFormatter())
for axis in [ax.xaxis, ax.yaxis]:
    axis.set_major_formatter(ScalarFormatter())
ax.set_xlabel('$M_\\mathrm{Hii}$ [$M_{\\odot}$]')
ax.set_ylabel('Kernel Density')
ax.tick_params(axis='both', which='major', direction='in', size=10,
               bottom=True, top=True, left=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', size=6,
               bottom=True, top=True, left=True, right=True)

df = pd.read_csv('./cores_006_em.csv')

mhii_all = mhii(s=df['flx_iso'],
                nu=6.,
                te=1e4,
                d=8.3e3,
                a=df['ella'] * 0.05,
                b=df['ellb'] * 0.05)

kde_t1 = kdeplot(data=mhii_all, kdeorder=0,
                 kdecolor=(0.5, 0.5, 0.5, 0.8), kdeax=ax,
                 text='$M_\\mathrm{Hii}$ (All)')
flx = df['flx_iso'].copy()[df['em'] > 1e7]
ella = df['ella'].copy()[df['em'] > 1e7]
ellb = df['ellb'].copy()[df['em'] > 1e7]
mhii_all = mhii(s=flx,
                nu=6.,
                te=1e4,
                d=8.3e3,
                a=ella * 0.05,
                b=ellb * 0.05)
kde_t2 = kdeplot(data=mhii_all, kdeorder=2,
                 kdecolor=(0.5, 0.5, 1.0, 0.8), kdeax=ax,
                 text='$M_\\mathrm{Hii}$ (UCHii)')


ax.set_xlim([0.002, 5])
kde_join = np.concatenate((kde_t1[1], kde_t2[1]))
ax.set_ylim([-0.08, 1.])

ax.axhline(y=0, lw=1, linestyle='-')

# ax.vlines(x=0.65, ymax=2, ymin=0, lw=1,
#           linestyle='--', color='k')
# ax.vlines(x=0.65 * 2, ymax=2, ymin=0, lw=1,
#           linestyle='--', color='gray')

# ax.set_aspect(np.log(ax.get_xlim()[
#               1] / ax.get_xlim()[0]) / np.log(ax.get_ylim()[1] / ax.get_ylim()[0]))
plt.savefig('../plot/mhii.pdf')
plt.clf()
