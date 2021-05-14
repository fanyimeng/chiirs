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
fig = plt.figure(figsize=(4.5, 4.0))

ax = plt.subplot(1, 1, 1)

ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')


ax.set_xlim([5e-3, 1e-1])
ax.set_ylim([5e-3, 1e-1])

ax.set_aspect(np.log(ax.get_xlim()[
              1] / ax.get_xlim()[0]) /
              np.log(ax.get_ylim()[1] /
                     ax.get_ylim()[0]))

df = pd.read_csv('006.csv')
r_calc = df['R'].copy() / 25.
r_calc[df['w022'] < 0.5] = np.nan

df = pd.read_csv('006.csv')
r_obs = df['r_img'].copy() * 0.05 / 25.
ax.scatter(r_obs, r_calc,
           color=(0.3, 0.3, 0.8, 0.5),
           edgecolor='none',
           s=15)

xx = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], 0.001)
yy = np.sqrt(xx**2 - 0.62 * 0.28 / 25.**2)
ax.plot(xx, yy,
        linewidth=1,
        linestyle='--',
        color=(0.3, 0.3, 0.8, 0.5),
        label = 'Deconvolved $r_{\\rm obs6}$')

# df022 = pd.read_csv('022.csv')
# df['r_img_022'] = np.ones(df.shape[0])
# for i in range(df.shape[0]):
#     if df.iloc[i,:]['w022'] > 0.5:
#         row = df022.iloc[int(df.iloc[i,:]['022_idx']), :]
#         df['r_img_022'][i] = row['r_img']
# df['r_img_022'][df['w022'] < 0.5] = np.nan

# r_obs = df['r_img_022'].copy() * 0.05 / 25.
# ax.scatter(r_obs, r_calc,
#            color=(0.8, 0.3, 0.3, 0.5),
#            edgecolor='none',
#            s=15)

# xx = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], 0.001)
# yy = np.sqrt(xx**2 - 0.27 * 0.25 / 25.**2)
# ax.plot(xx, yy,
#         linewidth=1,
#         linestyle='--',
#         color=(0.8, 0.3, 0.3, 0.5),
#         label = 'Deconvolved $r_{\\rm obs22}$')

ax.legend()

ax.set_xlabel('$r_{\\rm obs6}$ [pc]')
ax.set_ylabel('$r_{\\rm calc}$ [pc]')

plt.savefig('../plot/R.pdf')
