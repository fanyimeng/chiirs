import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import pandas as pd
from tqdm import tqdm
from matplotlib.ticker import *
import scipy.optimize as opt


def tau(te=1e4,
        nu=6,
        em=1e7):
    return 8.235e-2 * te**(-1.35) * nu**(-2.1) * em


def em(te=1e4, nu=6,
        tau=10):
    return tau / (8.235e-2 * te**(-1.35) * nu**(-2.1))


def tb(tau=1, te=1e4):
    return te * (1 - np.exp(-1 * tau))


def s2t(s=10, a=0.65, b=0.65, nu=6.):
    return 1.222e3 * s / (nu**2 * a * b)


def t2s(t=10, a=0.65, b=0.65, nu=6.):
    return t * (nu**2 * a * b) / 1.222e3


def nly(nu=6,
        s=1,  # mJy
        te=1e4,
        d=8.34e3):
    mynly = 4.771e39 * s * te**-0.45 * nu**0.1 * d**2
    return mynly


def solveSed(f6,
             f22,
             parameters,
             init_val=[1e9, 0.3]):
    def f(variables, parameters=parameters, f6=f6, f22=f22):
        s06 = t2s(t=tb(tau=tau(te=parameters[0],
                               nu=6.,
                               em=variables[0]),
                       te=parameters[0]),
                  a=variables[1], b=variables[1],
                  nu=6.) - f6
        s22 = t2s(t=tb(tau=tau(te=parameters[0],
                               nu=22.4,
                               em=variables[0]),
                       te=parameters[0]),
                  a=variables[1], b=variables[1],
                  nu=22.4) - f22
        return [s06, s22]
    solution = opt.fsolve(f, init_val)
    return solution


def solveSed_s6(f6,
                parameters):
    def f(variables, parameters=parameters, f6=f6):
        s06 = t2s(t=tb(tau=tau(te=parameters[0],
                               nu=6.,
                               em=variables),
                       te=parameters[0]),
                  a=parameters[1], b=parameters[1],
                  nu=6.) - f6
        return s06
    solution = opt.fsolve(f, 1e7)
    return solution


def ionsed(nu, s40):
    return s40 * (nu / 40.)**-0.1


def dustSed(nu,  # in GHz
            md,
            td,  # in K
            k0,  # k100GHz here
            nu0,  # 100GHz here
            beta):
    h = 6.626e-27  # cgs
    c = 3e10  # cgs
    kb = 1.38e-16  # erg/k

    snu = 2 * h * (nu * 1e9)**3 / (c**2)
    snu = snu / (np.exp(h * nu * 1e9 / kb / td) - 1)
    # snu = snu * (1 - np.exp(-1 * k0 * (nu / nu0)**beta * md)) * omega
    snu = snu * k0 * (nu / nu0)**beta * md / (8.34e3 * 3e18)**2
    # snu = md/0.00325/nu**(-3)/(k0 * (nu / nu0)**beta)
    return snu * 1e26


def solveSed_dust(f96,
                  parameters):
    #solveSed(flist[0], flist[5],  [100., 0.091, 100., 1.75])
    def f(variables, parameters=parameters, f96=f96):
        s96 = dustSed(96, variables, *parameters) \
            - f96
        return s96
    solution = opt.fsolve(f, 2e34)
    return solution


def photometry(df,
               fitsimage,
               rmsfits='',
               colname=['flx', 'err']):
    hdu = fits.open(fitsimage)[0]
    data = hdu.data
    data[data < 0] = 0
    header = hdu.header

    flx_factor = 4 * np.log(2.0) * header['CDELT2'] * header['CDELT1'] / (
        3.14159 * header['BMAJ'] * header['BMIN'])
    flx_factor = abs(flx_factor * 1e3)

    x1d = np.arange(0, data.shape[0])
    y1d = np.arange(0, data.shape[1])
    x2d = x1d[np.newaxis, :]
    y2d = y1d[:, np.newaxis]

    df[colname[0]] = np.ones(df.shape[0]) * -1
    if rmsfits != '':
        rms_data = fits.open(rmsfits)[0].data
        df[colname[1]] = np.ones(df.shape[0]) * -1

    for i in tqdm(range(df.shape[0])):
        row = df.iloc[i, :]
        mask = ((x2d - row['x_img'])**2 +
                (y2d - row['y_img'])**2 <= row['r_img']**2)
        df[colname[0]][i] = data[mask].sum() * flx_factor
        if rmsfits != '':
            num_pix = 3.14 * row['r_img']**2
            df[colname[1]][i] = rms_data[mask].sum() * flx_factor / \
                np.sqrt(num_pix)
    return df


'''
# for band in ['006', '022', '096']:
#     df = pd.read_csv('%s.csv' % band)
#     df = photometry(df=df,
#                     fitsimage='../data/img_%s.fits' % band,
#                     rmsfits='../data/rms_%s.fits' % band,
#                     colname=['flx_%s' % band, 'err_%s' % band]
#                     )

#     df.to_csv('%s.csv' % band, index=False)

#     df = pd.read_csv('%s.csv' % band)
#     ax = plt.subplot(1, 1, 1)

#     ax.set_xscale("log", nonposx='clip')
#     ax.xaxis.set_major_formatter(LogFormatter())
#     ax.xaxis.set_minor_formatter(LogFormatter())
#     # ax.set_yscale("log", nonposy='clip')
#     # ax.yaxis.set_major_formatter(LogFormatter())
#     # ax.yaxis.set_minor_formatter(LogFormatter())

#     # plt.scatter(df['flx_%s' % band], df['y_img'])
#     plt.errorbar(x=df['flx_%s' % band],
#                  y=df['y_img'],
#                  xerr=df['err_%s' % band],
#                  yerr=0,
#                  ls='none',
#                  marker='o', markersize='0.5', lw=1,)

#     plt.savefig('../plot/photometry_test_%s.pdf' % band)












df006 = pd.read_csv('006.csv')
df022 = pd.read_csv('022.csv')

df006['EM'] = np.zeros(df006.shape[0]) - 1
df006['R'] = np.zeros(df006.shape[0]) - 1
df006['f022'] = np.zeros(df006.shape[0]) - 1
df006['w022'] = np.zeros(df006.shape[0])
for i in range(df006.shape[0]):
    row006 = df006.iloc[i, :]
    if int(row006['022_idx']) > -0.5:
        row022 = df022.iloc[int(row006['022_idx']), :]
        cond1 = row022['flx_022'] >= row006['flx_006'] * (22.4 / 6.)**-0.1
        cond2 = row022['flx_022'] <= row006['flx_006'] * (22.4 / 6.)**2
        if cond1 and cond2:
            init_r = np.sqrt(
                (row006['r_img'] * 0.05)**2 - (0.62 / 2) * (0.28 / 2))
            init_em = solveSed_s6(
                row006['flx_006'], parameters=[1e4, init_r])[0]
            print(init_r, init_em)
            df006['EM'][i] = solveSed(
                row006['flx_006'], row022['flx_022'], [1e4], [init_em, init_r])[0]
            df006['R'][i] = solveSed(
                row006['flx_006'], row022['flx_022'], [1e4], [init_em, init_r])[1]
            df006['f022'][i] = row022['flx_022']
            df006['w022'][i] = 1
        else:
            df006['R'][i] = np.sqrt(
                (row006['r_img'] * 0.05)**2 - (0.62 / 2) * (0.28 / 2))
            df006['EM'][i] = solveSed_s6(
                f6=row006['flx_006'], parameters=[1e4, df006['R'][i]])
            print('R AND EM', df006['R'][i], df006['EM'][i])
            df006['f022'][i] = row022['flx_022']
    else:
        df006['R'][i] = np.sqrt(
            (row006['r_img'] * 0.05)**2 - (0.62 / 2) * (0.28 / 2))
        df006['EM'][i] = solveSed_s6(
            f6=row006['flx_006'], parameters=[1e4, df006['R'][i]])[0]
        print('R AND EM', df006['R'][i], df006['EM'][i])

df006.to_csv('006.csv', index=False)

df = pd.read_csv('006.csv')

ax = plt.subplot(1, 1, 1)

ax.set_xscale("log", nonposx='clip')
ax.xaxis.set_major_formatter(LogFormatter())
ax.xaxis.set_minor_formatter(LogFormatter())
ax.set_yscale("log", nonposy='clip')
ax.yaxis.set_major_formatter(LogFormatter())
ax.yaxis.set_minor_formatter(LogFormatter())

df['R'][df['R'] < 0] = np.nan

ax.set_xlim([1e-1, 10])
ax.set_ylim([1e-1, 10])
ax.set_aspect(np.log(ax.get_xlim()[
              1] / ax.get_xlim()[0]) / np.log(ax.get_ylim()[1] / ax.get_ylim()[0]))

for i in range(df.shape[0]):
    plt.text(df['r_img'][i] * 0.05, df['R'][i] + i * 0, df['006_idx'][i], size=10,
             color=(0.3, 0.3, 0.3, 0.3))
# plt.errorbar(x=df['flx_%s' % band],
#              y=df['y_img'],
#              xerr=df['err_%s' % band],
#              yerr=0,
#              ls='none',
#              marker='o', markersize='0.5', lw=1,)

xx = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], 0.001)
yy = np.sqrt(xx**2 - 0.62 * 0.28)

ax.plot(xx, yy)

plt.savefig('../plot/R_test.pdf')

plt.clf()
ax = plt.subplot(1, 1, 1)
ax.set_xlim([4e3, 6e3])
ax.set_ylim([4e3, 6e3])

for i in range(df.shape[0]):
    plt.text(df['x_img'][i], df['y_img'][i], df['006_idx'][i], size=1)
plt.savefig('../plot/Core_pos_test.pdf')
print(np.nanmax(df['EM']))
plt.clf()


plt.clf()
ax2 = plt.subplot(1, 1, 1)
df = pd.read_csv('006.csv')
ax2.set_xlim([5, 11])
ax2.set_ylim([0, 1e4])
for i in range(df.shape[0]):
    print('LOG:', np.log10(df['EM'][i]))
    ax2.text(np.log10(df['EM'][i]), df['y_img'][i], df['006_idx'][i], size=10)
plt.savefig('../plot/em_test.pdf')
print(np.nanmax(df['EM']))
plt.clf()


df = pd.read_csv('006.csv')
df['NLY'] = np.ones(df.shape[0]) - 1
df['f080'] = np.ones(df.shape[0]) - 1
for i in range(df.shape[0]):
    row = df.iloc[i, :]
    tau80 = tau(nu=80, em=row['EM'])
    tb80 = tb(tau=tau80)
    s80 = t2s(t=tb80, a=row['R'], b=row['R'], nu=80)
    df['f080'][i] = s80
    nly80 = nly(nu=80, s=s80)
    print(np.log10(nly80))
    df['NLY'][i] = np.log10(nly80)
df.to_csv('006.csv', index=False)

'''
df = pd.read_csv('006.csv')
df096 = pd.read_csv('096.csv')
df['w096'] = np.ones(df.shape[0])
df['obsf_096'] = np.ones(df.shape[0]) * -1
df['dustf_096'] = np.ones(df.shape[0]) * -1
df['gasmass'] = np.ones(df.shape[0]) * -1
df['gasn'] = np.ones(df.shape[0]) * -1
df['r_i'] = np.ones(df.shape[0]) * -1
df['t'] = np.ones(df.shape[0]) * -1
for i in range(df.shape[0]):
    row = df.iloc[i, :]
    row096 = df096.iloc[int(row['096_idx']), :]
    df['obsf_096'][i] = row096['flx_096']
    tau96 = tau(nu=96, em=row['EM'])
    tb96 = tb(tau=tau96)
    s96 = t2s(t=tb96, a=row['R'], b=row['R'], nu=96)
    if row096['flx_096'] - s96 > 0 and row096['flx_096'] - s96 > 3 * row096['err_096']:
        dust_r = np.sqrt((row096['r_img'] * 0.05)**2)
        df['dustf_096'][i] = row096['flx_096'] - s96
        df['gasmass'][i] = solveSed_dust(
            df['dustf_096'][i], [100, 2.631, 100., 1.05]) * 100 / 2e33
        volume = 3.14 * (dust_r / 50 * 2 * 3e18)**3
        density = df['gasmass'][i] * 2e33 / volume
        density = density / 2.33 / 1.7e-24
        df['gasn'][i] = density
        df['r_i'][i] = 1.99e-2 * \
            (10**row['NLY'] / 1e49)**(1 / 3.) * \
            (df['gasn'][i] / 1e5)**(-2 / 3.)
        lifetime = (row['R'] / 50 * 2 / df['r_i'][i])**(7 / 4.)
        lifetime = lifetime - 1
        ci = 10 / 3e13 * 3.14e7
        lifetime = lifetime * 4 * df['r_i'][i] / 7 / ci
    if df['gasn'][i] < 1e7:
        df['r_i'][i] = 1.99e-2 * \
            (10**row['NLY'] / 1e49)**(1 / 3.) * \
            (1e7 / 1e5)**(-2 / 3.)
        lifetime = (row['R'] / 50 * 2 / df['r_i'][i])**(7 / 4.)
        lifetime = lifetime - 1
        ci = 10 / 3e13 * 3.14e7
        lifetime = lifetime * 4 * df['r_i'][i] / 7 / ci
    df['t'][i] = lifetime

df.to_csv('006.csv', index=False)


fig = plt.figure(figsize=(4, 16 / 3))
plt.clf()
ax2 = plt.subplot(1, 1, 1)
df = pd.read_csv('006.csv')
df['t'][df['t'] < 100] = np.nan
df['x_img'][df['t'] < 100] = np.nan
ax2.set_xlim([3e3, 7e3])
ax2.set_ylim([3e3, 7e3])
im = ax2.scatter(df['x_img'], df['y_img'],
                 c=np.log10(df['t']),
                 s=0.1,
                 cmap=plt.cm.jet
                 )
pcbar = ax2.get_position().get_points().flatten()
cbar_ax = fig.add_axes(
    [pcbar[2] + 0.01 * (pcbar[2] - pcbar[0]),
     pcbar[1],
     0.02,
     pcbar[3] - pcbar[1]])
cb = plt.colorbar(im, cax=cbar_ax, orientation='vertical')
cbar_ax.yaxis.set_ticks_position('right')
cbar_ax.yaxis.set_label_position('right')
cbar_ax.tick_params(axis='both', which='both', direction='in')
cb.set_label("log(Time/yr)")
# for i in range(df.shape[0]):

# ax2.text(np.log10(df['t'][i]), df['y_img'][i], df['006_idx'][i], size=10)
plt.savefig('../plot/time_test.pdf')
plt.clf()


fig = plt.figure(figsize=(4.5, 4.5))
plt.clf()
ax2 = plt.subplot(1, 1, 1)
ax2.set_xscale("log", nonposx='clip')
# ax2.xaxis.set_major_formatter(LogFormatter())
# ax2.xaxis.set_minor_formatter(LogFormatter())
ax2.set_yscale("log", nonposy='clip')
# ax2.yaxis.set_major_formatter(LogFormatter())
# ax2.yaxis.set_minor_formatter(LogFormatter())

df = pd.read_csv('006.csv')

# df['r_i'][df['gasn'] < 1e7]=np.nan
# ax2.set_ylim([np.nanmin(df['R']) / 25 / 2, np.nanmax(df['R']) / 25 * 2])
# ax2.set_xlim([(4 * 1e4 / 100)**(2 / 3.) * np.nanmin(df['r_i']) / 2,
#               (4 * 1e4 / 100)**(2 / 3.) * np.nanmax(df['r_i']) * 2])
# im = ax2.scatter((4 * 1e4 / 100)**(2 / 3.) * df['r_i'], df['R'] / 25,
#                  c='r',
#                  s=1,
#                  # cmap=plt.cm.jet
#                  )
# xx = np.arange(-6, 0, 0.01)
# xx = 10**xx
# yy = xx
# ax2.plot(xx,yy)
# pcbar = ax2.get_position().get_points().flatten()
# cbar_ax = fig.add_axes(
#     [pcbar[2] + 0.01 * (pcbar[2] - pcbar[0]),
#      pcbar[1],
#      0.02,
#      pcbar[3] - pcbar[1]])
# cb = plt.colorbar(im, cax=cbar_ax, orientation='vertical')
# cbar_ax.yaxis.set_ticks_position('right')
# cbar_ax.yaxis.set_label_position('right')
# cbar_ax.tick_params(axis='both', which='both', direction='in')
# cb.set_label("log(EM)")

x = df['gasn']
x[x < 0] = np.nan
y = np.sqrt(df['EM'] / (df['R'] / 25.))
ax2.set_xlim([np.nanmin(x) / 2, np.nanmax(x) * 2])
ax2.set_ylim([np.nanmin(y) / 2, np.nanmax(y) * 2])

ax2.set_aspect(np.log(ax2.get_xlim()[
    1] / ax2.get_xlim()[0]) / np.log(ax2.get_ylim()[1] / ax2.get_ylim()[0]))


im = ax2.scatter(x, y, color='r')

xx = np.arange(1e6, 1e11, 1e6)
yy = xx / 100 / 4.

ax2.plot(xx, yy, label = 'Equilibrium')
ax2.legend()
ax2.set_xlabel('$n_{\\rm H_2}$ [$\\rm cm^{-3}$]')
ax2.set_ylabel('$n_{\\rm e}$ [$\\rm cm^{-3}$]')

# xx = np.arange(-6, 0, 0.01)
# for i in range(df.shape[0]):

# ax2.text(np.log10(df['t'][i]), df['y_img'][i], df['006_idx'][i], size=10)
plt.savefig('../plot/pressure_qual_test.pdf')
plt.clf()


'''
(nu=6,
        s=1,  # mJy
        te=1e4,
        d=8.34e3)
def tau(te=1e4,
        nu=6,
        em=1e7):
    return 8.235e-2 * te**(-1.35) * nu**(-2.1) * em


def em(te=1e4, nu=6,
        tau=10):
    return tau / (8.235e-2 * te**(-1.35) * nu**(-2.1))


def tb(tau=1, te=1e4):
    return te * (1 - np.exp(-1 * tau))


def s2t(s=10, a=0.65, b=0.65, nu=6.):
    return 1.222e3 * s / (nu**2 * a * b)


def t2s(t=10, a=0.65, b=0.65, nu=6.):
    return t * (nu**2 * a * b) / 1.222e3
'''
