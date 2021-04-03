import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.ticker import *
from scipy.optimize import fsolve
from scipy.optimize import least_squares

# f = interp1d(data[1], data[0])

df = pd.read_csv('cores006.csv')


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


plt.rcParams['figure.dpi'] = 75.
plt.rcParams['savefig.dpi'] = 75.
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.clf()
fig = plt.figure(figsize=(6.0, 6.0))

df = pd.read_csv('cores006.csv')

ax = plt.subplot(1, 1, 1)


ax.set_xscale("log", nonposx='clip')
ax.xaxis.set_major_formatter(LogFormatter())
ax.xaxis.set_minor_formatter(LogFormatter())
ax.set_yscale("log", nonposy='clip')
ax.yaxis.set_major_formatter(LogFormatter())
ax.yaxis.set_minor_formatter(LogFormatter())

# ax.set_xlim([np.nanmin(xarr) * 0.5, np.nanmax(xarr) * 4])


# tau = np.arange(1e-4, 1000, 0.1)
s = df['FLUX_ISO']
r = np.sqrt(df['A_IMAGE'] * df['B_IMAGE']) * 0.05
r[r < 0.65 / 2] = 0.65 / 2

# yy = t2s(t = tb(tau = tau(nu = xx)), nu  = xx)

emlist = []
for s in [df['FLUX_ISO'],
          df['FLUX_ISO'] + df['FLUXERR_ISO'],
          df['FLUX_ISO'] - df['FLUXERR_ISO']]:
    r = np.sqrt(df['A_IMAGE'] * df['B_IMAGE']) * 0.05
    r[r > 0] = 0.65 / 2

    # yy = t2s(t = tb(tau = tau(nu = xx)), nu  = xx)

    taulist = np.ones_like(s)
    for i in range(len(taulist)):
        def taufunc(x):
            t = s2t(s[i], a=r[i], b=r[i])
            return tb(tau=x) - t
        taulist[i] = fsolve(taufunc, 0.2)
        # print(taulist[i])
    emlist.append(em(tau=taulist))
df['EM_6_065'] = emlist[0]
df['EM_6_065_UP'] = emlist[1]
df['EM_6_065_DOWN'] = emlist[2]


emlist = []
for s in [df['FLUX_ISO'],
          df['FLUX_ISO'] + df['FLUXERR_ISO'],
          df['FLUX_ISO'] - df['FLUXERR_ISO']]:
    r = np.sqrt(df['A_IMAGE'] * df['B_IMAGE']) * 0.05
    # r[r < 0.65 / 2] = 0.65 / 2

    # yy = t2s(t = tb(tau = tau(nu = xx)), nu  = xx)

    taulist = np.ones_like(s)
    for i in range(len(taulist)):
        def taufunc(x):
            t = s2t(s[i], a=r[i], b=r[i])
            return tb(tau=x) - t
        taulist[i] = fsolve(taufunc, 0.2)
        # print(taulist[i])
    emlist.append(em(tau=taulist))
df['EM_6'] = emlist[0]
df['EM_6_UP'] = emlist[1]
df['EM_6_DOWN'] = emlist[2]


# plt.scatter(s, df['EM_6'])
# plt.scatter(s, df['EM_6_DOWN'])
# plt.scatter(s, df['EM_6_UP'])
ax.errorbar(s, df['EM_6'], xerr=0,
            yerr=np.array([np.abs(df['EM_6'] - df['EM_6_DOWN']),
                           np.abs(df['EM_6'] - df['EM_6_UP'])]),
            ls='none',
            marker='o', markersize='2', lw=2,
            capsize=2,
            zorder=10000)

ax.errorbar(s, df['EM_6_065'], xerr=0,
            yerr=np.array([np.abs(df['EM_6_065'] - df['EM_6_065_DOWN']),
                           np.abs(df['EM_6_065'] - df['EM_6_065_UP'])]),
            ls='none',
            marker='o', markersize='2', lw=2,
            capsize=2,
            zorder=10000)

def emfit(frequency, flux):
    def fun(x, freq, y):
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

        tau = tau(te=x[0],
                  nu=freq,
                  em=x[1])
        tb = tb(tau=tau, te=x[0])
        s = t2s(t=tb, nu=freq)
        return s - y
    x0 = np.array([1e4, 1e7])
    res_soft_l1 = least_squares(fun, x0, loss='soft_l1',
                                args=(frequency, flux))
    '''
    https://stackoverflow.com/questions/42388139
    '''
    res = []
    for i in range(len(frequency)):
        res.append(fun(res_soft_l1.x, frequency[i], flux[i]))
    J = res_soft_l1.jac
    # print np.dot(J.T, J)
    cov = np.linalg.inv(np.dot(J.T, J)) * (np.square(res).mean())
    err = np.sqrt(np.diag(cov))
    # print err
    return res_soft_l1.x[1], err[-1]


def emfit_fixte(frequency, flux):
    def fun(x, freq, y):
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

        tau = tau(te=1e4,
                  nu=freq,
                  em=x[0])
        tb = tb(tau=tau, te=1e4)
        s = t2s(t=tb, nu=freq)
        return s - y
    x0 = np.array([1e4, 1e7])
    res_soft_l1 = least_squares(fun, x0, loss='soft_l1',
                                args=(frequency, flux))
    # '''
    # https://stackoverflow.com/questions/42388139
    # '''
    # res = []
    # for i in range(len(frequency)):
    #     res.append(fun(res_soft_l1.x, frequency[i], flux[i]))
    # J = res_soft_l1.jac
    # # print np.dot(J.T, J)
    # cov = np.linalg.inv(np.dot(J.T, J)) * (np.square(res).mean())
    # err = np.sqrt(np.diag(cov))
    # # print err
    return res_soft_l1.x[0]


# s6 = df['FLX_PEAK3SIGMA_006']
# s10 = df['FLX_PEAK3SIGMA_010']

def em_610(s6, s10):
    em2 = np.ones_like(s6)
    for i in range(len(em2)):
        if s10[i] > s6[i] * ((10 / 6)**(-0.1)):
            em2[i] = emfit_fixte(
                np.array([6, 10]), np.array([s6[i], s10[i]]))
            # print(em2[i])
        else:
            em2[i] = -1
    return em2


s6 = df['FLX_PEAK3SIGMA_006']
s10 = df['FLX_PEAK3SIGMA_010']
s6_up = df['FLX_PEAK3SIGMA_006'] + df['ERR_PEAK3SIGMA_006']
s10_up = df['FLX_PEAK3SIGMA_010'] + df['ERR_PEAK3SIGMA_010']
s6_down = df['FLX_PEAK3SIGMA_006'] - df['ERR_PEAK3SIGMA_006']
s10_down = df['FLX_PEAK3SIGMA_010'] - df['ERR_PEAK3SIGMA_010']


df['EM_610'] = em_610(s6=s6, s10=s10)
df['EM_610_UP'] = em_610(s6=s6_down, s10=s10_up)
df['EM_610_DOWN'] = em_610(s6=s6_up, s10=s10_down)

ax.errorbar(s, df['EM_610'],
            xerr=0,
            yerr=np.array([np.abs(df['EM_610'] - df['EM_610_DOWN']),
                           np.abs(df['EM_610'] - df['EM_610_UP'])]),
            alpha = 0.5,
            ls='none',
            marker='o', markersize='2', lw=2,
            capsize=2,
            zorder=100)



'''
df['EM_610'], df['EMERR_610'] = em_610(s6=df['FLX_PEAK3SIGMA_006'],
                                       s10=df['FLX_PEAK3SIGMA_010'],
                                       s6_up=df['FLX_PEAK3SIGMA_006'] +
                                       df['ERR_PEAK3SIGMA_006'],
                                       s10_up=df['FLX_PEAK3SIGMA_010'] +
                                       df['ERR_PEAK3SIGMA_010'],
                                       s6_down=df['FLX_PEAK3SIGMA_006'] -
                                       df['ERR_PEAK3SIGMA_006'],
                                       s10_down=df['FLX_PEAK3SIGMA_010'] -
                                       df['ERR_PEAK3SIGMA_010'])
'''
# plt.scatter(s, em2)

ax.set_xlim([np.nanmin(s) * 0.5, np.nanmax(s) * 4])

print(np.exp(-tau(em=5e8)))

df.to_csv('cores006.csv', index=False)

plt.savefig('../plot/em_test.pdf')
