import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
from matplotlib.ticker import *
from scipy.optimize import fsolve
from scipy.optimize import least_squares
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


plt.rcParams['figure.dpi'] = 75.
plt.rcParams['savefig.dpi'] = 75.
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.clf()
fig = plt.figure(figsize=(6.0, 6.0))

ax = plt.subplot(3, 1, 1)


ax.set_xscale("log", nonposx='clip')
ax.xaxis.set_major_formatter(LogFormatter())
ax.xaxis.set_minor_formatter(LogFormatter())
ax.set_yscale("log", nonposy='clip')
ax.yaxis.set_major_formatter(LogFormatter())
ax.yaxis.set_minor_formatter(LogFormatter())

nu = np.arange(1, 80, 0.1)

xx = [6., 10., 22.4]
yy = [25, 25, 170]
te = 1e4

mytau = tau(te=te,
            nu=nu,
            em=1.4e9)

mytb = tb(tau=mytau, te=te)

mys = t2s(t=mytb,
          a=0.29, b=0.29,
          nu=nu)


# ax.scatter(xx, yy)
# ax.plot(nu, mys)

ax.set_xlim([np.nanmin(nu) * 0.5, np.nanmax(nu) * 4])


def solveSed(f6,
             f22,
             parameters):
    def f(variables, parameters=parameters, f6=f6, f22=f22):
        s06 = t2s(t=tb(tau=tau(te=parameters[0],
                               nu=6.,
                               em=variables[0]),
                       te=parameters[0]),
                  a=variables[1], b=variables[1],
                  nu=6.) - f6
        s22 = t2s(t=tb(tau=tau(te=parameters[0],
                               nu=22.,
                               em=variables[0]),
                       te=parameters[0]),
                  a=variables[1], b=variables[1],
                  nu=22.) - f22
        return [s06, s22]
    solution = opt.fsolve(f, (1e7, 0.7))
    return solution


res = solveSed(1000, 1000, [1e4])
print(res)


size = 1000
xx = range(1, size)
yy = range(1, size)

EM, RAD = np.meshgrid(xx, yy)
EM = EM.astype(float)
RAD = RAD.astype(float)
print(EM.shape, RAD.shape)

simu_file = open('s-s_em-r', 'w')
header_text = 's6,s22,em,r\n'
simu_file.write(header_text)

for x in range(EM.shape[0]):
    for y in range(EM.shape[1]):
        s6 = x
        s22 = y
        if s22 > (22 / 6.)**-0.1 * s6 and s22 < (22 / 6.)**2 * s6:
            EM[x, y] = solveSed(s6, s22, [1e4])[0]
            RAD[x, y] = solveSed(s6, s22, [1e4])[1]
            # mystr = '%f,%f,%f,%f\n' % (s6, s22, *solveSed(s6, s22, [1e4]))
            # simu_file.write(mystr)
        else:
            EM[x, y] = 0
            RAD[x, y] = 0
            # mystr = '%f,%f,%f,%f\n' % (s6, s22, 0, 0)
            # simu_file.write(mystr)
simu_file.close()
f_em=open('EM.npy', 'wb')
f_r=open('R.npy', 'wb')
np.save(f_em, EM)
np.save(f_r, RAD)
f_em.close()
f_r.close()
r = np.load('R.npy').T
em = np.load('EM.npy').T

ax = plt.subplot(1, 2, 1)
ax.set_xlim([0, size])
ax.set_ylim([0, size])
ax.imshow(r,
          vmin=0.01,
          vmax=1)


ax = plt.subplot(1, 2, 2)
ax.set_xlim([0, size])
ax.set_ylim([0, size])
ax.imshow(np.log10(em),
          vmin=7,
          vmax=9)

plt.show()

# plt.savefig('../plot/sed_test.pdf')
