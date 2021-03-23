import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from scipy.optimize import fsolve


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
    for i in range(36):
        x, y = ell_rad(ell2, i * 10.)
        if ell_in(ell1, x, y) < 1:
            return 1
    return 0


fig = plt.figure(figsize=(8.0, 4.0))
ax = plt.subplot(1, 1, 1)
ax.set_aspect('equal')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
paraset = []
for i in range(100):
    para = np.random.rand(5)
    scaled = np.array([100, 100, 4, 4, 360])
    paraset.append(scaled * para)
# paraset = np.array(paraset)
collist = np.zeros(len(paraset))
print(collist.shape)
for i in range(len(paraset)):
    # pointx = []
    # pointy = []
    # for ang in range(12):
    #     pointx.append(ell_rad(paraset[i], ang*30.)[0])
    #     pointy.append(ell_rad(paraset[i], ang*30.)[1])
    # ax.scatter(pointx,pointy)

    coll = 0
    for j in range(i + 1, len(paraset)):
        if(ell_collide(paraset[i], paraset[j])):
            collist[i] = 1
            collist[j] = 1
for i in range(len(collist)):
    e1 = patches.Ellipse((paraset[i][0], paraset[i][1]),
                         paraset[i][2] * 2, paraset[i][3] * 2,
                         angle=paraset[i][4],
                         linewidth=0.5,
                         fill=False,
                         zorder=2,
                         edgecolor=(collist[i], 0, 1 - collist[i]))
    ax.add_patch(e1)
# ax.axvline(x=ell_rad([30, 20, 300]) * np.cos(np.radians(120 + 300)) + 20)
# ax.axhline(y=ell_rad([30, 20, 300]) * np.sin(np.radians(120 + 300)) + 30)
fig.suptitle('ellipse test', fontsize=11)
plt.savefig('ell_collide_test.pdf', bbox_inches='tight')
plt.clf()
plt.close()
