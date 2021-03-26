import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
import pandas as pd
from tqdm import tqdm


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

    # collide_degree = 0
    ang_list = np.arange(0, 360, 1)
    x, y = ell_rad(ell2, ang_list)
    ell_in_list = ell_in(ell1, x, y)
    # print(ell_in_list)
    collide_degree = (ell_in_list < 1).sum()
    # print(collide_degree)
    # for i in range(36):
    #     x, y = ell_rad(ell2, i * 10.)
    #     if ell_in(ell1, x, y) < 1:
    #         collide_degree = collide_degree + 1
    return collide_degree


def collidelist(coresa, coresb, outputcsv='cores_006_coll.csv'):
    collist = []
    for i in tqdm(range(len(coresa['cenx'].to_numpy()))):
        colldeg = 0
        coll_obj = -1
        for j in range(len(coresb['cenx'].to_numpy())):
            newcolldeg = ell_collide(coresa.iloc[i].to_numpy()[:5],
                                     coresb.iloc[j].to_numpy()[:5])
            if newcolldeg > colldeg:
                # print(newcolldeg)
                coll_obj = j
        collist.append(coll_obj)
    collist = np.array(collist)
    coresa['collide'] = collist
    coresa.to_csv(outputcsv)
    print((collist >= 0).sum(), '/', collist.shape)


headerlist = ['cenx', 'ceny', 'ella', 'ellb',
              'ellt', 'ellr', 'flx_iso', 'flx_err']
cores006 = pd.read_csv('006_sex_test.txt',
                       delim_whitespace=True,
                       header=None,
                       names=headerlist)
cores096 = pd.read_csv('096_sex_test.txt',
                       delim_whitespace=True,
                       header=None,
                       names=headerlist)

collidelist(cores006, cores096, 'cores_006_coll.csv')
collidelist(cores096, cores006, 'cores_096_coll.csv')


'''
fig = plt.figure(figsize=(8.0, 4.0))
ax = plt.subplot(1, 1, 1)
ax.set_aspect('equal')
ax.set_xlim(0, 9600)
ax.set_ylim(0, 9600)
paraset = []
for i in range(100):
    para = np.random.rand(5)
    scaled = np.array([100, 100, 4, 4, 360])
    paraset.append(scaled * para)
# paraset = np.array(paraset)
collist = np.zeros(len(paraset))
print(collist.shape)
for i in range(len(paraset)):
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
fig.suptitle('ellipse test', fontsize=11)
plt.savefig('ell_collide_test.pdf', bbox_inches='tight')
plt.clf()
plt.close()
'''
