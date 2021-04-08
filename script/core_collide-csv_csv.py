import pandas as pd
import numpy as np


def collide_dist(para1, para2):
    dist = np.array(para1[:2]) - np.array(para2[:2])
    dist = np.sqrt(np.sum(dist**2))
    normalized_dist = dist / abs(para1[2] + para2[2])
    return normalized_dist


def core_collide(df1, df2, colname):
    df1[colname] = range(df1.shape[0])
    for i in range(df1.shape[0]):
        dist = 2
        row1 = df1.iloc[i, :]
        para1 = [row1['x_img'], row1['y_img'], row1['r_img']]
        for j in range(df2.shape[0]):
            row2 = df2.iloc[j, :]
            para2 = [row2['x_img'], row2['y_img'], row2['r_img']]
            newdist = collide_dist(para1, para2)
            if newdist < dist:
                dist = newdist
                df1[colname][i] = j
                # print(dist)
        if dist > 1:
            df1[colname][i] = -1
    return df1


df_dic = {}
df_dic['df006'] = pd.read_csv('006.csv')
df_dic['df022'] = pd.read_csv('022.csv')
df_dic['df096'] = pd.read_csv('096.csv')

for band1 in ['006',
              '022',
              '096']:
    for band2 in ['006',
                  '022',
                  '096']:
        df_dic['df%s' % band1] = core_collide(
            df_dic['df%s' % band1],
            df_dic['df%s' % band2],
            '%s_idx' % band2)
    # df = df.sort_values('y_img')
    df_dic['df%s' % band1].to_csv('%s.csv' % band1, index=False)
