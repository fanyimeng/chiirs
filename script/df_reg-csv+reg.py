import pandas as pd
import numpy as np


def df2reg(input_df,
           output_reg,
           cover_old='x',
           col_list=['XPEAK_IMAGE',
                     'YPEAK_IMAGE',
                     'A_IMAGE',
                     'B_IMAGE',
                     'THETA_IMAGE'],
           shape='ellipse'):
    header_text = '%s %s %s %s %s %s' % (
        '# Region file format: DS9 version 4.1\n',
        'global color=green dashlist=8 3 width=2',
        'font="helvetica 10 normal roman"',
        'select=1 highlite=1 dash=0 fixed=0 edit=1',
        'move=1 delete=1 include=1 source=1\n',
        'image\n')

    regfile = open('%s' % (output_reg), cover_old)
    regfile.write(header_text)
    if shape == 'ellipse':
        reg_fmt = 'ellipse(%7.1f,%7.1f,%7.1f,%7.1f,%+7.1f)\n'
        for i in range(input_df.shape[0]):
            row = input_df.iloc[i, :]
            ell_para_list = [row[col_list[i]] for i in range(len(col_list))]
            print(ell_para_list)
            regfile.write(reg_fmt % tuple(ell_para_list))
    if shape == 'circle':
        reg_fmt = 'circle(%7.1f,%7.1f,%7.1f)\n'
        for i in range(input_df.shape[0]):
            row = input_df.iloc[i, :]
            circ_para_list = [row[col_list[i]] for i in range(4)]
            print(circ_para_list)
            regfile.write(reg_fmt % (circ_para_list[0],
                                     circ_para_list[1],
                                     circ_para_list[3]))
    regfile.close()
    return 0


def reg2df(regionfile, circ=True):
    length = 0
    f = open(regionfile, 'r')
    para_dic = {
        'x_img': [],
        'y_img': [],
        'r_img': []
    }
    for line in f:
        if 'circ' in line:
            parastr = line[line.index('(') + 1:line.index(')')]
            columns = np.array(parastr.split(','))
            columns = columns.astype(np.float)
            print(columns)
            para_dic['x_img'].append(columns[0])
            para_dic['y_img'].append(columns[1])
            para_dic['r_img'].append(np.ceil(columns[2]))
    df = pd.DataFrame(para_dic)
    return df


# df = pd.read_csv('gaume1995.csv')
# df2reg(input_df=df,
#        output_reg='gaume_circ.reg',
#        col_list=['G95_X', 'G95_Y', 'G95_A', 'G95_B', 'G95_T'],
#        cover_old='w',
#        shape='circle')

# df = pd.read_csv('cores006_w096.csv')
# df2reg(input_df=df,
#        output_reg='meng_circ.reg',
#        cover_old='w',
#        shape='circle')


df = reg2df('/Users/meng/Downloads/reg/adjusted_006_2.reg')
df = df.sort_values('y_img')
df.to_csv('006.csv', index=True)
df = reg2df('/Users/meng/Downloads/reg/adjusted_022.reg')
df = df.sort_values('y_img')
df.to_csv('022.csv', index=True)
df = reg2df('/Users/meng/Downloads/reg/adjusted_096.reg')
df = df.sort_values('y_img')
df.to_csv('096.csv', index=True)

