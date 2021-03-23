from techlib import *
text = {}
text['6'] = ['(a)', '${\\rm psfmode=clark}$']
text['7'] = ['(b)', '${\\rm psfmode=hogbom}$']

for i in [6, 7]:
    techplot_sp(fitsfile='/Users/meng/Dropbox/astro/images/tech/cnb/xband/0630_2.%d.image.fits' % (i),
                outfile='psfmode_%s.pdf' % text[str(i)][0][1],
                plottext=[text[str(i)][0], text[str(i)][1]],
                maplim=[-2, 50],
                prolim=[-20, 30],
                fromimage=False)

# for i in [8, 9]:
#     techplot_sp(fitsfile='/Users/meng/Dropbox/astro/images/tech/cnb/xband/0630_2.%d.image.fits' % (i),
#                 outfile='gain_%s.pdf' % text[str(i)][0][1],
#                 plottext=[text[str(i)][0], text[str(i)][1]],
#                 maplim=[-2, 50],
#                 prolim=[-20, 30],
#                 fromimage=False)
