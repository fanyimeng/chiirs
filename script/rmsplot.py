from techlib import *
text = {}
# text['cbC_mos_mtmfs_deep.image.tt0'] = ['(a)', 'tt0']
text['rg_006_rms_100_int.fits'] = ['6 GHz (VLA A Array)', 'RMS']
text['dustonly_rms_100_cubic.fits'] = ['Dust Only', 'RMS']
text['rg_096_rms_100_int.fits'] = ['96 GHz (ALMA)', 'RMS']

rootdir = '/Users/meng/Dropbox/astro/chronicle/small/data/'
for i in ['dustonly_rms_100_cubic.fits']:
  techplot_rms(fitsfile=rootdir + '%s' % (i),
               outfile='rms_%s.pdf' % 'dust',
               plottext=[text['dustonly_rms_100_cubic.fits'][0], text['dustonly_rms_100_cubic.fits'][1]],
               start_coor2=["17:47:20.5", "-28:26:29.780"],
               end_coor2=["17:47:20.5", "-28:20:29.780"],
               maplim=[0, 0.8],
               prolim=[0, 0.8],
               fromimage=False)
rootdir = '/Users/meng/Dropbox/astro/chronicle/small/data/'
for i in ['abcd_006_lpf_new_rms_100_1e-3_cubic.fits']:
  techplot_rms(fitsfile=rootdir + '%s' % (i),
               outfile='rms_%s.pdf' % '006',
               plottext=[text['rg_006_rms_100_int.fits'][0], text['rg_006_rms_100_int.fits'][1]],
               start_coor2=["17:47:20.5", "-28:26:29.780"],
               end_coor2=["17:47:20.5", "-28:20:29.780"],
               maplim=[0, 0.8],
               prolim=[0, 0.8],
               fromimage=False)
for i in ['rg_096_rms_100_1e-3_cubic.fits']:
  techplot_rms(fitsfile=rootdir + '%s' % (i),
               outfile='rms_%s.pdf' % '096',
               plottext=[text['rg_096_rms_100_int.fits'][0], text['rg_096_rms_100_int.fits'][1]],
               start_coor2=["17:47:20.5", "-28:26:29.780"],
               end_coor2=["17:47:20.5", "-28:20:29.780"],
               maplim=[0, 0.8],
               prolim=[0, 0.8],
               fromimage=False)
# for i in ['rg_010_rms_100_int.fits']:
#   techplot_rms(fitsfile=rootdir + '%s' % (i),
#                outfile='rms_%s.pdf' % '010',
#                plottext=[text[str(i)][0], text[str(i)][1]],
#                start_coor2=["17:47:20.5", "-28:26:29.780"],
#                end_coor2=["17:47:20.5", "-28:20:29.780"],
#                maplim=[0, 1],
#                prolim=[-0.1, 0.5],
#                fromimage=False)
# for i in ['rg_096_rms_100_int.fits']:
#   techplot_rms(fitsfile=rootdir + '%s' % (i),
#                outfile='rms_%s.pdf' % '096',
#                plottext=[text[str(i)][0], text[str(i)][1]],
#                start_coor2=["17:47:20.5", "-28:26:29.780"],
#                end_coor2=["17:47:20.5", "-28:20:29.780"],
#                maplim=[0, 1],
#                prolim=[-0.1, 0.5],
#                fromimage=False)

# for i in ['cbC_mos_mtmfs_deep.image.tt0']:
#     techplot_sp(fitsfile=rootdir + '%s.fits' % (i),
#                 outfile='mtmfs_%s.pdf' % text[str(i)][0][1],
#                 plottext=[text[str(i)][0], text[str(i)][1]],
#                 maplim=[-1, 25],
#                 prolim=[-3, 5],
#                 fromimage=False)
