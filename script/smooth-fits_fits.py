originalfiles = ['abcd_006_lpf_new.fits']

rootdir = '../data/'
# for file in originalfiles:
#     imsmooth(imagename=rootdir + file,
#              major='0.65arcsec',
#              minor='0.65arcsec',
#              pa='0deg',
#              targetres=True,
#              outfile=rootdir + 'sm_' + file.replace('fits', 'im'))
#     exportfits(imagename=rootdir + 'sm_' + file.replace('fits', 'im'),
#                fitsimage=rootdir + 'sm_' + file,
#                dropdeg=True)
rootdir = '../data/'
for file in originalfiles:
    imsmooth(imagename=rootdir + file,
             major='2.724arcsec',
             minor='2.523arcsec',
             pa='-81.685deg',
             targetres=True,
             outfile=rootdir + 'DCsm_' + file.replace('fits', 'im'))
    exportfits(imagename=rootdir + 'DCsm_' + file.replace('fits', 'im'),
               fitsimage=rootdir + 'DCsm_' + file,
               dropdeg=True)
