originalfiles = ['abcd_006_lpf_new.fits',
                 'abcd_010_lpf_new.fits']

rootdir = '../data/'
for file in originalfiles:
    imsmooth(imagename='%s%s' % (rootdir, file),
             major='0.65arcsec',
             minor='0.65arcsec',
             pa='0deg',
             targetres=True,
             outfile='%ssm_%s.im' % (rootdir, file))
    exportfits(imagename='%ssm_%s.im' % (rootdir, file),
               fitsimage='%ssm_%s' % (rootdir, file),
               dropdeg=True)
