originalfiles = ['022.4.fits']
rootdir = '../data/'

for file in originalfiles:
    filename = rootdir + file
    importfits(fitsimage=filename,
               imagename=filename.replace('.fits', 'im'))
    filename = filename.replace('.fits', 'im')
    imregrid(imagename=filename,
             template=rootdir + 'rg_096.im',
             output=rootdir + 'rg_' + file.replace('.fits', '.im'))
    exportfits(imagename=rootdir + 'rg_' + file.replace('.fits', '.im'),
               fitsimage=rootdir + 'rg_' + file,
               dropdeg=True)
