# importfits(
#     fitsimage='dcbC_sc_1.8.image.tt0.pbcor.fits',
#     imagename='dcbC.im')
# importfits(
#     fitsimage='rg_010.fits',
#     imagename='aX.im')
feather(
    imagename='abcd_006_Nolpf_sd0.1.im',
    highres='rg_006.im',
    lowres='dcbC.im',
    lowpassfiltersd=False,
    sdfactor = 0.2)
exportfits(
    imagename='abcd_006_Nolpf_sd0.1.im',
    fitsimage='abcd_006_Nolpf_sd0.1.fits',
    dropdeg=True)
