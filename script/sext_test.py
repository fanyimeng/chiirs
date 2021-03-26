
# Detecting the compact sources in NGC 6334-I and I(N) automatically using Source Extractor

import matplotlib
import matplotlib.pylab as plt
# import pyfits
from astropy.io import fits
import astropy.io.ascii as ascii
import os
import numpy as np
from astropy import units as u

matplotlib.get_backend()

import aplpy


def sext_find(filename, rmsfile, backgroundfile, outputname):
    filename = filename
    header = fits.open(filename)[0].header

    fig = plt.figure()
    gc3 = aplpy.FITSFigure(filename, figure=fig, subplot=(1, 1, 1))
    gc3.recenter(x=266.8340, y=-28.3850, radius=3. / 60.)

    gc3.show_colorscale(cmap='gray', vmin=0.0, vmax=0.01, stretch='asinh')
    gc3.add_colorbar()
    gc3.colorbar.set_axis_label_text('Intensity (Jy/beam)')

    import sewpy

    sew = sewpy.SEW(params=["X_IMAGE", "Y_IMAGE", "FLUX_RADIUS", "FLUX_MAX", "FLUX_ISO", "FLUXERR_ISO", "FLAGS",
                            "A_IMAGE",
                            "B_IMAGE",
                            "THETA_IMAGE"],
                    config={"DETECT_MINAREA": 60,
                            "DETECT_MAXAREA": 2500 * 4,
                            "DETECT_THRESH": 5,
                            "BACK_SIZE": 32,
                            # "BACK_FILTERSIZE": 32,
                            "BACKPHOTO_TYPE": 'LOCAL',
                            "BACKPHOTO_THICK": 24,
                            "THRESH_TYPE": "RELATIVE",
                            "DEBLEND_MINCONT": 0.00001,
                            "DEBLEND_NTHRESH": 64,
                            # "CHECKIMAGE_NAME": outputname + '_check_rms.fits',
                            # "CHECKIMAGE_TYPE": "BACKGROUND"
                            })
    sizefactor = 2.0
    fluxfactor = 4e3 * np.log(2.0) * header['CDELT2'] * header['CDELT1'] / (
        3.14159 * header['BMAJ'] * header['BMIN'])
    # fluxfactor = 1e3
    fluxfactor = np.abs(fluxfactor)
    out = sew(filename)
    ot = out["table"]  # this is an astropy table.
    x = ot["X_IMAGE"]
    y = ot["Y_IMAGE"]
    a = ot["A_IMAGE"]
    b = ot["B_IMAGE"]
    theta = ot["THETA_IMAGE"]
    r = ot["FLUX_RADIUS"]
    flx_max = ot["FLUX_MAX"]
    flx_iso = ot["FLUX_ISO"]
    flx_err = ot["FLUXERR_ISO"]
    ra, dec = gc3.pixel2world(x, y)
    gc3.show_ellipses(x, y, width=a * sizefactor, height=b * sizefactor, angle=theta,
                      lw=0.05, facecolor='none', edgecolor='blue', coords_frame='pixel')

    print(len(r))
    # rmshdu = fits.open('dustonly_rms_100_cubic.fits')
    rmshdu = fits.open(rmsfile)
    rms = rmshdu[0].data
    bgdhdu = fits.open(backgroundfile)
    bgd = bgdhdu[0].data
    newx = []
    newy = []
    newa = []
    newb = []
    newt = []
    newr = []
    newflx_max = []
    newflx_iso = []
    newflx_err = []
    for i in range(len(x)):
        eccen = a[i] < 3 * b[i]
        size = a[i] < 50 / 2.
        aboverms = flx_max[i] > 10. * rms[int(x[i]), int(y[i])]
        not_v = not (6230 < x[i] < 7000 and 2500 < y[i] < 3227)
        if eccen and size and aboverms and not_v:
            newx.append(x[i])
            newy.append(y[i])
            newa.append(a[i])
            newb.append(b[i])
            newt.append(theta[i])
            newr.append(r[i])
            newflx_max.append(flx_max[i])
            newflx_iso.append(flx_iso[i])
            newflx_err.append(flx_err[i])
    print(len(newx))
    gc3.add_label(0.8, 0.92, str(len(newx)) + ' ' + outputname.replace('_sex', ''),
                  relative='True', color='gray', size='small')

    newx = np.array(newx)
    newy = np.array(newy)
    newa = np.array(newa)
    newb = np.array(newb)
    newt = np.array(newt)
    newr = np.array(newr)
    newflx_max = np.array(newflx_max)
    newflx_iso = np.array(newflx_iso)
    newflx_err = np.array(newflx_err)
    gc3.show_ellipses(newx, newy, width=newa * sizefactor, height=newb * sizefactor, angle=newt,
                      lw=0.1, facecolor='none', edgecolor=(1, 0.5, 0.5, 1), coords_frame='pixel')
    fig.savefig(outputname + '.pdf', dpi=900)
    np.savetxt(outputname + '.txt',
               np.c_[newx, newy, newa, newb, newt, newr, newflx_iso * fluxfactor, newflx_err * fluxfactor])
    return 0


# sext_find('../data/abcd_006_lpf_new.fits',
#           '../data/abcd_006_lpf_new_rms_100_1e-3_cubic.fits', '006_sex_test')
sext_find('../data/abcd_006_lpf_new.fits',
          '../data/abcd_006_lpf_new_rms_100_1e-3_cubic.fits',
          '../data/006_sex_test_check.fits',
          '006_sex_test')
sext_find('../data/rg_096.fits',
          '../data/rg_096_rms_100_1e-3_cubic.fits',
          '../data/096_sex_test_check.fits',
          '096_sex_test')
