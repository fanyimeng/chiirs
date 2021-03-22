from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates
from visualization import make_scalebar
from constants import distance
import matplotlib.patheffects as PathEffects
from matplotlib.colors import PowerNorm
from matplotlib.ticker import *
import numpy as np
from matplotlib import gridspec
import os
from mpl_plot_templates import asinh_norm


def techplot(fitsfile='',
             outfile='',
             plottext='',
             maplim=[-1.0, 25],
             prolim=[-3, 5],
             fromimage=False):

    if fromimage:
        temp_str = '''exportfits(imagename = '%s',
        fitsimage = '%s',
        dropdeg = True,
        overwrite = True)
        '''

        exfits_str = temp_str % (fitsfile, fitsfile + '.fits')
        print(exfits_str)

        f = open('exportfits.py', "w")
        f.write(exfits_str)
        f.close()

        os.system('casa --log2term --nogui -c exportfits.py')

        fitsfile = fitsfile + '.fits'

    def region_anno(ax, physcoor, txtcoor, text):
        sgrb2s_coor = coordinates.SkyCoord(physcoor[0],
                                           physcoor[1],
                                           unit=(u.h, u.deg),
                                           frame='fk5')
        sgrb2s_x, sgrb2s_y = mywcs.wcs_world2pix([[sgrb2s_coor.ra.deg,
                                                   sgrb2s_coor.dec.deg]],
                                                 0)[0]
        sgrb2s_txt_coor = coordinates.SkyCoord(txtcoor[0],
                                               txtcoor[1],
                                               unit=(u.h, u.deg),
                                               frame='fk5')
        sgrb2s_txt_x, sgrb2s_txt_y = mywcs.wcs_world2pix([[sgrb2s_txt_coor.ra.deg,
                                                           sgrb2s_txt_coor.dec.deg]],
                                                         0)[0]
        anno = ax0.annotate(text, xy=(sgrb2s_x, sgrb2s_y),
                            xytext=(sgrb2s_txt_x, sgrb2s_txt_y),
                            color='#FFFFFF',
                            arrowprops=dict(edgecolor='#FFFFFF',
                                            arrowstyle='-', linewidth=0.8),
                            path_effects=[PathEffects.withStroke(linewidth=2,
                                                                 foreground="k")])
        anno.arrow_patch.set_path_effects([
            PathEffects.Stroke(linewidth=2, foreground="k"),
            PathEffects.Normal()])
        return 0

    def getProfileData(ax1, ax2,
                       header,
                       data,
                       start_coor=["17:47:20.5", "-28:23:06"],
                       end_coor=["17:47:24.0", "-28:23:10"],
                       lsty='--',
                       color='b',
                       order=0):
        mywcs = wcs.WCS(header).celestial
        bl, tr = (coordinates.SkyCoord(start_coor[0],
                                       start_coor[1],
                                       unit=(u.h, u.deg),
                                       frame='fk5'),
                  coordinates.SkyCoord(end_coor[0],
                                       end_coor[1],
                                       unit=(u.h, u.deg),
                                       frame='fk5'))
        (zx1, zy1), (zx2, zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                        bl.dec.deg]], 0)[0],
                                  mywcs.wcs_world2pix([[tr.ra.deg,
                                                        tr.dec.deg]], 0)[0])
        vert = np.abs(zy2 - zy1) > np.abs(zx2 - zx1)
        data = data.transpose()
        if vert:
            fluxlist = data[int(zx1), int(zy1):int(zy2)]
            limit = [zy1, zy2]
            ax1.vlines(x=int(zx1), ymin=int(zy1), ymax=int(zy2), lw=2,
                       linestyle=lsty, color=color)
        else:
            fluxlist = data[int(zx1):int(zx2), int(zy1)]
            limit = [zx1, zx2]
            ax1.hlines(y=int(zy1), xmin=int(zx1), xmax=int(zx2), lw=2,
                       linestyle=lsty, color=color)
        ax2.axhline(y=0, ls='--', lw=0.8, color='grey')
        position = range(len(fluxlist))
        position = position * np.abs(header['CDELT2']) * 3600
        position = np.array(position) - np.ones_like(position) * \
            np.nanmean(position)
        fluxlist = fluxlist
        ax2.plot(position, fluxlist, lw='0.9', color=color)
        ax2.fill_between(position, fluxlist, 0,
                         facecolor=color,  # The fill color
                         alpha=0.2,
                         zorder=100)
        return fluxlist, position

    plt.rcParams['figure.dpi'] = 75.
    plt.rcParams['savefig.dpi'] = 300.
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['axes.titlesize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.clf()
    tick_fontsize = 9
    fig = plt.figure(figsize=(4, 16 / 3))
    data_path = fitsfile
    hdu = fits.open(data_path)[0]
    mywcs = wcs.WCS(hdu.header).celestial
    bl, tr = (coordinates.SkyCoord("17:47:33.07874",
                                   "-28:26:29.780",
                                   unit=(u.h, u.deg),
                                   frame='fk5'),
              coordinates.SkyCoord("17:47:05.79726",
                                   "-28:20:29.780",
                                   unit=(u.h, u.deg),
                                   frame='fk5'))
    (zx1, zy1), (zx2, zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                    bl.dec.deg]], 0)[0],
                              mywcs.wcs_world2pix([[tr.ra.deg,
                                                    tr.dec.deg]], 0)[0])
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0], projection=mywcs)
    ax0.set_xlim(zx1, zx2)
    ax0.set_ylim(zy1, zy2)
    im = ax0.imshow(hdu.data.squeeze() * 1e3,
                    transform=ax0.get_transform(mywcs),
                    vmin=maplim[0],
                    vmax=maplim[1],
                    cmap=plt.cm.binary,
                    interpolation='nearest',
                    origin='lower',
                    aspect='equal',
                    norm=PowerNorm(gamma=0.5))
    ra = ax0.coords[0]
    dec = ax0.coords[1]
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm')
    ra.set_axislabel('RA (J2000)', minpad=0.5)
    dec.set_axislabel('Dec (J2000)', minpad=-0.9)
    ra.display_minor_ticks(True)
    ra.set_minor_frequency(6)
    dec.display_minor_ticks(True)
    dec.set_minor_frequency(6)
    ra.set_ticks_position('all')
    ra.set_ticklabel_position('t')
    ra.set_axislabel_position('t')
    dec.set_ticks_position('all')
    dec.set_ticklabel_position('l')
    dec.set_axislabel_position('l')
    ra.set_ticks(size=7, direction='in', color='k')
    dec.set_ticks(size=7, direction='in', color='k')
    ax0.tick_params(which='minor', length=4)

    ax0.text((zx2 - zx1) * 0.95 + zx1,
             (zy2 - zy1) * 0.92 + zy1,
             plottext,
             fontsize=11,
             ha='right',
             color='white',
             path_effects=[PathEffects.withStroke(linewidth=2.0,
                                                  foreground="k")])

    pat_bl, pat_tr = (coordinates.SkyCoord("17:47:25.0",
                                           "-28:26:00",
                                           unit=(u.h, u.deg),
                                           frame='fk5'),
                      coordinates.SkyCoord("17:47:18",
                                           "-28:24:20",
                                           unit=(u.h, u.deg),
                                           frame='fk5'))
    (px1, py1), (px2, py2) = (mywcs.wcs_world2pix([[pat_bl.ra.deg,
                                                    pat_bl.dec.deg]],
                                                  0)[0],
                              mywcs.wcs_world2pix([[pat_tr.ra.deg,
                                                    pat_tr.dec.deg]],
                                                  0)[0])

    region_anno(ax=ax0,
                physcoor=["17:47:20.5", "-28:23:06"],
                txtcoor=["17:47:24.0", "-28:23:10"],
                text='M')
    region_anno(ax=ax0,
                physcoor=["17:47:20.2", "-28:22:21"],
                txtcoor=["17:47:16.5", "-28:21:50"],
                text='N')
    region_anno(ax=ax0,
                physcoor=["17:47:20.430", "-28:23:45.06"],
                txtcoor=["17:47:17.5", "-28:23:55"],
                text='S')
    region_anno(ax=ax0,
                physcoor=["17:47:21.0", "-28:25:10"],
                txtcoor=["17:47:18", "-28:25:10"],
                text='DS')
    region_anno(ax=ax0,
                physcoor=["17:47:19.3", "-28:24:38"],
                txtcoor=["17:47:17", "-28:24:20"],
                text='AA')
    region_anno(ax=ax0,
                physcoor=["17:47:13.0", "-28:24:40"],
                txtcoor=["17:47:10", "-28:24:20"],
                text='V')

    scalebarpos = coordinates.SkyCoord(
        "17:47:11", "-28:26:00", unit=(u.h, u.deg), frame='fk5')
    scalebarlen = (2 * u.pc / distance).to(u.arcsec,
                                           u.dimensionless_angles())
    make_scalebar(ax0, scalebarpos,
                  length=scalebarlen,
                  color='k',
                  label='50" 2 pc',
                  text_offset=1.0 * u.arcsec,
                  fontsize=8
                  )

    ax1 = plt.subplot(gs[1])

    prof = getProfileData(ax1=ax0, ax2=ax1,
                          header=hdu.header,
                          data=hdu.data.squeeze() * 1e3,
                          start_coor=["17:47:33.07874", "-28:25:10.0"],
                          end_coor=["17:47:05.79726", "-28:25:10.0"],
                          lsty='--',
                          color=(0.5, 0.5, 1.0, 1.0),
                          order=0)

    prof = getProfileData(ax1=ax0, ax2=ax1,
                          header=hdu.header,
                          data=hdu.data.squeeze() * 1e3,
                          start_coor=["17:47:26.4", "-28:26:29.780"],
                          end_coor=["17:47:26.4", "-28:20:29.780"],
                          lsty='--',
                          color=(1.0, 0.5, 0.5, 1.0),
                          order=0)

    ax1.set_xlim([np.nanmin(prof[1]), np.nanmax(prof[1])])
    ax1.set_ylim(prolim)
    asp = np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0]
    ax1.set_aspect(asp / 3)

    # ax1.yaxis.set_major_locator(MultipleLocator(2))
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.xaxis.set_major_locator(MultipleLocator(60))
    ax1.xaxis.set_minor_locator(MultipleLocator(10))

    ax1.tick_params(axis='both', which='major', direction='in', size=6,
                    bottom=True, top=True, left=True, right=True)
    ax1.tick_params(axis='both', which='minor', direction='in', size=4,
                    bottom=True, top=True, left=True, right=True)
    ax1.set_xlabel('Offset [arcsec]')
    ax1.set_ylabel("$S$ [mJy beam$^{-1}$]")

    # ax1.grid(True)

    plt.subplots_adjust(hspace=0.05)
    pcbar = ax0.get_position().get_points().flatten()
    cbar_ax = fig.add_axes(
        [pcbar[2] + 0.01 * (pcbar[2] - pcbar[0]),
         pcbar[1],
         0.02,
         pcbar[3] - pcbar[1]])
    cb = plt.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar_ax.yaxis.set_ticks_position('right')
    cbar_ax.yaxis.set_label_position('right')
    cbar_ax.tick_params(axis='both', which='both', direction='in')
    cb.set_label("$S$ [mJy beam$^{-1}$]")

    plt.savefig(outfile, bbox_inches='tight')
    plt.clf()
    return 0


'''

























'''


def techplot_sp(fitsfile='',
                outfile='',
                plottext=['', ''],
                maplim=[-1.0, 25],
                prolim=[-3, 5],
                fromimage=False,
                img_bl=["17:47:33.07874",
                        "-28:26:29.780"],
                img_tr=["17:47:05.79726",
                        "-28:20:29.780"],
                start_coor1=["17:47:33.07874", "-28:25:10.0"],
                end_coor1=["17:47:05.79726", "-28:25:10.0"],
                start_coor2=["17:47:26.4", "-28:26:29.780"],
                end_coor2=["17:47:26.4", "-28:20:29.780"],
                sb_pos=["17:47:11", "-28:26:00"],
                anno_region=True,
                writebeam=True,
                normalizedprofile=False
                ):

    if fromimage:
        temp_str = '''exportfits(imagename = '%s',
        fitsimage = '%s',
        dropdeg = True,
        overwrite = True)
        '''

        exfits_str = temp_str % (fitsfile, fitsfile + '.fits')
        print(exfits_str)

        f = open('exportfits.py', "w")
        f.write(exfits_str)
        f.close()

        os.system('casa --log2term --nogui -c exportfits.py')

        fitsfile = fitsfile + '.fits'

    def region_anno(ax, physcoor, txtcoor, text):
        sgrb2s_coor = coordinates.SkyCoord(physcoor[0],
                                           physcoor[1],
                                           unit=(u.h, u.deg),
                                           frame='fk5')
        sgrb2s_x, sgrb2s_y = mywcs.wcs_world2pix([[sgrb2s_coor.ra.deg,
                                                   sgrb2s_coor.dec.deg]],
                                                 0)[0]
        sgrb2s_txt_coor = coordinates.SkyCoord(txtcoor[0],
                                               txtcoor[1],
                                               unit=(u.h, u.deg),
                                               frame='fk5')
        sgrb2s_txt_x, sgrb2s_txt_y = mywcs.wcs_world2pix([[sgrb2s_txt_coor.ra.deg,
                                                           sgrb2s_txt_coor.dec.deg]],
                                                         0)[0]
        anno = ax0.annotate(text, xy=(sgrb2s_x, sgrb2s_y),
                            xytext=(sgrb2s_txt_x, sgrb2s_txt_y),
                            color='#FFFFFF',
                            arrowprops=dict(edgecolor='#FFFFFF',
                                            arrowstyle='-', linewidth=0.8),
                            path_effects=[PathEffects.withStroke(linewidth=2,
                                                                 foreground="k")])
        anno.arrow_patch.set_path_effects([
            PathEffects.Stroke(linewidth=2, foreground="k"),
            PathEffects.Normal()])
        return 0

    def getProfileData(ax1, ax2,
                       header,
                       data,
                       start_coor=["17:47:20.5", "-28:23:06"],
                       end_coor=["17:47:24.0", "-28:23:10"],
                       lsty='--',
                       color='b',
                       order=0):
        mywcs = wcs.WCS(header).celestial
        bl, tr = (coordinates.SkyCoord(start_coor[0],
                                       start_coor[1],
                                       unit=(u.h, u.deg),
                                       frame='fk5'),
                  coordinates.SkyCoord(end_coor[0],
                                       end_coor[1],
                                       unit=(u.h, u.deg),
                                       frame='fk5'))
        (zx1, zy1), (zx2, zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                        bl.dec.deg]], 0)[0],
                                  mywcs.wcs_world2pix([[tr.ra.deg,
                                                        tr.dec.deg]], 0)[0])
        vert = np.abs(zy2 - zy1) > np.abs(zx2 - zx1)
        data = data.transpose()
        if vert:
            fluxlist = data[int(zx1), int(zy1):int(zy2)]
            limit = [zy1, zy2]
            ax1.vlines(x=int(zx1), ymin=int(zy1), ymax=int(zy2), lw=2,
                       linestyle=lsty, color=color)
        else:
            fluxlist = data[int(zx1):int(zx2), int(zy1)]
            limit = [zx1, zx2]
            ax1.hlines(y=int(zy1), xmin=int(zx1), xmax=int(zx2), lw=2,
                       linestyle=lsty, color=color)
        ax2.axhline(y=0, ls='--', lw=0.8, color='grey')
        position = range(len(fluxlist))
        position = position * np.abs(header['CDELT2']) * 3600
        position = np.array(position) - np.ones_like(position) * \
            np.nanmean(position)
        fluxlist = fluxlist
        ax2.plot(position, fluxlist, lw='0.9', color=color)
        ax2.fill_between(position, fluxlist, 0,
                         facecolor=color,  # The fill color
                         alpha=0.2,
                         zorder=100)
        return fluxlist, position

    plt.rcParams['figure.dpi'] = 75.
    plt.rcParams['savefig.dpi'] = 300.
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['axes.titlesize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.clf()
    tick_fontsize = 9
    fig = plt.figure(figsize=(4, 16 / 3))
    data_path = fitsfile
    hdu = fits.open(data_path)[0]
    mywcs = wcs.WCS(hdu.header).celestial
    bl, tr = (coordinates.SkyCoord(img_bl[0],
                                   img_bl[1],
                                   unit=(u.h, u.deg),
                                   frame='fk5'),
              coordinates.SkyCoord(img_tr[0],
                                   img_tr[1],
                                   unit=(u.h, u.deg),
                                   frame='fk5'))
    (zx1, zy1), (zx2, zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                    bl.dec.deg]], 0)[0],
                              mywcs.wcs_world2pix([[tr.ra.deg,
                                                    tr.dec.deg]], 0)[0])
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0], projection=mywcs)
    ax0.set_xlim(zx1, zx2)
    ax0.set_ylim(zy1, zy2)
    im = ax0.imshow(hdu.data.squeeze() * 1e3,
                    transform=ax0.get_transform(mywcs),
                    vmin=maplim[0],
                    vmax=maplim[1],
                    cmap=plt.cm.binary,
                    interpolation='nearest',
                    origin='lower',
                    aspect='equal',
                    norm=PowerNorm(gamma=0.5))
    ra = ax0.coords[0]
    dec = ax0.coords[1]
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm')
    ra.set_axislabel('RA (J2000)', minpad=0.5)
    dec.set_axislabel('Dec (J2000)', minpad=-0.9)
    ra.display_minor_ticks(True)
    ra.set_minor_frequency(6)
    dec.display_minor_ticks(True)
    dec.set_minor_frequency(6)
    ra.set_ticks_position('all')
    ra.set_ticklabel_position('t')
    ra.set_axislabel_position('t')
    dec.set_ticks_position('all')
    dec.set_ticklabel_position('l')
    dec.set_axislabel_position('l')
    ra.set_ticks(size=7, direction='in', color='k')
    dec.set_ticks(size=7, direction='in', color='k')
    ax0.tick_params(which='minor', length=4)

    ax0.text((zx2 - zx1) * 0.95 + zx1,
             (zy2 - zy1) * 0.92 + zy1,
             plottext[1],
             fontsize=11,
             ha='right',
             color='white',
             path_effects=[PathEffects.withStroke(linewidth=2.0,
                                                  foreground="k")])
    ax0.text((zx2 - zx1) * 0.05 + zx1,
             (zy2 - zy1) * 0.92 + zy1,
             plottext[0],
             fontsize=11,
             ha='left',
             color='white',
             path_effects=[PathEffects.withStroke(linewidth=2.0,
                                                  foreground="k")])

    if writebeam:
        beamstr = '$(%4.2f^{\\prime\\prime},\\ %4.2f^{\\prime\\prime},\\ %4.1f^{\\circ})$' % (hdu.header['BMAJ'] * 3600,
                                                                                              hdu.header['BMIN'] *
                                                                                              3600,
                                                                                              hdu.header['BPA'])
        beamstr = '(%4.2f", %4.2f", %4.1f°)' % (hdu.header['BMAJ'] * 3600,
                                                hdu.header['BMIN'] * 3600,
                                                hdu.header['BPA'])
        ax0.text((zx2 - zx1) * 0.05 + zx1,
                 (zy2 - zy1) * 0.08 + zy1,
                 beamstr,
                 fontsize=9,
                 ha='left',
                 color='k')

    pat_bl, pat_tr = (coordinates.SkyCoord("17:47:25.0",
                                           "-28:26:00",
                                           unit=(u.h, u.deg),
                                           frame='fk5'),
                      coordinates.SkyCoord("17:47:18",
                                           "-28:24:20",
                                           unit=(u.h, u.deg),
                                           frame='fk5'))
    (px1, py1), (px2, py2) = (mywcs.wcs_world2pix([[pat_bl.ra.deg,
                                                    pat_bl.dec.deg]],
                                                  0)[0],
                              mywcs.wcs_world2pix([[pat_tr.ra.deg,
                                                    pat_tr.dec.deg]],
                                                  0)[0])
    if anno_region:
        region_anno(ax=ax0,
                    physcoor=["17:47:20.5", "-28:23:06"],
                    txtcoor=["17:47:24.0", "-28:23:10"],
                    text='M')
        region_anno(ax=ax0,
                    physcoor=["17:47:20.2", "-28:22:21"],
                    txtcoor=["17:47:16.5", "-28:21:50"],
                    text='N')
        region_anno(ax=ax0,
                    physcoor=["17:47:20.430", "-28:23:45.06"],
                    txtcoor=["17:47:17.5", "-28:23:55"],
                    text='S')
        region_anno(ax=ax0,
                    physcoor=["17:47:21.0", "-28:25:10"],
                    txtcoor=["17:47:18", "-28:25:10"],
                    text='DS')
        region_anno(ax=ax0,
                    physcoor=["17:47:19.3", "-28:24:38"],
                    txtcoor=["17:47:17", "-28:24:20"],
                    text='AA')
        region_anno(ax=ax0,
                    physcoor=["17:47:13.0", "-28:24:40"],
                    txtcoor=["17:47:10", "-28:24:20"],
                    text='V')

    scalebarpos = coordinates.SkyCoord(
        sb_pos[0], sb_pos[1], unit=(u.h, u.deg), frame='fk5')
    scalebarlen = (2 * u.pc / distance).to(u.arcsec,
                                           u.dimensionless_angles())
    make_scalebar(ax0, scalebarpos,
                  length=scalebarlen,
                  color='k',
                  label='50" 2 pc',
                  text_offset=1.0 * u.arcsec,
                  fontsize=8
                  )

    ax1 = plt.subplot(gs[1])

    prof = getProfileData(ax1=ax0, ax2=ax1,
                          header=hdu.header,
                          data=hdu.data.squeeze() * 1e3,
                          start_coor=start_coor1,
                          end_coor=end_coor1,
                          lsty='--',
                          color=(0.5, 0.5, 1.0, 1.0),
                          order=0)

    prof2 = getProfileData(ax1=ax0, ax2=ax1,
                           header=hdu.header,
                           data=hdu.data.squeeze() * 1e3,
                           start_coor=start_coor2,
                           end_coor=end_coor2,
                           lsty='--',
                           color=(1.0, 0.5, 0.5, 1.0),
                           order=0)

    ax1.set_xlim([np.nanmin(prof[1]), np.nanmax(prof[1])])
    ax1.set_ylim(prolim)
    asp = np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0]
    ax1.set_aspect(asp / 3)

    # ax1.yaxis.set_major_locator(MultipleLocator(2))
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.xaxis.set_major_locator(MultipleLocator(60))
    ax1.xaxis.set_minor_locator(MultipleLocator(10))

    ax1.tick_params(axis='both', which='major', direction='in', size=6,
                    bottom=True, top=True, left=True, right=True)
    ax1.tick_params(axis='both', which='minor', direction='in', size=4,
                    bottom=True, top=True, left=True, right=True)
    ax1.set_xlabel('Offset [arcsec]')
    ax1.set_ylabel("$S$ [mJy beam$^{-1}$]")

    # if normalizedprofile:
    #     ax1.set_ylim(np.nanmax([prof[0], prof[2]]) * -
    #                  0.24, np.nanmax([prof[0], prof[2]]) * 1.2)
    #     asp = np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0]
    #     ax1.set_aspect(asp / 3)

    # ax1.grid(True)

    plt.subplots_adjust(hspace=0.05)
    pcbar = ax0.get_position().get_points().flatten()
    cbar_ax = fig.add_axes(
        [pcbar[2] + 0.01 * (pcbar[2] - pcbar[0]),
         pcbar[1],
         0.02,
         pcbar[3] - pcbar[1]])
    cb = plt.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar_ax.yaxis.set_ticks_position('right')
    cbar_ax.yaxis.set_label_position('right')
    cbar_ax.tick_params(axis='both', which='both', direction='in')
    cb.set_label("$S$ [mJy beam$^{-1}$]")

    plt.savefig(outfile, bbox_inches='tight')
    plt.clf()
    return 0


'''

















































'''


def techplot_alpha(fitsfile='',
                   outfile='',
                   plottext=['', ''],
                   maplim=[-1.0, 25],
                   prolim=[-3, 5],
                   fromimage=False,
                   img_bl=["17:47:33.07874",
                           "-28:26:29.780"],
                   img_tr=["17:47:05.79726",
                           "-28:20:29.780"],
                   start_coor1=["17:47:33.07874", "-28:25:10.0"],
                   end_coor1=["17:47:05.79726", "-28:25:10.0"],
                   start_coor2=["17:47:26.4", "-28:26:29.780"],
                   end_coor2=["17:47:26.4", "-28:20:29.780"],
                   sb_pos=["17:47:11", "-28:26:00"],
                   anno_region=True,
                   writebeam=True,
                   normalizedprofile=False
                   ):

    if fromimage:
        temp_str = '''exportfits(imagename = '%s',
        fitsimage = '%s',
        dropdeg = True,
        overwrite = True)
        '''

        exfits_str = temp_str % (fitsfile, fitsfile + '.fits')
        print(exfits_str)

        f = open('exportfits.py', "w")
        f.write(exfits_str)
        f.close()

        os.system('casa --log2term --nogui -c exportfits.py')

        fitsfile = fitsfile + '.fits'

    def region_anno(ax, physcoor, txtcoor, text):
        sgrb2s_coor = coordinates.SkyCoord(physcoor[0],
                                           physcoor[1],
                                           unit=(u.h, u.deg),
                                           frame='fk5')
        sgrb2s_x, sgrb2s_y = mywcs.wcs_world2pix([[sgrb2s_coor.ra.deg,
                                                   sgrb2s_coor.dec.deg]],
                                                 0)[0]
        sgrb2s_txt_coor = coordinates.SkyCoord(txtcoor[0],
                                               txtcoor[1],
                                               unit=(u.h, u.deg),
                                               frame='fk5')
        sgrb2s_txt_x, sgrb2s_txt_y = mywcs.wcs_world2pix([[sgrb2s_txt_coor.ra.deg,
                                                           sgrb2s_txt_coor.dec.deg]],
                                                         0)[0]
        anno = ax0.annotate(text, xy=(sgrb2s_x, sgrb2s_y),
                            xytext=(sgrb2s_txt_x, sgrb2s_txt_y),
                            color='#FFFFFF',
                            arrowprops=dict(edgecolor='#FFFFFF',
                                            arrowstyle='-', linewidth=0.8),
                            path_effects=[PathEffects.withStroke(linewidth=2,
                                                                 foreground="k")])
        anno.arrow_patch.set_path_effects([
            PathEffects.Stroke(linewidth=2, foreground="k"),
            PathEffects.Normal()])
        return 0

    def getProfileData(ax1, ax2,
                       header,
                       data,
                       start_coor=["17:47:20.5", "-28:23:06"],
                       end_coor=["17:47:24.0", "-28:23:10"],
                       lsty='--',
                       color='b',
                       order=0):
        mywcs = wcs.WCS(header).celestial
        bl, tr = (coordinates.SkyCoord(start_coor[0],
                                       start_coor[1],
                                       unit=(u.h, u.deg),
                                       frame='fk5'),
                  coordinates.SkyCoord(end_coor[0],
                                       end_coor[1],
                                       unit=(u.h, u.deg),
                                       frame='fk5'))
        (zx1, zy1), (zx2, zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                        bl.dec.deg]], 0)[0],
                                  mywcs.wcs_world2pix([[tr.ra.deg,
                                                        tr.dec.deg]], 0)[0])
        vert = np.abs(zy2 - zy1) > np.abs(zx2 - zx1)
        data = data.transpose()
        if vert:
            fluxlist = data[int(zx1), int(zy1):int(zy2)]
            limit = [zy1, zy2]
            ax1.vlines(x=int(zx1), ymin=int(zy1), ymax=int(zy2), lw=2,
                       linestyle=lsty, color=color)
        else:
            fluxlist = data[int(zx1):int(zx2), int(zy1)]
            limit = [zx1, zx2]
            ax1.hlines(y=int(zy1), xmin=int(zx1), xmax=int(zx2), lw=2,
                       linestyle=lsty, color=color)
        ax2.axhline(y=0, ls='--', lw=0.8, color='grey')
        position = range(len(fluxlist))
        position = position * np.abs(header['CDELT2']) * 3600
        position = np.array(position) - np.ones_like(position) * \
            np.nanmean(position)
        fluxlist = fluxlist
        ax2.plot(position, fluxlist, lw='0.9', color=color)
        ax2.fill_between(position, fluxlist, 0,
                         facecolor=color,  # The fill color
                         alpha=0.2,
                         zorder=100)
        return fluxlist, position

    plt.rcParams['figure.dpi'] = 75.
    plt.rcParams['savefig.dpi'] = 300.
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['axes.titlesize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.clf()
    tick_fontsize = 9
    fig = plt.figure(figsize=(4, 16 / 3))
    data_path = fitsfile
    hdu = fits.open(data_path)[0]
    mywcs = wcs.WCS(hdu.header).celestial
    bl, tr = (coordinates.SkyCoord(img_bl[0],
                                   img_bl[1],
                                   unit=(u.h, u.deg),
                                   frame='fk5'),
              coordinates.SkyCoord(img_tr[0],
                                   img_tr[1],
                                   unit=(u.h, u.deg),
                                   frame='fk5'))
    (zx1, zy1), (zx2, zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                    bl.dec.deg]], 0)[0],
                              mywcs.wcs_world2pix([[tr.ra.deg,
                                                    tr.dec.deg]], 0)[0])
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0], projection=mywcs)
    ax0.set_xlim(zx1, zx2)
    ax0.set_ylim(zy1, zy2)
    im = ax0.imshow(hdu.data.squeeze(),
                    transform=ax0.get_transform(mywcs),
                    vmin=maplim[0],
                    vmax=maplim[1],
                    cmap=plt.cm.jet,
                    interpolation='nearest',
                    origin='lower',
                    aspect='equal',
                    norm=asinh_norm.AsinhNorm())
    ra = ax0.coords[0]
    dec = ax0.coords[1]
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm')
    ra.set_axislabel('RA (J2000)', minpad=0.5)
    dec.set_axislabel('Dec (J2000)', minpad=-0.9)
    ra.display_minor_ticks(True)
    ra.set_minor_frequency(6)
    dec.display_minor_ticks(True)
    dec.set_minor_frequency(6)
    ra.set_ticks_position('all')
    ra.set_ticklabel_position('t')
    ra.set_axislabel_position('t')
    dec.set_ticks_position('all')
    dec.set_ticklabel_position('l')
    dec.set_axislabel_position('l')
    ra.set_ticks(size=7, direction='in', color='k')
    dec.set_ticks(size=7, direction='in', color='k')
    ax0.tick_params(which='minor', length=4)

    ax0.text((zx2 - zx1) * 0.95 + zx1,
             (zy2 - zy1) * 0.92 + zy1,
             plottext[1],
             fontsize=11,
             ha='right',
             color='white',
             path_effects=[PathEffects.withStroke(linewidth=2.0,
                                                  foreground="k")])
    ax0.text((zx2 - zx1) * 0.05 + zx1,
             (zy2 - zy1) * 0.92 + zy1,
             plottext[0],
             fontsize=11,
             ha='left',
             color='white',
             path_effects=[PathEffects.withStroke(linewidth=2.0,
                                                  foreground="k")])

    if writebeam:
        beamstr = '$(%4.2f^{\\prime\\prime},\\ %4.2f^{\\prime\\prime},\\ %4.1f^{\\circ})$' % (hdu.header['BMAJ'] * 3600,
                                                                                              hdu.header['BMIN'] *
                                                                                              3600,
                                                                                              hdu.header['BPA'])
        beamstr = '(%4.2f", %4.2f", %4.1f°)' % (hdu.header['BMAJ'] * 3600,
                                                hdu.header['BMIN'] * 3600,
                                                hdu.header['BPA'])
        ax0.text((zx2 - zx1) * 0.05 + zx1,
                 (zy2 - zy1) * 0.08 + zy1,
                 beamstr,
                 fontsize=9,
                 ha='left',
                 color='k')

    pat_bl, pat_tr = (coordinates.SkyCoord("17:47:25.0",
                                           "-28:26:00",
                                           unit=(u.h, u.deg),
                                           frame='fk5'),
                      coordinates.SkyCoord("17:47:18",
                                           "-28:24:20",
                                           unit=(u.h, u.deg),
                                           frame='fk5'))
    (px1, py1), (px2, py2) = (mywcs.wcs_world2pix([[pat_bl.ra.deg,
                                                    pat_bl.dec.deg]],
                                                  0)[0],
                              mywcs.wcs_world2pix([[pat_tr.ra.deg,
                                                    pat_tr.dec.deg]],
                                                  0)[0])
    if anno_region:
        region_anno(ax=ax0,
                    physcoor=["17:47:20.5", "-28:23:06"],
                    txtcoor=["17:47:24.0", "-28:23:10"],
                    text='M')
        region_anno(ax=ax0,
                    physcoor=["17:47:20.2", "-28:22:21"],
                    txtcoor=["17:47:16.5", "-28:21:50"],
                    text='N')
        region_anno(ax=ax0,
                    physcoor=["17:47:20.430", "-28:23:45.06"],
                    txtcoor=["17:47:17.5", "-28:23:55"],
                    text='S')
        region_anno(ax=ax0,
                    physcoor=["17:47:21.0", "-28:25:10"],
                    txtcoor=["17:47:18", "-28:25:10"],
                    text='DS')
        region_anno(ax=ax0,
                    physcoor=["17:47:19.3", "-28:24:38"],
                    txtcoor=["17:47:17", "-28:24:20"],
                    text='AA')
        region_anno(ax=ax0,
                    physcoor=["17:47:13.0", "-28:24:40"],
                    txtcoor=["17:47:10", "-28:24:20"],
                    text='V')

    scalebarpos = coordinates.SkyCoord(
        sb_pos[0], sb_pos[1], unit=(u.h, u.deg), frame='fk5')
    scalebarlen = (2 * u.pc / distance).to(u.arcsec,
                                           u.dimensionless_angles())
    make_scalebar(ax0, scalebarpos,
                  length=scalebarlen,
                  color='k',
                  label='50" 2 pc',
                  text_offset=1.0 * u.arcsec,
                  fontsize=8
                  )

    ax1 = plt.subplot(gs[1])

    prof = getProfileData(ax1=ax0, ax2=ax1,
                          header=hdu.header,
                          data=hdu.data.squeeze(),
                          start_coor=start_coor1,
                          end_coor=end_coor1,
                          lsty='--',
                          color=(0.5, 0.5, 1.0, 1.0),
                          order=0)

    prof2 = getProfileData(ax1=ax0, ax2=ax1,
                           header=hdu.header,
                           data=hdu.data.squeeze(),
                           start_coor=start_coor2,
                           end_coor=end_coor2,
                           lsty='--',
                           color=(1.0, 0.5, 0.5, 1.0),
                           order=0)

    ax1.set_xlim([np.nanmin(prof[1]), np.nanmax(prof[1])])
    ax1.set_ylim(prolim)
    asp = np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0]
    ax1.set_aspect(asp / 3)

    # ax1.yaxis.set_major_locator(MultipleLocator(2))
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.xaxis.set_major_locator(MultipleLocator(60))
    ax1.xaxis.set_minor_locator(MultipleLocator(10))

    ax1.tick_params(axis='both', which='major', direction='in', size=6,
                    bottom=True, top=True, left=True, right=True)
    ax1.tick_params(axis='both', which='minor', direction='in', size=4,
                    bottom=True, top=True, left=True, right=True)
    ax1.set_xlabel('Offset [arcsec]')
    ax1.set_ylabel("tt1")

    # if normalizedprofile:
    #     ax1.set_ylim(np.nanmax([prof[0], prof[2]]) * -
    #                  0.24, np.nanmax([prof[0], prof[2]]) * 1.2)
    #     asp = np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0]
    #     ax1.set_aspect(asp / 3)

    # ax1.grid(True)

    plt.subplots_adjust(hspace=0.05)
    pcbar = ax0.get_position().get_points().flatten()
    cbar_ax = fig.add_axes(
        [pcbar[2] + 0.01 * (pcbar[2] - pcbar[0]),
         pcbar[1],
         0.02,
         pcbar[3] - pcbar[1]])
    cb = plt.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar_ax.yaxis.set_ticks_position('right')
    cbar_ax.yaxis.set_label_position('right')
    cbar_ax.tick_params(axis='both', which='both', direction='in')
    cb.set_label("tt1")

    plt.savefig(outfile, bbox_inches='tight')
    plt.clf()
    return 0


'''




































































'''


def techplot_line(fitsfile='',
                  outfile='',
                  plottext=['', ''],
                  maplim=[-1.0, 25],
                  prolim=[-3, 5],
                  fromimage=False,
                  img_bl=["17:47:33.07874",
                          "-28:26:29.780"],
                  img_tr=["17:47:05.79726",
                          "-28:20:29.780"],
                  start_coor1=["17:47:33.07874", "-28:25:10.0"],
                  end_coor1=["17:47:05.79726", "-28:25:10.0"],
                  start_coor2=["17:47:26.4", "-28:26:29.780"],
                  end_coor2=["17:47:26.4", "-28:20:29.780"],
                  sb_pos=["17:47:11", "-28:26:00"],
                  anno_region=True,
                  writebeam=True,
                  normalizedprofile=False,
                  velocity = 65
                  ):

    if fromimage:
        temp_str = '''exportfits(imagename = '%s',
        fitsimage = '%s',
        dropdeg = True,
        overwrite = True)
        '''

        exfits_str = temp_str % (fitsfile, fitsfile + '.fits')
        print(exfits_str)

        f = open('exportfits.py', "w")
        f.write(exfits_str)
        f.close()

        os.system('casa --log2term --nogui -c exportfits.py')

        fitsfile = fitsfile + '.fits'

    def region_anno(ax, physcoor, txtcoor, text):
        sgrb2s_coor = coordinates.SkyCoord(physcoor[0],
                                           physcoor[1],
                                           unit=(u.h, u.deg),
                                           frame='fk5')
        sgrb2s_x, sgrb2s_y = mywcs.wcs_world2pix([[sgrb2s_coor.ra.deg,
                                                   sgrb2s_coor.dec.deg]],
                                                 0)[0]
        sgrb2s_txt_coor = coordinates.SkyCoord(txtcoor[0],
                                               txtcoor[1],
                                               unit=(u.h, u.deg),
                                               frame='fk5')
        sgrb2s_txt_x, sgrb2s_txt_y = mywcs.wcs_world2pix([[sgrb2s_txt_coor.ra.deg,
                                                           sgrb2s_txt_coor.dec.deg]],
                                                         0)[0]
        anno = ax0.annotate(text, xy=(sgrb2s_x, sgrb2s_y),
                            xytext=(sgrb2s_txt_x, sgrb2s_txt_y),
                            color='#FFFFFF',
                            arrowprops=dict(edgecolor='#FFFFFF',
                                            arrowstyle='-', linewidth=0.8),
                            path_effects=[PathEffects.withStroke(linewidth=2,
                                                                 foreground="k")])
        anno.arrow_patch.set_path_effects([
            PathEffects.Stroke(linewidth=2, foreground="k"),
            PathEffects.Normal()])
        return 0

    def getProfileData(ax1, ax2,
                       header,
                       data,
                       pix_coor=["17:47:20.5", "-28:23:06"],
                       lsty='--',
                       color='b',
                       order=0):
        mywcs = wcs.WCS(header).celestial
        bl = coordinates.SkyCoord(pix_coor[0],
                                  pix_coor[1],
                                  unit=(u.h, u.deg),
                                  frame='fk5')
        (zx1, zy1) = mywcs.wcs_world2pix([[bl.ra.deg,
                                           bl.dec.deg]], 0)[0]
        data = data.transpose()
        fluxlist = data[int(zx1), int(zy1), :]
        ax2.axhline(y=0, ls='--', lw=0.8, color='grey')
        ax1.plot(zx1, zy1, '+', color = color, path_effects=[PathEffects.withStroke(linewidth=2,
                                                                 foreground="k")])
        position = np.array(range(len(fluxlist)))
        # for i in range(len(position)):
        #     position[i] = i-
        position = position - (np.ones_like(position) * header['CRPIX3'])
        position = position * header['CDELT3']
        position = position + np.ones_like(position) * header['CRVAL3']
        position = position / 1e3
        fluxlist = fluxlist
        ax2.plot(position, fluxlist, lw='0.9', color=color)
        ax2.fill_between(position, fluxlist, 0,
                         facecolor=color,  # The fill color
                         alpha=0.2,
                         zorder=100)
        return fluxlist, position

    plt.rcParams['figure.dpi'] = 75.
    plt.rcParams['savefig.dpi'] = 300.
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['axes.titlesize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.clf()
    tick_fontsize = 9
    fig = plt.figure(figsize=(4, 16 / 3))
    data_path = fitsfile
    hdu = fits.open(data_path)[0]
    header = hdu.header
    channel = velocity*1e3-header['CRVAL3']
    channel = channel/header['CDELT3']
    channel = int(channel-header['CRPIX3'])
    mywcs = wcs.WCS(hdu.header).celestial
    bl, tr = (coordinates.SkyCoord(img_bl[0],
                                   img_bl[1],
                                   unit=(u.h, u.deg),
                                   frame='fk5'),
              coordinates.SkyCoord(img_tr[0],
                                   img_tr[1],
                                   unit=(u.h, u.deg),
                                   frame='fk5'))
    (zx1, zy1), (zx2, zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                    bl.dec.deg]], 0)[0],
                              mywcs.wcs_world2pix([[tr.ra.deg,
                                                    tr.dec.deg]], 0)[0])
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0], projection=mywcs)
    ax0.set_xlim(zx1, zx2)
    ax0.set_ylim(zy1, zy2)
    im = ax0.imshow(hdu.data.squeeze()[channel] * 1e3,
                    transform=ax0.get_transform(mywcs),
                    vmin=maplim[0],
                    vmax=maplim[1],
                    cmap=plt.cm.binary,
                    interpolation='nearest',
                    origin='lower',
                    aspect='equal',
                    norm=PowerNorm(gamma=0.5))
    ra = ax0.coords[0]
    dec = ax0.coords[1]
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm')
    ra.set_axislabel('RA (J2000)', minpad=0.5)
    dec.set_axislabel('Dec (J2000)', minpad=-0.9)
    ra.display_minor_ticks(True)
    ra.set_minor_frequency(6)
    dec.display_minor_ticks(True)
    dec.set_minor_frequency(6)
    ra.set_ticks_position('all')
    ra.set_ticklabel_position('t')
    ra.set_axislabel_position('t')
    dec.set_ticks_position('all')
    dec.set_ticklabel_position('l')
    dec.set_axislabel_position('l')
    ra.set_ticks(size=7, direction='in', color='k')
    dec.set_ticks(size=7, direction='in', color='k')
    ax0.tick_params(which='minor', length=4)

    ax0.text((zx2 - zx1) * 0.95 + zx1,
             (zy2 - zy1) * 0.92 + zy1,
             plottext[1],
             fontsize=11,
             ha='right',
             color='white',
             path_effects=[PathEffects.withStroke(linewidth=2.0,
                                                  foreground="k")])
    ax0.text((zx2 - zx1) * 0.05 + zx1,
             (zy2 - zy1) * 0.92 + zy1,
             plottext[0],
             fontsize=11,
             ha='left',
             color='white',
             path_effects=[PathEffects.withStroke(linewidth=2.0,
                                                  foreground="k")])

    if writebeam:
        beamstr = '$(%4.2f^{\\prime\\prime},\\ %4.2f^{\\prime\\prime},\\ %4.1f^{\\circ})$' % (hdu.header['BMAJ'] * 3600,
                                                                                              hdu.header['BMIN'] *
                                                                                              3600,
                                                                                              hdu.header['BPA'])
        beamstr = '(%4.2f", %4.2f", %4.1f°)' % (hdu.header['BMAJ'] * 3600,
                                                hdu.header['BMIN'] * 3600,
                                                hdu.header['BPA'])
        ax0.text((zx2 - zx1) * 0.05 + zx1,
                 (zy2 - zy1) * 0.08 + zy1,
                 beamstr,
                 fontsize=9,
                 ha='left',
                 color='k')

    pat_bl, pat_tr = (coordinates.SkyCoord("17:47:25.0",
                                           "-28:26:00",
                                           unit=(u.h, u.deg),
                                           frame='fk5'),
                      coordinates.SkyCoord("17:47:18",
                                           "-28:24:20",
                                           unit=(u.h, u.deg),
                                           frame='fk5'))
    (px1, py1), (px2, py2) = (mywcs.wcs_world2pix([[pat_bl.ra.deg,
                                                    pat_bl.dec.deg]],
                                                  0)[0],
                              mywcs.wcs_world2pix([[pat_tr.ra.deg,
                                                    pat_tr.dec.deg]],
                                                  0)[0])
    if anno_region:
        region_anno(ax=ax0,
                    physcoor=["17:47:20.5", "-28:23:06"],
                    txtcoor=["17:47:24.0", "-28:23:10"],
                    text='M')
        region_anno(ax=ax0,
                    physcoor=["17:47:20.2", "-28:22:21"],
                    txtcoor=["17:47:16.5", "-28:21:50"],
                    text='N')
        region_anno(ax=ax0,
                    physcoor=["17:47:20.430", "-28:23:45.06"],
                    txtcoor=["17:47:17.5", "-28:23:55"],
                    text='S')
        region_anno(ax=ax0,
                    physcoor=["17:47:21.0", "-28:25:10"],
                    txtcoor=["17:47:18", "-28:25:10"],
                    text='DS')
        region_anno(ax=ax0,
                    physcoor=["17:47:19.3", "-28:24:38"],
                    txtcoor=["17:47:17", "-28:24:20"],
                    text='AA')
        region_anno(ax=ax0,
                    physcoor=["17:47:13.0", "-28:24:40"],
                    txtcoor=["17:47:10", "-28:24:20"],
                    text='V')

    scalebarpos = coordinates.SkyCoord(
        sb_pos[0], sb_pos[1], unit=(u.h, u.deg), frame='fk5')
    scalebarlen = (2 * u.pc / distance).to(u.arcsec,
                                           u.dimensionless_angles())
    make_scalebar(ax0, scalebarpos,
                  length=scalebarlen,
                  color='k',
                  label='50" 2 pc',
                  text_offset=1.0 * u.arcsec,
                  fontsize=8
                  )

    ax1 = plt.subplot(gs[1])

    prof = getProfileData(ax1=ax0, ax2=ax1,
                          header=hdu.header,
                          data=hdu.data.squeeze() * 1e3,
                          pix_coor=start_coor1,
                          lsty='--',
                          color=(0.5, 0.5, 1.0, 1.0),
                          order=0)

    prof2 = getProfileData(ax1=ax0, ax2=ax1,
                           header=hdu.header,
                           data=hdu.data.squeeze() * 1e3,
                           pix_coor=start_coor2,
                           lsty='--',
                           color=(1.0, 0.5, 0.5, 1.0),
                           order=0)

    ax1.set_xlim([np.nanmin(prof[1]), np.nanmax(prof[1])])
    ax1.set_ylim(prolim)
    asp = np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0]
    ax1.set_aspect(asp / 3)
    ax1.vlines(x = velocity, ymin = prolim[0], ymax = prolim[1], color = (0.5,0.5,0.5,1), zorder = 1000)

    # ax1.yaxis.set_major_locator(MultipleLocator(2))
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.xaxis.set_major_locator(MultipleLocator(60))
    ax1.xaxis.set_minor_locator(MultipleLocator(10))

    ax1.tick_params(axis='both', which='major', direction='in', size=6,
                    bottom=True, top=True, left=True, right=True)
    ax1.tick_params(axis='both', which='minor', direction='in', size=4,
                    bottom=True, top=True, left=True, right=True)
    ax1.set_xlabel('$V_{\\rm LSR}$ [km s$^{-1}$]')
    ax1.set_ylabel("$S$ [mJy beam$^{-1}$]")

    # if normalizedprofile:
    #     ax1.set_ylim(np.nanmax([prof[0], prof[2]]) * -
    #                  0.24, np.nanmax([prof[0], prof[2]]) * 1.2)
    #     asp = np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0]
    #     ax1.set_aspect(asp / 3)

    # ax1.grid(True)

    plt.subplots_adjust(hspace=0.05)
    pcbar = ax0.get_position().get_points().flatten()
    cbar_ax = fig.add_axes(
        [pcbar[2] + 0.01 * (pcbar[2] - pcbar[0]),
         pcbar[1],
         0.02,
         pcbar[3] - pcbar[1]])
    cb = plt.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar_ax.yaxis.set_ticks_position('right')
    cbar_ax.yaxis.set_label_position('right')
    cbar_ax.tick_params(axis='both', which='both', direction='in')
    cb.set_label("$S$ [mJy beam$^{-1}$]")

    plt.savefig(outfile, bbox_inches='tight')
    plt.clf()
    return 0
































def techplot_rms(fitsfile='',
                   outfile='',
                   plottext=['', ''],
                   maplim=[-1.0, 25],
                   prolim=[-3, 5],
                   fromimage=False,
                   img_bl=["17:47:33.07874",
                           "-28:26:29.780"],
                   img_tr=["17:47:05.79726",
                           "-28:20:29.780"],
                   start_coor1=["17:47:33.07874", "-28:25:10.0"],
                   end_coor1=["17:47:05.79726", "-28:25:10.0"],
                   start_coor2=["17:47:26.4", "-28:26:29.780"],
                   end_coor2=["17:47:26.4", "-28:20:29.780"],
                   sb_pos=["17:47:11", "-28:26:00"],
                   anno_region=True,
                   writebeam=True,
                   normalizedprofile=False
                   ):

    if fromimage:
        temp_str = '''exportfits(imagename = '%s',
        fitsimage = '%s',
        dropdeg = True,
        overwrite = True)
        '''

        exfits_str = temp_str % (fitsfile, fitsfile + '.fits')
        print(exfits_str)

        f = open('exportfits.py', "w")
        f.write(exfits_str)
        f.close()

        os.system('casa --log2term --nogui -c exportfits.py')

        fitsfile = fitsfile + '.fits'

    def region_anno(ax, physcoor, txtcoor, text):
        sgrb2s_coor = coordinates.SkyCoord(physcoor[0],
                                           physcoor[1],
                                           unit=(u.h, u.deg),
                                           frame='fk5')
        sgrb2s_x, sgrb2s_y = mywcs.wcs_world2pix([[sgrb2s_coor.ra.deg,
                                                   sgrb2s_coor.dec.deg]],
                                                 0)[0]
        sgrb2s_txt_coor = coordinates.SkyCoord(txtcoor[0],
                                               txtcoor[1],
                                               unit=(u.h, u.deg),
                                               frame='fk5')
        sgrb2s_txt_x, sgrb2s_txt_y = mywcs.wcs_world2pix([[sgrb2s_txt_coor.ra.deg,
                                                           sgrb2s_txt_coor.dec.deg]],
                                                         0)[0]
        anno = ax0.annotate(text, xy=(sgrb2s_x, sgrb2s_y),
                            xytext=(sgrb2s_txt_x, sgrb2s_txt_y),
                            color='#FFFFFF',
                            arrowprops=dict(edgecolor='#FFFFFF',
                                            arrowstyle='-', linewidth=0.8),
                            path_effects=[PathEffects.withStroke(linewidth=2,
                                                                 foreground="k")])
        anno.arrow_patch.set_path_effects([
            PathEffects.Stroke(linewidth=2, foreground="k"),
            PathEffects.Normal()])
        return 0

    def getProfileData(ax1, ax2,
                       header,
                       data,
                       start_coor=["17:47:20.5", "-28:23:06"],
                       end_coor=["17:47:24.0", "-28:23:10"],
                       lsty='--',
                       color='b',
                       order=0):
        mywcs = wcs.WCS(header).celestial
        bl, tr = (coordinates.SkyCoord(start_coor[0],
                                       start_coor[1],
                                       unit=(u.h, u.deg),
                                       frame='fk5'),
                  coordinates.SkyCoord(end_coor[0],
                                       end_coor[1],
                                       unit=(u.h, u.deg),
                                       frame='fk5'))
        (zx1, zy1), (zx2, zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                        bl.dec.deg]], 0)[0],
                                  mywcs.wcs_world2pix([[tr.ra.deg,
                                                        tr.dec.deg]], 0)[0])
        vert = np.abs(zy2 - zy1) > np.abs(zx2 - zx1)
        data = data.transpose()
        if vert:
            fluxlist = data[int(zx1), int(zy1):int(zy2)]
            limit = [zy1, zy2]
            ax1.vlines(x=int(zx1), ymin=int(zy1), ymax=int(zy2), lw=2,
                       linestyle=lsty, color=color)
        else:
            fluxlist = data[int(zx1):int(zx2), int(zy1)]
            limit = [zx1, zx2]
            ax1.hlines(y=int(zy1), xmin=int(zx1), xmax=int(zx2), lw=2,
                       linestyle=lsty, color=color)
        ax2.axhline(y=0, ls='--', lw=0.8, color='grey')
        position = range(len(fluxlist))
        position = position * np.abs(header['CDELT2']) * 3600
        position = np.array(position) - np.ones_like(position) * \
            np.nanmean(position)
        fluxlist = fluxlist
        ax2.plot(position, fluxlist, lw='0.9', color=color)
        ax2.fill_between(position, fluxlist, 0,
                         facecolor=color,  # The fill color
                         alpha=0.2,
                         zorder=100)
        return fluxlist, position

    plt.rcParams['figure.dpi'] = 75.
    plt.rcParams['savefig.dpi'] = 300.
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['axes.titlesize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.clf()
    tick_fontsize = 9
    fig = plt.figure(figsize=(4, 16 / 3))
    data_path = fitsfile
    hdu = fits.open(data_path)[0]
    mywcs = wcs.WCS(hdu.header).celestial
    bl, tr = (coordinates.SkyCoord(img_bl[0],
                                   img_bl[1],
                                   unit=(u.h, u.deg),
                                   frame='fk5'),
              coordinates.SkyCoord(img_tr[0],
                                   img_tr[1],
                                   unit=(u.h, u.deg),
                                   frame='fk5'))
    (zx1, zy1), (zx2, zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                    bl.dec.deg]], 0)[0],
                              mywcs.wcs_world2pix([[tr.ra.deg,
                                                    tr.dec.deg]], 0)[0])
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0], projection=mywcs)
    ax0.set_xlim(zx1, zx2)
    ax0.set_ylim(zy1, zy2)
    im = ax0.imshow(hdu.data.squeeze()*1e3,
                    transform=ax0.get_transform(mywcs),
                    vmin=maplim[0],
                    vmax=maplim[1],
                    cmap=plt.cm.jet,
                    interpolation='nearest',
                    origin='lower',
                    aspect='equal',
                    norm=asinh_norm.AsinhNorm())
    ra = ax0.coords[0]
    dec = ax0.coords[1]
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm')
    ra.set_axislabel('RA (J2000)', minpad=0.5)
    dec.set_axislabel('Dec (J2000)', minpad=-0.9)
    ra.display_minor_ticks(True)
    ra.set_minor_frequency(6)
    dec.display_minor_ticks(True)
    dec.set_minor_frequency(6)
    ra.set_ticks_position('all')
    ra.set_ticklabel_position('t')
    ra.set_axislabel_position('t')
    dec.set_ticks_position('all')
    dec.set_ticklabel_position('l')
    dec.set_axislabel_position('l')
    ra.set_ticks(size=7, direction='in', color='k')
    dec.set_ticks(size=7, direction='in', color='k')
    ax0.tick_params(which='minor', length=4)

    ax0.text((zx2 - zx1) * 0.95 + zx1,
             (zy2 - zy1) * 0.92 + zy1,
             plottext[1],
             fontsize=11,
             ha='right',
             color='white',
             path_effects=[PathEffects.withStroke(linewidth=2.0,
                                                  foreground="k")])
    ax0.text((zx2 - zx1) * 0.05 + zx1,
             (zy2 - zy1) * 0.92 + zy1,
             plottext[0],
             fontsize=11,
             ha='left',
             color='white',
             path_effects=[PathEffects.withStroke(linewidth=2.0,
                                                  foreground="k")])

    if 0:
        beamstr = '$(%4.2f^{\\prime\\prime},\\ %4.2f^{\\prime\\prime},\\ %4.1f^{\\circ})$' % (hdu.header['BMAJ'] * 3600,
                                                                                              hdu.header['BMIN'] *
                                                                                              3600,
                                                                                              hdu.header['BPA'])
        beamstr = '(%4.2f", %4.2f", %4.1f°)' % (hdu.header['BMAJ'] * 3600,
                                                hdu.header['BMIN'] * 3600,
                                                hdu.header['BPA'])
        ax0.text((zx2 - zx1) * 0.05 + zx1,
                 (zy2 - zy1) * 0.08 + zy1,
                 beamstr,
                 fontsize=9,
                 ha='left',
                 color='k')

    pat_bl, pat_tr = (coordinates.SkyCoord("17:47:25.0",
                                           "-28:26:00",
                                           unit=(u.h, u.deg),
                                           frame='fk5'),
                      coordinates.SkyCoord("17:47:18",
                                           "-28:24:20",
                                           unit=(u.h, u.deg),
                                           frame='fk5'))
    (px1, py1), (px2, py2) = (mywcs.wcs_world2pix([[pat_bl.ra.deg,
                                                    pat_bl.dec.deg]],
                                                  0)[0],
                              mywcs.wcs_world2pix([[pat_tr.ra.deg,
                                                    pat_tr.dec.deg]],
                                                  0)[0])
    if anno_region:
        region_anno(ax=ax0,
                    physcoor=["17:47:20.5", "-28:23:06"],
                    txtcoor=["17:47:24.0", "-28:23:10"],
                    text='M')
        region_anno(ax=ax0,
                    physcoor=["17:47:20.2", "-28:22:21"],
                    txtcoor=["17:47:16.5", "-28:21:50"],
                    text='N')
        region_anno(ax=ax0,
                    physcoor=["17:47:20.430", "-28:23:45.06"],
                    txtcoor=["17:47:17.5", "-28:23:55"],
                    text='S')
        region_anno(ax=ax0,
                    physcoor=["17:47:21.0", "-28:25:10"],
                    txtcoor=["17:47:18", "-28:25:10"],
                    text='DS')
        region_anno(ax=ax0,
                    physcoor=["17:47:19.3", "-28:24:38"],
                    txtcoor=["17:47:17", "-28:24:20"],
                    text='AA')
        region_anno(ax=ax0,
                    physcoor=["17:47:13.0", "-28:24:40"],
                    txtcoor=["17:47:10", "-28:24:20"],
                    text='V')

    scalebarpos = coordinates.SkyCoord(
        sb_pos[0], sb_pos[1], unit=(u.h, u.deg), frame='fk5')
    scalebarlen = (2 * u.pc / distance).to(u.arcsec,
                                           u.dimensionless_angles())
    make_scalebar(ax0, scalebarpos,
                  length=scalebarlen,
                  color='w',
                  label='50" 2 pc',
                  text_offset=1.0 * u.arcsec,
                  fontsize=8
                  )

    ax1 = plt.subplot(gs[1])

    prof = getProfileData(ax1=ax0, ax2=ax1,
                          header=hdu.header,
                          data=hdu.data.squeeze()*1e3,
                          start_coor=start_coor1,
                          end_coor=end_coor1,
                          lsty='--',
                          color=(0.5, 0.5, 1.0, 1.0),
                          order=0)

    prof2 = getProfileData(ax1=ax0, ax2=ax1,
                           header=hdu.header,
                           data=hdu.data.squeeze()*1e3,
                           start_coor=start_coor2,
                           end_coor=end_coor2,
                           lsty='--',
                           color=(1.0, 0.5, 0.5, 1.0),
                           order=0)

    ax1.set_xlim([np.nanmin(prof[1]), np.nanmax(prof[1])])
    ax1.set_ylim(prolim)
    asp = np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0]
    ax1.set_aspect(asp / 3)

    # ax1.yaxis.set_major_locator(MultipleLocator(2))
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.xaxis.set_major_locator(MultipleLocator(60))
    ax1.xaxis.set_minor_locator(MultipleLocator(10))

    ax1.tick_params(axis='both', which='major', direction='in', size=6,
                    bottom=True, top=True, left=True, right=True)
    ax1.tick_params(axis='both', which='minor', direction='in', size=4,
                    bottom=True, top=True, left=True, right=True)
    ax1.set_xlabel('Offset [arcsec]')
    ax1.set_ylabel("RMS [mJy beam$^{-1}$]")

    # if normalizedprofile:
    #     ax1.set_ylim(np.nanmax([prof[0], prof[2]]) * -
    #                  0.24, np.nanmax([prof[0], prof[2]]) * 1.2)
    #     asp = np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0]
    #     ax1.set_aspect(asp / 3)

    # ax1.grid(True)

    plt.subplots_adjust(hspace=0.05)
    pcbar = ax0.get_position().get_points().flatten()
    cbar_ax = fig.add_axes(
        [pcbar[2] + 0.01 * (pcbar[2] - pcbar[0]),
         pcbar[1],
         0.02,
         pcbar[3] - pcbar[1]])
    cb = plt.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar_ax.yaxis.set_ticks_position('right')
    cbar_ax.yaxis.set_label_position('right')
    cbar_ax.tick_params(axis='both', which='both', direction='in')
    cb.set_label("RMS [mJy beam$^{-1}$]")

    plt.savefig(outfile, bbox_inches='tight')
    plt.clf()
    return 0




