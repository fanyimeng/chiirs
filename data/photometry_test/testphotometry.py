from astropy.io import fits
import numpy as np


def corePhoto_round(pixlist, header):
    flux = []
    factor = 1.133 * 0.65**2 / (0.05**2)
    for pix in pixlist:
        flux.append(pix[2] / factor * 1e3)
    return np.sum(flux), np.nanmax(flux) * factor


def circPixList(center, rad):
    pixlist = []
    x1 = int(center[0]) - int(rad) - 1
    x2 = int(center[0]) + int(rad) + 1
    y1 = int(center[1]) - int(rad) - 1
    y2 = int(center[1]) + int(rad) + 1
    for i in range(x1, x2):
        for j in range(y1, y2):
            if (i - center[0])**2 + (j - center[1])**2 <= rad**2:
                pixlist.append([i, j])
    if rad < 1:
        pixlist = [center]
    return pixlist


radlist = np.arange(0, 2.0, 0.01)
output = []
# for rad in radlist:
#     file = 'sm_photomodel_%.2f.fits' % rad
#     hdu = fits.open(file % rad)[0]
#     data = hdu.data
#     header = hdu.header
#     rad_inpix = 0.65 / 2.355 * 3. / 3600 / header['CDELT2']
#     factor = 1.133 * 0.65**2 / (0.005**2)
#     # print(rad_inpix, rad)
#     sum = []
#     for pix in circPixList([1000, 1000], rad_inpix):
#         sum.append(data[pix[0], pix[1]])
#     # print(np.nanmax(sum), np.nansum(sum) / factor,
#         # np.nansum(sum) / factor / np.nanmax(sum))
#     print([rad,
#                    np.nansum(sum) / factor / np.nanmax(sum),
#                    np.nansum(sum) / factor,
#                    np.nanmax(sum)])
#     print('%2.5f %2.5f' % (rad, np.nansum(sum) / factor / np.nanmax(sum)))
# # np.save('photo_test_result.npy', np.array(output))

output = []
for rad in radlist:
    file = 'photomodel_%.2f.fits' % rad
    hdu = fits.open(file % rad)[0]
    data = hdu.data
    header = hdu.header
    rad_inpix = rad*2 / 2.355 * 3. / 3600 / header['CDELT2']
    factor = 1.133 * 0.005**2 / (0.005**2)
    sum = data
    # print(rad_inpix, rad)
    # sum = []
    # for pix in circPixList([1000, 1000], rad_inpix):
    #     sum.append(data[pix[0], pix[1]])
    # # print(np.nanmax(sum), np.nansum(sum) / factor,
    #     # np.nansum(sum) / factor / np.nanmax(sum))
    output.append([rad,
                   np.nansum(sum) / factor / np.nanmax(sum),
                   np.nansum(sum) / factor,
                   np.nanmax(sum)])
    print('%2.5f %2.5f' % (rad, np.nansum(sum) / factor / np.nanmax(sum)))
    print([rad,
                   np.nansum(sum) / factor / np.nanmax(sum),
                   np.nansum(sum) / factor,
                   np.nanmax(sum)])
# np.save('photo_test_result_abs.npy', np.array(output))