import numpy as np
from scipy import stats as stats
import matplotlib.pyplot as plt
from scipy import fftpack
import pandas as pd
from astropy.io import fits
from skimage import restoration

# https://stackoverflow.com/questions/17473917/is-there-a-equivalent-of-scipy-signal-deconvolve-for-2d-arrays


def convolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft * psf_fft)))


def deconvolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft / psf_fft)))


df = pd.read_csv('../../script/cores006.csv')
corenum = 24

data = fits.open('../sm_abcd_006_lpf_new.fits')[0].data.T
subimage = data[int(df['XPEAK_IMAGE'][corenum]) - 50:int(df['XPEAK_IMAGE'][corenum]) + 50,
                int(df['YPEAK_IMAGE'][corenum]) - 50:int(df['YPEAK_IMAGE'][corenum]) + 50]

subimage = subimage.astype(float)

sx, sy = 100, 100
x2d, y2d = np.ogrid[0:sx, 0:sy]
# beam_ang = np.radians(0)
beam_a = 0.65
beam_b = 0.65
cx, cy = 50, 50

star_ang = np.radians(90-df['THETA_IMAGE'][corenum])
star_a = df['A_IMAGE'][corenum]
star_b = df['A_IMAGE'][corenum]


mask = ((x2d-cx)**2 + (y2d -cy)**2) / (3* np.sqrt(df['A_IMAGE'][corenum] * df['A_IMAGE'][corenum]))**2 > 1

# subimage[mask] = 0


# star = stats.norm.pdf(((((x2d - cx) * np.cos(star_ang) -
#                          (y2d - cy) * np.sin(star_ang)) /
#                         (star_a * 1  / 2.355))**2 +
#                        (((x2d - cx) * np.sin(star_ang) -
#                          (y2d - cy) * np.cos(star_ang)) /
#                         (star_b * 1  / 2.355))**2
#                        ), 0, 1)
# psf = stats.norm.pdf(((((x2d - cx) * np.cos(beam_ang) -
#                         (y2d - cy) * np.sin(beam_ang)) /
#                        (beam_a * 1 / 0.05 / 2.355))**2 +
#                       (((x2d - cx) * np.sin(beam_ang) -
#                         (y2d - cy) * np.cos(beam_ang)) /
#                        (beam_b * 1 / 0.05 / 2.355))**2
#                       ), 0, 1)

psf =  stats.norm.pdf( np.sqrt(((x2d-cx)**2 + (y2d -cy)**2) / (0.65 /0.05/2.355)**2) , 0, 1)
psf_small =  stats.norm.pdf( ((x2d-cx)**2 + (y2d -cy)**2) / (1.85 /0.05/2.355)**2 , 0, 1)

psf = psf / np.nanmax(psf)

star =  stats.norm.pdf( ((x2d-cx)**2 + (y2d -cy)**2) / (np.sqrt(df['A_IMAGE'][corenum] * df['B_IMAGE'][corenum])/2.355)**2 , 0, 1)

# subimage = np.real(psf_small)
star_conv = convolve(subimage, psf)

# star_conv = fftpack.fftshift(star_conv)
# star_deconv = restoration.richardson_lucy(star_conv, psf, iterations=100)

# star_deconv[mask] = 0

star_deconv = deconvolve(subimage, psf)
# star_deconv = fftpack.fftshift(star_deconv)

f, axes = plt.subplots(2, 2)
axes[0, 0].imshow(subimage)
axes[0, 0].set_title('real %i' % (corenum))
axes[0, 1].imshow(psf)
axes[0, 1].set_title('beam %i' % (corenum))
# axes[1, 0].imshow(np.real(star_conv))
axes[1, 0].imshow(np.real(star_conv))
axes[1, 0].set_title('model %i' % (corenum))
axes[1, 1].imshow(np.real(star_deconv))
axes[1, 1].set_title('result %i' % (corenum))
# plt.show()
plt.savefig('conv_test.pdf')
