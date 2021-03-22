# chiirs
**C**ompact **H**ii **R**egions in **S**grB2


### Review of RMS map

- Original data.
    - 6 GHz: `rg_006.fits` and `dcbC_sc_1.8.image.tt0.pbcor.fits`, [feather](./script/feather.py) into `abcd_006_lpf_new.fits`. 
    - 96 GHz: `rm_096.fits`.
- In [./script/rms_map.py](rms_map.py). Following masking, using `ndimage.generic_filter(im, rmsfunc, size=(100, 100))` to apply a rms filter to the image. The `rmsfunc` is `np.std(x)`.
    - Output: `*_rms_100_1e-3.fits`. `1e-3` is the masking level.
- To fill in the masked-out-pixels, use [fill_sample.py](./script/fill_sample.py). Firstly, the full-resolution images (9600^2 px) here should be resampled into 1/100^2 resolution. Then utilize `scipy.interpolate.griddata((x1, y1),...,method='cubic')` to interpolate cubically to fill in the masked hollow regions in the images. 
    - Output: `*_rms_100_1e-3_cubic.fits`
- Plot the rmsmap using [rmsplot.py](./script/rmsplot.py) `<--` [techlib.py](./script/techlib.py).
	- Output: [rms_006.pdf](./plot/rms_006.pdf) and [rms_096.pdf](./plot/rms_096.pdf)


### Ellipses collision

- write a package to detect collision of ellipses
- apply the ell_col to sextractor to give a mutual-mapping of 6 GHz and 96 GHz cores.

