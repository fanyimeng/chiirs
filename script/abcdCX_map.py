from techlib import *

techplot_sp(fitsfile='../data/abcd_006_lpf_new.fits',
            outfile='../plot/abcdC_map.pdf',
            plottext=['C Band', ''],
            maplim=[-0.1, 2],
            prolim=[-0.3, 0.8],
            fromimage=False)
techplot_sp(fitsfile='../data/abcd_010_lpf_new.fits',
            outfile='../plot/abcdX_map.pdf',
            plottext=['X Band', ''],
            maplim=[-0.1, 2],
            prolim=[-0.3, 0.8],
            fromimage=False)
