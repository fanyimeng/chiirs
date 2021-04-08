from techlib import *

# techplot_rms(fitsfile='../data/abcd_006_lpf_new_rms_100_1e-3_cubic.fits',
#              outfile='../plot/abcd_006_rms.pdf',
#              plottext=['6 cm', 'RMS map'],
#              start_coor2=["17:47:20.5", "-28:26:29.780"],
#              end_coor2=["17:47:20.5", "-28:20:29.780"],
#              maplim=[0, 0.5],
#              prolim=[0, 0.5],
#              fromimage=False)
# techplot_rms(fitsfile='../data/rg_096_rms_100_1e-3_cubic.fits',
#              outfile='../plot/096_rms.pdf',
#              plottext=['3 mm', 'RMS map'],
#              start_coor2=["17:47:20.5", "-28:26:29.780"],
#              end_coor2=["17:47:20.5", "-28:20:29.780"],
#              maplim=[0, 0.5],
#              prolim=[0, 0.5],
#              fromimage=False)
techplot_rms(fitsfile='../data/abcd_010_lpf_new_rms_100_5e-4_cubic.fits',
             outfile='../plot/010_rms.pdf',
             plottext=['3 cm', 'RMS map'],
             start_coor2=["17:47:20.5", "-28:26:29.780"],
             end_coor2=["17:47:20.5", "-28:20:29.780"],
             maplim=[0, 0.5],
             prolim=[0, 0.5],
             fromimage=False)