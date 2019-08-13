# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pylab import *
import scipy
import aplpy
import pyparsing
import pyregion
#import pyfits
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib import gridspec
from astropy import wcs
import astropy

#h = fits.getheader('F:/fits_bum/n4522_hi_new.fits')
#f = fits.open('F:/fits_bum/n4522_hi_new.fits')

#w = wcs.WCS(f[0].header)

#print w.wcs.naxis



#print h

#print f

#ha = fits.open('F:/fits_bum/n4522_ha.fits')
#co = fits.open('F:/fits_bum/12n4522.mom0.pbcor.fits')
#ir = fits.open('F:/fits_bum/n4522_24.fits')
#hi = fits.open('F:/fits_bum/n4522_hi.mom0.fits')

#hah = ha[0].header
#ha = ha[0].data
#coh = co[0].header
#co = co[0].data
#irh = ir[0].header
#ir = ir[0].data
#hih = hi[0].header
#hi = hi[0].data


#base = '/media/woorak/data/fits_bum/'
base = 'F:/fits_bum/'
#fig=plt.figure(figsize=(7,10))
fig=plt.figure(figsize=(4.5,4.5))
#fig=plt.subplot(2,1,1)
#fig, (gc1, gc2) = plt.subplots(1,2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#gc3=aplpy.FITSFigure(base+'n4522_hi_region_neww.fits', figure=fig,subplot=[0.1,0.15,0.4,0.7])
gc1=aplpy.FITSFigure(base+'n4522_co_convol.fits', figure=fig,subplot=[0.02,0.02,0.96,0.96])
#gc1.set_frame_color('black')
gc1.axis_labels.hide()
gc1.tick_labels.hide()
#gc3.axis_labels.hide()
#gc3.tick_labels.hide()
#gc3.set_tick_labels_font(size='7')
#gc3.set_axis_labels_font(size='7')
#gc1.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm')
#gc3.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm')
#gc3.set_tick_color('black')
#gc3.set_system_latex(True)
#gc1.show_grayscale(vmin=4000, vmax=10000, invert='true')#n4522_ha_newhead_convol.fits
gc1.show_contour(base+'n4522_opt_r.fits', colors="grey", levels=[15,300], filled='Ture')
gc1.show_contour(base+'n4522_hi_region.fits', colors="k", levels=[0.05,0.2,0.3,0.4,0.5])#, filled='Ture')
gc1.show_contour(base+'n4522_co_mom0.pbcor.fits', colors="red", levels=[0.065,9*0.065,25*0.065,36*0.065,49*0.065],linewidths=1.5)
gc1.show_contour(base+'NGC4522_I.fits', colors="cyan", levels=[0.000028*3,0.000028*5,0.000028*7,0.000028*10])
gc1.add_beam()
gc1.beam.set_color('cyan')
gc1.beam.set_facecolor('none')
gc1.show_contour(base+'n4522_ha_newhead2.fits', colors="blue", levels=[25*3,25*6,25*9,25*15,3000], filled='Ture')
#gc1.set_title('NGC 4522')
#gc1.show_regions(base+'new_region_cross_final2_ds9.reg')

#gc1.xlim(188.433,188.5)
#gc3.show_contour(base+'n4522_opt_r.fits', colors="black", levels=[10,200], filled='Ture')
#gc3.show_contour(base+'n4522_hi_region.fits', colors="blue", levels=[0.05,0.9], filled='Ture')
#gc3.show_contour(base+'n4522_ha_newhead2.fits', colors="red", levels=[80,10000], filled='Ture')
#gc3.show_contour(base+'n4522_co_mom0.pbcor.fits', colors="y", levels=[0.17,0.5,1,3],linewidths=1.5)
#gc3.show_regions(base+'new_region_cross_final2_ds9.reg')
#gc3.add_beam()
#gc3.beam.set_color('blue')
#gc3.show_contour('12n4330_aver_10_convolve.fits', colors="white", levels=[5,15,25,35,45],linewidths=0.7)
#levels=[1.5,3,6,9,12,15,18,21,23]
#levels=[1.5,6,12,18,23]

#gc3.show_beam(major=.0073222, minor=.0066611,angle=-56., fill=False, color="red", corner='bottom left')
#gc3.show_beam(major=.0073222, minor=.0066611,angle=-56., fill=False, color="red", corner='bottom left', hatch='/')
#gc3.show_beam(major=.0017639, minor=.0012417,angle=-28.4, fill=True, color="blue", corner='bottom left', pad=1)
gc1.add_scalebar(0.016666/6.)# , label='10 arcsec', corner='top right') #0.016666
gc1.scalebar.set_corner('top right')
gc1.scalebar.set_label('10 arcsec')
gc1.scalebar.set_linewidth(3)

#gc3.show_scalebar(length=0.016666 , label='1 arcmin', corner='top left') #0.016666

#gc1.show_arrows(188.443, 9.163, -0.012, 0.004250, width=1, head_width=3, head_length=2, color="black")
#gc3.show_arrows(188.443, 9.163, -0.012, 0.004250, width=1, head_width=3, head_length=2, color="black")
#gc1.show_arrows(188.99008, 9.35006, -0.00066, -0.007042, width=3, head_width=7, head_length=5, color="black")
#gc1.recenter()
plt.savefig('combine_new2.png',dpi=500)
plt.show()
