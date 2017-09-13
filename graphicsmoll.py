import matplotlib
matplotlib.use('TkAgg')
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import pylab
import pyfits
import sys
from PIL import Image
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.image as mpimg
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import math


image = hp.read_map('COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits');
#wmap_map_I = hp.read_map('../healpy/test/data/wmap_band_imap_r9_7yr_W_v4.fits')
#image1=hp.mrdfits('COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
#cmb = image[1].data;
#wmap_map_I = image[1].nest2ring(2048, 50331648);
fig, ax = plt.subplots()
cmb=hp.cartview(image, coord=['G','C'], title='Histogram equalized Ecliptic', unit='mK', norm='hist', min=-1,max=1, xsize=2000,return_projected_map=True)

ax.imshow(cmb,origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
#axes.set_xlim([xmin,xmax])
print(ax.get_ylim())

plt.show()



#masked
'''wmap_map_I = hp.read_map('COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
mask = hp.read_map('COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits').astype(np.bool)
wmap_map_I_masked = hp.ma(wmap_map_I)
wmap_map_I_masked.mask = np.logical_not(mask)
cmb=hp.mollview(wmap_map_I_masked.filled(),return_projected_map=True)
print(cmb.shape);
show()'''