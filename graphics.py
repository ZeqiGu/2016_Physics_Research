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
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.coordinates import get_sun
'''def sunpos(ctime1,ctime2,interval):
    start = Time('2012-7-13 00:00:00') - utcoffset
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times = midnight + delta_midnight
    sunequa = get_sun(times).transform_to('icrs')
    return sunequa

image1 = pyfits.open('cmb_day_pa2_150_2way_1_sun_map0050.fits');
pxls1 = image1[0].data;
print(pxls1[1][1])

#print(image1.info())
(row,col)=pxls1.shape;
        
f=open('PA2_Center_Sidelobes_06012017.txt');
ll=0;
#radial=0;
azimuth=0;
radius=90;
for l in f:
    ll=ll+1;
    print("line %s" %ll)
    a=l.split();
    while azimuth<360:
        
        for radial in np.arange(0,91,0.5):
            r=(radius*radial)/90;
            pix[r*math.cos(azimuth)][r*math.sin(azimuth)]=a[radial*2];
    azimuth+=2;'''
    
    
f=open('PA2_Center_Sidelobes_06012017.txt');
#pix=np.zeros((200,200));
pix=np.zeros((180,181));
az=-2;
for l in f:
    a=l.split();
    print(l)
    if len(a) == 0: break;
    az+=2;
    print(az)
    for radial in range(1,182):
        #pix[az][radial-1]=float(a[radial]);
        pix[int((radial-1)*0.5*math.sin(az))+90][int((radial-1)*0.5*math.cos(az))+90]=float(a[radial]);
        #pix[int((az)*0.5*math.sin(radial))+90][int((az)*0.5*math.cos(radial))+90]=float(a[radial]);
        #print(radial)
        




fig, ax = plt.subplots()
plt.imshow(pix,origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
#plt.imshow(y4,origin='lower',cmap=mpl.cm.get_cmap('jet'));
plt.colorbar(ax=ax)
plt.show()    
    
    
    
    