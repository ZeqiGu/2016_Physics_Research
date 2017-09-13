#!/usr/bin/env python3
from __future__ import division
'''Draw a plot showing the observed area (in RA, Dec)
of the Atacama Cosmology Telescope in 2016. The color
is based on the length of time the telecsope spent
in that area. Data file needed.'''
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.image as mpimg
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import healpy as hp
import pyfits
from astropy import units as u
from astropy.coordinates import SkyCoord
from PIL import Image as img
import math

#f1=open('PA2_Center_Sidelobes_06012017.txt');
rad=0.01745329252;

'''def sidelobeAtPix(pix,y0,x0):
    azimuth=0;
    for l1 in f1:
    #line+=1
        a=l1.split();
        print("azimuth %f" % azimuth)
        if len(a) == 0: break;
        cosmin=min(math.cos(azimuth),math.cos(azimuth+2*rad));
        cosmax=max(math.cos(azimuth),math.cos(azimuth+2*rad));
        #print(az)
        for radial in range(1,182):
            for r in np.arange((radial-1)*0.5,(radial)*0.5,0.3):
                for x in range(int(r*cosmin),int(r*cosmax)):
                    if (r*r-x*x)<0:
                        y=0;
                        #print(line)
                    if azimuth>=180*rad:y=-int(math.sqrt(r*r-x*x));
                    else: y=int(math.sqrt(r*r-x*x));
                    pix[int(y+y0)][int(x+x0)]+=float(a[radial])*100;
        azimuth+=2*rad;'''
            
    
    
 
    
    

pixel=1;
pic= np.zeros((int(180/pixel),int(360/pixel)));
f=open('2017_final.txt');
ll=0;
for l in f:
    ll=ll+1;
    print("line %s" %ll)
    a=l.split();
    if a[0]=="askans_schedule":
        print (l.strip())
        continue;
    length=int((a[3]))-int(float(a[2]));
    alt=int(float(a[4]));
    throw=int(float(a[6]));
    az=int(float(a[5]));
    #print(az)
    ctime1=int(float(a[2]));
    ctime2=int(float(a[3]));
    ctime=ctime1;
    
    for b in range(int(az-throw), int(az+throw)):
        #sidelobeAtPix(pic,alt,b);
        pic[alt][b]+=(100000/throw);
        azimuth=0;
        f1=open('PA2_Center_Sidelobes_06012017.txt');  
        for l1 in f1:
        #line+=1
            a1=l1.split();
            #print("azimuth %f" % azimuth)
            if len(a1) == 0: break;
            cosmin=min(math.cos(azimuth),math.cos(azimuth+2*rad));
            cosmax=max(math.cos(azimuth),math.cos(azimuth+2*rad));
            #print(az)
            for radial in range(1,182):
                for r in np.arange((radial-1)*0.5,(radial)*0.5,0.3):
                    for x in range(int(r*cosmin),int(r*cosmax)):
                        if (r*r-x*x)<0:
                            y=0;
                            #print(line)
                        if azimuth>=180*rad:y=-int(math.sqrt(r*r-x*x));
                        else: y=int(math.sqrt(r*r-x*x));
                        pic[int(y+alt)][int(x+b)]+=float(a1[radial]);
            azimuth+=2*rad;
                  
        
        
            
fig, ax = plt.subplots();
#ax.set_ylim([0,90])
plt.imshow(pic,origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.ylabel('Alt (pixel)')
plt.xlabel('Az (pixel)')
plt.colorbar(ax=ax)
plt.show()
savefig('sim+depth.jpg');