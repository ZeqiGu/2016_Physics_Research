#!/usr/bin/env python3

#Combine the depth map with Planck cmb map. Some rounding problems cause discrete grids.
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
'''import profile
import pstats
import cProfile
import re
cProfile.run('re.compile("graphics1.py")', 'restats')
p = pstats.Stats('restats')
p.strip_dirs().sort_stats(-1).print_stats()

p.sort_stats('time').print_stats(10)'''

fig, ax = plt.subplots()
latitude = -22.9586111111111;
longitude = -67.7875;
rad=0.01745329252;
pixel=1;
print(pixel)
# create a numpy array
pic= np.zeros((int(180/pixel),int(360/pixel)));
'''image = hp.read_map('COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits');
cmb=hp.cartview(image, coord=['G','C'], title='Histogram equalized Ecliptic', unit='mK', norm='hist', min=-1,max=1, xsize=2000,return_projected_map=True)
#old
image1=hp.mrdfits('COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits',1)
def cmbval(alpha,beta):
    c_icrs = SkyCoord(ra=alpha*u.degree, dec=beta*u.degree, frame='icrs');
    pixindex=hp.ang2pix(2048,math.radians(90-c_icrs.galactic.b.value),math.radians(c_icrs.galactic.l.value));
    pix=image1[0][pixindex];           
    return float(pix);'''
##        
f=open('2017.txt');
ll=0;
for l in f:
    ll=ll+1;
    #print("line %s" %ll)
    a=l.split();
    if a[0]=="askans_schedule":
        print (l.strip())
        continue;
    length=int((a[3]))-int(float(a[2]));
    alt=int(float(a[4]));
    throw=int(float(a[6]));
    az=int(float(a[5]));
    
    
    ctime1=int(float(a[2]));
    ctime2=int(float(a[3]));
    ctime=ctime1;
    while ctime<ctime2:
        julian= (ctime/86400.0) + 2440587.5;
        gmst = 18.697394558 + 24.06570982441908*(julian- 2451545.0)
        gmst -= floor(gmst/24.0)*24.0;
        lst = gmst+(longitude/15.0);
        lst -= floor(lst/24.0)*24.0;
        b=int(az-throw);
        while b<int(az+throw):
            dec1 = math.asin(sin(alt*rad)*sin(latitude*rad) + cos(alt*rad)*cos(latitude*rad)*cos(b*rad));
            LHA=math.degrees(math.acos((sin(alt*rad)-sin(dec1)*sin(latitude*rad))/(cos(dec1)*cos(latitude*rad))));
            LHA1=math.asin(-sin(b*rad)*cos(alt*rad)/cos(dec1));
                      
            if LHA1<0:
                LHA = 360-LHA;
                        
            RA=int(lst*15-LHA)%360;
            while RA<0:
                RA+=360;
            dec1=int(math.degrees(dec1));
            print("dec %f" % dec1)
            print("RA %f" % RA)
            '''cmb=cmbval(RA,dec1);
            if cmb is None: cmb=0;
            #print(cmb)'''
            
            if pixel<1:
                for t1 in range(int(1/pixel)):
                    for t2 in range(int(1/pixel)):
                        #if t1/2==0:
                        pic[int((90)/pixel+(dec1+t1*pixel)/pixel)][int((RA+t2*pixel)/pixel)]+=(1000/throw);
                       
            else:
                pic[int((90)/pixel+(dec1/pixel))][int(RA/pixel)]=pic[int((90)/pixel+((dec1)/pixel))][int((RA)/pixel)]+(1000/throw);
            #print ("dec1 %s" % pic[dec1+90][RA+90])
            b=b+0.2;
        
        ctime=ctime+120;
        
#max=np.amax(pic);
#print (max);
    
'''[row,col]=pic.shape;
pic1=np.zeros((row, col));
for r in range(row):
    for c in range (col):
        r1=2*r; c1=2*c;
        
        #if c1>col-1: c1=-1;
        if pic[r][c] !=0:
            pic1[r][c]=pic[r][c]*100*np.mean([cmb[r1][c1],cmb[max(0,r1-1)][c1],cmb[max(0,r1-1)][max(0,c1-1)],cmb[max(0,r1-1)][min(col-1,c1+1)],
                                                        cmb[r1][max(0,c1-1)],cmb[r1][min(col-1,c1+1)],cmb[min(row-1,r1+1)][c1],cmb[min(row-1,r1+1)][max(0,c1-1)],cmb[min(row-1,r1+1)][min(col-1,c1+1)]]);
        #if pic[r][c] !=0: pic1[r][c]=pic[r][c]*cmb[r][c];'''


'''ax.imshow(pic,origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.ylabel('Dec (deg)')
plt.xlabel('RA (deg)')
plt.colorbar(ax=ax);
plt.show()'''
fig, ax = plt.subplots()
plt.imshow(pic,origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.colorbar(ax=ax)
plt.show()    
#savefig('obscmb1.jpg');
