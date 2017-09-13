#!/usr/bin/env python3
import matplotlib
matplotlib.use('TkAgg')
import healpy as hp
import matplotlib.pyplot as plt
import pyfits
import sys
from PIL import Image as img
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.image as mpimg
import numpy as np
from pylab import *
import math
from astropy import units as u
from astropy.coordinates import SkyCoord,EarthLocation, GCRS,AltAz
from astropy.time import Time
from scipy import signal
from scipy import misc
from astropy.io import fits
import time
import numba
from numba import jit
# Useful constants
hdulist = fits.open('LFI_SkyMap_030_1024_R2.01_full.fits')
print(hdulist.info())
rad=0.01745329252;
latitude = -22.9586111111111;
longitude = -67.7875;
epsilon=23.43333
raknot=282.25;
decknot=62.6;
longknot=33;
act = EarthLocation(lat=-22.9586111111111*u.deg, lon=-67.7875*u.deg, height=5190*u.m)

@jit(nopython=True)
def func0(alt,rad,latitude,b,lst):
            dec1 = math.asin(sin(alt*rad)*sin(latitude*rad) + cos(alt*rad)*cos(latitude*rad)*cos(b*rad));
            LHA=math.degrees(math.acos((sin(alt*rad)-sin(dec1)*sin(latitude*rad))/(cos(dec1)*cos(latitude*rad))));
            LHA1=math.asin(-sin(b*rad)*cos(alt*rad)/cos(dec1));                    
            if LHA1<0: LHA = 360-LHA;
            RA=(lst*15-LHA) % 360;
            while RA<0: RA+=360;
            dec1=math.degrees(dec1);
            return(RA,dec1)


@jit(nopython=True)
def func1(cmblat,rad,cmblong,decknot,longknot,raknot,lst):
    # Another approach: convert from galactic to equatorial
                        cmbdecsin=cos(cmblat*rad)*sin((cmblong-longknot)*rad)*sin(decknot*rad)+sin(cmblat*rad)*cos(decknot*rad);
                        cmbdec=math.degrees(math.asin(cmbdecsin));                       
                        cmbrasin=(cos(cmblat*rad)*sin((cmblong-longknot)*rad)*cos(decknot*rad)-sin(cmblat*rad)*sin(decknot*rad))/cos(cmbdec)
                        cmbrasin=-(sin(cmblat*rad)-sin(cmbdec*rad)*cos(decknot*rad))/(cos(cmbdec*rad)*sin(decknot*rad));
                        cmbracos=cos(cmblat*rad)*cos((cmblong-longknot)*rad)/cos(cmbdec);
                        cmbra=math.degrees(math.atan2(cmbrasin,cmbracos))+raknot;
                        if cmbra<0: cmbra+=360;
                        if cmbra>360: cmbra-=360;
                        #print("cmbra %f" % cmbra, "cmbdec %f" % cmbdec)
                        # calculate each cmb's position from equtorial to horizontal
                        cmbLHA=(lst*15-cmbra+360) % 360;                 
                        cmbaltitudesin=math.sin(cmbdec*rad)*math.sin((latitude)*rad)+math.cos(cmbdec*rad)*math.cos(latitude*rad)*math.cos(cmbLHA*rad)
                        cmbaltitude=math.degrees(math.asin(cmbaltitudesin))                    
                        cmbazimuthsin=-(math.sin(cmbLHA*rad)*math.cos(cmbdec*rad))/math.cos(cmbaltitude*rad)
                        cmbazimuthcos=(math.sin(cmbdec*rad)-math.sin(latitude*rad)*math.sin(cmbaltitude*rad))/(math.cos(latitude*rad)*math.cos(cmbaltitude*rad))
                        cmbazimuth=math.degrees(math.atan2(cmbazimuthsin,cmbazimuthcos));
                        if cmbazimuth<0: cmbazimuth+=360;
                        #print("cmbaltitude %f" % cmbaltitude, "az %f" % cmbazimuth)
                        #print("galactic to horizontal2--- %s seconds ---" % (time.time() - start_time))
                        return (cmbazimuth,cmbaltitude)

@jit(nopython=True)
def func2(cmbaltitude,b,alt,rad,cmbazimuth):
    # calculate altdelta and azdelta between (cmbaltitude, cmbazimuth) and (alt, b)
                            altdeltacos=sin(cmbaltitude*rad)*sin(alt*rad)+cos(cmbaltitude*rad)*cos(alt*rad)*cos(abs(cmbazimuth-b)*rad)
                            altdeltasin=math.sqrt((cos(cmbaltitude*rad)*sin(abs(cmbazimuth-b)*rad))**2+(cos(alt*rad)*sin(cmbaltitude*rad)-
                                                  sin(alt*rad)*cos(cmbaltitude*rad)*cos(abs(cmbazimuth-b)*rad))**2);
                            #altdelta=math.degrees(math.atan2(altdeltasin,altdeltacos));
                            altdelta=math.degrees(math.acos(altdeltacos))
                            
                            azdeltasin=sin((90-cmbaltitude)*rad)*sin((b-cmbazimuth)*rad)/sin(altdelta*rad);
                            azdeltacos=(cos((90-cmbaltitude)*rad)-cos(altdelta*rad)*cos((90-alt)*rad))/(sin(altdelta*rad)*sin((90-alt)*rad));
                            azdelta=math.degrees(math.atan2(azdeltasin,azdeltacos))                   
                            azdelta=(azdelta+90) % 360;
                            #print("altdelta %f" % altdelta, "azdelta %f" % azdelta)
                            return (azdelta,altdelta)
                            
'''@jit
def func3(convolarr,cmbrow,cmb):
                reversedconvolarr=np.fliplr(np.flipud(convolarr));
                if cmbrow <500: convolution=signal.convolve(cmb, reversedconvolarr, mode='valid');
                else: convolution=signal.convolve2d(cmb, reversedconvolarr, mode='valid');
                return(convolution)'''

# Draw the sun sidelobe map without conversion
f=open('PA2_Center_Sidelobes_06012017.txt');
pix1=np.zeros((181,181));
pix0=np.zeros((181,181));
rad=0.01745329252;
az=0;
for l in f:
    a=l.split();
    if len(a) == 0: break;
    cosmin=min(math.cos(az),math.cos(az+2*rad));
    cosmax=max(math.cos(az),math.cos(az+2*rad));
    # Draw the map line by line directly 
    for radial in range(1,182):
        if az<180:
            pix0[az][radial-1]=float(a[radial]);
    az+=1;

# Get the cmb array
image = hp.read_map('LFI_SkyMap_030_1024_R2.01_full.fits',field=0);
yes=0;
no=0;
# Read map in galactic coordinate. (Don't know whether it is in ra/dec cartesian, upside down or not)
pixel=1;
#cmb=hp.mollview(image, coord=['G','G'], title='Histogram equalized Ecliptic', unit='mK', norm='hist', xsize=int(360/pixel),return_projected_map=True)
cmb=np.ones((180,360))
[cmbrow,cmbcol]=cmb.shape;
cmb[cmb == -inf] = 0
print(np.amax(cmb),np.amin(cmb))
pix3=np.zeros((cmbrow,cmbcol));
# The array for convolution actually takes the depth map's coordinates.
f=open('2017.txt');
ll=0;
# Use a hashmap to increase efficiency. Key: boresight equatorial coordinate(dec1,RA); didn't use horizontal because then
# the data is also related to ctime. Value: convolution result.
convolDict={};
# Key: (boresight's alt, boresight's az, cmb's alt, cmb's az) because the sidelobe data to be extracted is determined
# by them. Value: cmb sidelobes at that pixel. 
cmbDict={};
convolarr= np.zeros((int(180/pixel),int(360/pixel)));
result= np.zeros((int(180/pixel),int(360/pixel))); 
for l in f:
    start_time = time.time()
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
    ctime1=int(float(a[2]));
    ctime2=int(float(a[3]));
    for ctime in np.arange(ctime1, ctime2, 180):
        julian= (ctime/86400.0) + 2440587.5;
        gmst = 18.697394558 + 24.06570982441908*(julian- 2451545.0)       
        gmst -= floor(gmst/24.0)*24.0;
        lst = gmst+(longitude/15.0);
        lst -= floor(lst/24.0)*24.0;       
        t0=Time(ctime, format='unix');
        '''t00=Time.to_datetime(t0)
        if t00.month>7: utcoffset = -4*u.hour;  # Originally for calculating local time to fill the "obstime" 
        else:utcoffset = -5*u.hour'''           # in following functions, but then found that it takes utc time.
        t=Time(t0.isot, format='isot',scale='utc')
        #print(t0.isot)
        start_time = time.time();    
        for b in np.arange(int(az-throw),int(az+throw),1): # for all the pixels in horizontal depth map
            # about the pixel number: two ways to make the dimensions of depth map and cmb map match: first is to
            # directly get a rougher cmb map from read_map function; second is to expand depth map to accomdate the
            # higher resolution of the cmb map.
            # calculate telescope's bore sight in equatorial coordinate
            (RA,dec1)=func0(alt,rad,latitude,b,lst)
            #print((RA,dec1))
            #print(dec1/pixel,RA/pixel)
            if int(RA)==100 and no<2:
                no+=1
                print("no")
                yes=0;
                
            #if convolDict.get((int(alt),int(b))) == None:              
                for cmblatpix in np.arange(cmbrow): #this should be in the read file loop because ctime is involved
                    cmblat=cmblatpix*180/cmbrow-90;
                    if yes==0:                  
                        for cmblongpix in range(cmbcol):
                                    
                                    # calculate each cmb's position from galactic to horizontal                       
                                    cmblong=cmblongpix*360/cmbcol;
                                    '''cmb_galactic = SkyCoord(l=cmblong*u.degree, b=cmblat*u.degree,frame='galactic')
                                    cmb_horizontal=cmb_galactic.transform_to(AltAz(obstime=t, location=act))
                                    cmbaltitude=cmb_horizontal.alt.value;
                                    cmbazimuth=cmb_horizontal.az.value;'''                       
                                    (cmbazimuth,cmbaltitude)=func1(cmblat,rad,cmblong,decknot,longknot,raknot,lst);                        
                                    index=(int(round(alt)),int(round(b)),int(round(cmbaltitude)),int(round(cmbazimuth)))
                                    #if cmbDict.get(index) == None:
                                    #if index not in cmbDict:
                                        # calculate altdelta and azdelta between (cmbaltitude, cmbazimuth) and (alt, b)                           
                                    '''boresight=SkyCoord(ra=RA*u.degree, dec=dec1*u.degree,frame='gcrs');
                                        altdelta=cmb_horizontal.separation(boresight).value'''                            
                                    (azdelta,altdelta)=func2(cmbaltitude,b,alt,rad,cmbazimuth)
                                    #print((cmbaltitude,cmbazimuth,alt,b))
                                        # extract the sidelobe data
                                    try: cmbsidelobe0=pix0[int(azdelta/2)][int(altdelta*2)]/throw;
                                    except IndexError: cmbsidelobe0=0;
                                    if cmbaltitude<0: cmbsidelobe0=0;                          
                                    #cmbDict[index]=cmbsidelobe0;
                                        #print("not use cmbDict")
                                    #else:
                                                #cmbsidelobe0=cmbDict.get(index);
                                                #print("use cmbDict")
                
                                    #put data on the correct position(s). The coordinate system should be the same with cmb map.
                                    #convolarr updates itself every loop
                                    #convolarr[int(cmbdec+90)][int(cmbra)]=cmbsidelobe0;
                                    convolarr[cmblatpix][cmblongpix]=cmbsidelobe0;
                                    yes=1;
                # convolution of cmb and convolarr
                #convolution=func3(convolarr,cmbrow,cmb) This instead slows it down
                                    reversedconvolarr=np.fliplr(np.flipud(convolarr));
                                    if cmbrow <500: convolution=signal.convolve(cmb, reversedconvolarr, mode='valid');
                                    else: convolution=signal.convolve2d(cmb, reversedconvolarr, mode='valid');
                #convolDict[(int(alt),int(b))]=convolution;
                
            #else:
                        
                        #convolution=convolDict.get((int(alt),int(b)));
                        #print("use convolDict")
            # put the result to the right position(s): no matter pixel >1 or <1 (higher or lowe resolution), since the result
            # array has the same size of cmb map, we can always assign the data one-to-one to the (dec1, RA) looping
            # coordinate of the result array. If the pixel<1, the difference of dec and RA will be enlarged themselves so the
            # overlapped indices will reduce.
                                    result[int((90)/pixel+(dec1/pixel))][int(RA/pixel)]+=convolution;
                                    pix3[int((90)/pixel+(dec1/pixel))][int(RA/pixel)]+=1000/throw;
        #print("loop--- %s seconds ---" % (time.time() - start_time))
                                    print(cmblatpix,cmblongpix,azdelta, altdelta,convolution)
            else: continue;
    ''' If use this for pixel<1, the resolution will increase although the map is enlarged.
            if pixel<1:
                for t1 in range(int((90+dec1)/pixel),int((90+dec1)/pixel+1/pixel)):
                    for t2 in range (int(RA/pixel),int(RA/pixel+1/pixel)):
                        result[t1][t2]+=convolution;
            else:
                result[int((90)/pixel+(dec1/pixel))][int(RA/pixel)]+=convolution;'''
            
            
#print(np.nonzero(result))           
fig, ax=plt.subplots(2,1)
plt.subplot(211)
result[result==-inf] = 0
plt.imshow(result,origin='lower',cmap=mpl.cm.get_cmap('jet'));
plt.colorbar(ax=ax[0])
plt.subplot(212)
plt.imshow(pix3,origin='lower',cmap=mpl.cm.get_cmap('jet'));
plt.colorbar(ax=ax[1])
plt.show()  

        
        
        
        
        
        
        
        
        
        
        
''' Self-written code for coordinate conversions; proved to be correct.
#Convert cmb from ecliptic to equatorial
                    cmbdecsin=sin(cmblat*rad)*cos(epsilon*rad)+cos(cmblat*rad)*sin(epsilon*rad)*sin(cmblong*rad);
                    cmbdec=math.degrees(math.asin(cmbdecsin));
                    
                    cmbrasin=-(sin(cmblat*rad)-sin(cmbdec*rad)*cos(epsilon*rad))/(cos(cmbdec*rad)*sin(epsilon*rad));
                    cmbracos=sin(cmblong*rad)*cos(cmblat*rad)/cos(cmbdec*rad);
                    cmbra=math.degrees(math.atan2(cmbrasin,cmbracos));
                    if cmbra<0: cmbra+=360;
                    #print("cmbra %f" % cmbra, "cmbdec %f" % cmbdec)
                    # Another approach: convert from galactic to equatorial
                    cmbdecsin=cos(cmblat*rad)*sin((cmblong-longknot)*rad)*sin(decknot*rad)+sin(cmblat*rad)*cos(decknot*rad);
                    cmbdec=math.degrees(math.asin(cmbdecsin));
                    
                    cmbrasin=(cos(cmblat*rad)*sin((cmblong-longknot)*rad)*cos(decknot*rad)-sin(cmblat*rad)*sin(decknot*rad))/cos(cmbdec)
                    cmbrasin=-(sin(cmblat*rad)-sin(cmbdec*rad)*cos(decknot*rad))/(cos(cmbdec*rad)*sin(decknot*rad));
                    cmbracos=cos(cmblat*rad)*cos((cmblong-longknot)*rad)/cos(cmbdec);
                    cmbra=math.degrees(math.atan2(cmbrasin,cmbracos))+raknot;
                    if cmbra<0: cmbra+=360;
                    if cmbra>360: cmbra-=360;
                    #print("cmbra %f" % cmbra, "cmbdec %f" % cmbdec)
                    # calculate each cmb's position from equtorial to horizontal
                    cmbLHA=(lst*15-cmbra+360) % 360;                 
                    cmbaltitudesin=math.sin(cmbdec*rad)*math.sin((latitude)*rad)+math.cos(cmbdec*rad)*math.cos(latitude*rad)*math.cos(cmbLHA*rad)
                    cmbaltitude=math.degrees(math.asin(cmbaltitudesin))                    
                    cmbazimuthsin=-(math.sin(cmbLHA*rad)*math.cos(cmbdec*rad))/math.cos(cmbaltitude*rad)
                    cmbazimuthcos=(math.sin(cmbdec*rad)-math.sin(latitude*rad)*math.sin(cmbaltitude*rad))/(math.cos(latitude*rad)*math.cos(cmbaltitude*rad))
                    cmbazimuth=math.degrees(math.atan2(cmbazimuthsin,cmbazimuthcos));
                    if cmbazimuth<0: cmbazimuth+=360;
                    print("cmbaltitude %f" % cmbaltitude, "az %f" % cmbazimuth)
                    # calculate altdelta and azdelta between (cmbaltitude, cmbazimuth) and (alt, b)
                    altdeltacos=sin(cmbaltitude*rad)*sin(alt*rad)+cos(cmbaltitude*rad)*cos(alt*rad)*cos(abs(cmbazimuth-b)*rad)
                    altdeltasin=math.sqrt((cos(cmbaltitude*rad)*sin(abs(cmbazimuth-b)*rad))**2+(cos(alt*rad)*sin(cmbaltitude*rad)-
                                          sin(alt*rad)*cos(cmbaltitude*rad)*cos(abs(cmbazimuth-b)*rad))**2);
                    #altdelta=math.degrees(math.atan2(altdeltasin,altdeltacos));
                    altdelta=math.degrees(math.acos(altdeltacos))'''

                    
                




