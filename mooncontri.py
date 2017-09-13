#!/usr/bin/env python3
import matplotlib
matplotlib.use('TkAgg')
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import pylab
import pyfits
import sys
from PIL import Image as img
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
from astropy.coordinates import get_moon

rad=0.01745329252;
latitude = -22.9586111111111;
longitude = -67.7875;

# Draw the sun sidelobe map
f=open('PA2_Center_Sidelobes_06012017.txt');
pix1=np.zeros((181,181));
pix0=np.zeros((181,181));
rad=0.01745329252;
az=0;
az0=0;
for l in f:
    a=l.split();
    if len(a) == 0: break;
    cosmin=min(math.cos(az),math.cos(az+2*rad));
    cosmax=max(math.cos(az),math.cos(az+2*rad));
    # Draw the map in the converted polar coordinate
    for radial in range(1,182):
        for r in np.arange((radial-1)*0.5,(radial)*0.5,0.3):
            for x in range(int(r*cosmin),int(r*cosmax)):
                if (r*r-x*x)<0:
                    y=0;
                y=int(math.sqrt(r*r-x*x));
                if az >= 180*rad:
                    y*=-1;
                pix1[int(y+90)][int(x+90)]=float(a[radial]);
    az+=2*rad;
    # Draw the map line by line directly 
    for radial in range(1,182):
        if az0<180:
            pix0[az0][radial-1]=float(a[radial]);
    az0+=1;
'''for r in range(181):
    for c in range(181):
        #pix0[r][c]=-((r-90)**2+(c-90)**2);
        pix0[r][c]=abs(90-r);
        
for l in range(181):
    print("az %f" % az)
    cosmin=min(math.cos(az),math.cos(az+2*rad));
    cosmax=max(math.cos(az),math.cos(az+2*rad));
    # Draw the map in the converted polar coordinate

    for radial in range(181):
        print("radial %f" % radial)
        for r in np.arange((radial-1)*0.5,(radial)*0.5,0.3):
            for x in range(int(r*cosmin),int(r*cosmax)):
                if (r*r-x*x)<0:
                    y=0;
                #if az>=180*rad:y=-int(math.sqrt(r*r-x*x));
                y=int(math.sqrt(r*r-x*x));
                if az >= 180*rad:
                    y*=-1;
                    print("ok")
                print("x %f" %x, "y %f" %y)
                pix1[int(y+90)][int(x+90)]=float(pix0[l][radial]);
    az+=2*rad;'''
        



# Major Part: First in horizontal coordinate, for each boresight point, calculate and get the right
# pixel of sun's position in the sidelobe map. Then convert the whole map into equatorial coordinate.
f=open('2017_final.txt');
ll=0;
pix=np.zeros((180,360))
[row,col]=pix.shape;
pix2=np.zeros((row, col));
pix3=np.zeros((row, col));
pix00=np.zeros((row, col));
pix4=np.zeros((row, col));
pix5=np.zeros((row, col));
pix6=np.zeros((row, col));
pix7=np.zeros((row, col));
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
    ctime1=int(float(a[2]));
    ctime2=int(float(a[3]));

    for ctime in np.arange(ctime1, ctime2, 180):
        
        t0=Time(ctime, format='unix');
        t=Time(t0.isot, format='isot',scale='utc')
        #print(t);
        x=get_moon(t);
        sunra=x.ra.value;
        sundec=x.dec.value;
        '''sunra=30;
        sundec=-23;'''
        #print("ra %f" % sunra, "dec %f" % sundec)
        julian= (ctime/86400.0) + 2440587.5;
        gmst = 18.697394558 + 24.06570982441908*(julian- 2451545.0)       
        gmst -= floor(gmst/24.0)*24.0;
        #print("now sidereal %f" %gmst)
        lst = gmst+(longitude/15.0);
        lst -= floor(lst/24.0)*24.0;       
        #print("now local sidereal %f" % lst)
        # Besides latitude, the sun's position only changes with time.
        sunLHA=(lst*15-sunra+360) % 360;
        #if sunLHA<0: sunLHA += 24;
        #print("sunLHA %f" % sunLHA)
        
        sunaltitudesin=math.sin(sundec*rad)*math.sin((latitude)*rad)+math.cos(sundec*rad)*math.cos(latitude*rad)*math.cos(sunLHA*rad)
        sunaltitude=math.degrees(math.asin(sunaltitudesin))
        
        sunazimuthsin=-(math.sin(sunLHA*rad)*math.cos(sundec*rad))/math.cos(sunaltitude*rad)
        sunazimuthcos=(math.sin(sundec*rad)-math.sin(latitude*rad)*math.sin(sunaltitude*rad))/(math.cos(latitude*rad)*math.cos(sunaltitude*rad))
        sunazimuth=math.degrees(math.atan2(sunazimuthsin,sunazimuthcos));
        if sunazimuth<0: sunazimuth+=360;
        #print("alt %f" % sunaltitude, "az %f" % sunazimuth)
        
        for b in range(int(az-throw),int(az+throw)):          
            dec1 = math.asin(sin(alt*rad)*sin(latitude*rad) + cos(alt*rad)*cos(latitude*rad)*cos(b*rad));
            LHA=math.degrees(math.acos((sin(alt*rad)-sin(dec1)*sin(latitude*rad))/(cos(dec1)*cos(latitude*rad))));
            LHA1=math.asin(-sin(b*rad)*cos(alt*rad)/cos(dec1));                    
            if LHA1<0: LHA = 360-LHA;
            #print("telescope alt %f" % alt, "az %f" % b)
            # 1 is telescope's; 2 is sun's. Calculate angle distance between boresight and sun, it's not simly subtract
            # their coordinates because this is on a sphere.
            altdeltacos=sin(sunaltitude*rad)*sin(alt*rad)+cos(sunaltitude*rad)*cos(alt*rad)*cos(abs(sunazimuth-b)*rad)
            altdeltasin=math.sqrt((cos(sunaltitude*rad)*sin(abs(sunazimuth-b)*rad))**2+(cos(alt*rad)*sin(sunaltitude*rad)-
                                  sin(alt*rad)*cos(sunaltitude*rad)*cos(abs(sunazimuth-b)*rad))**2);
            #altdelta=math.degrees(math.atan2(altdeltasin,altdeltacos));
            altdelta=math.degrees(math.acos(altdeltacos))
            
            azdeltasin=sin((90-sunaltitude)*rad)*sin((b-sunazimuth)*rad)/sin(altdelta*rad);
            azdeltacos=(cos((90-sunaltitude)*rad)-cos(altdelta*rad)*cos((90-alt)*rad))/(sin(altdelta*rad)*sin((90-alt)*rad));
            azdelta=math.degrees(math.atan2(azdeltasin,azdeltacos))
            if azdelta<0: azdelta+=360;
            azdelta+=90;
            if azdelta>360: azdelta-=360;
            #print("altdelta %f" % altdelta, "azdelta %f" % azdelta)
            
            sunazcoord=altdelta*math.cos(azdelta*rad);
            sunaltcoord=altdelta*math.sin(azdelta*rad)
            # Read directly from the unconverted sidelobe map. Since each row is 2 deg azimuthal
            # and col is 0.5 deg radial.
            try: sunsidelobe=pix1[int(90+sunaltcoord)][int(90+sunazcoord)]/throw;               
            except IndexError: sunsidelobe=0;            
            try: sunsidelobe0=pix0[int(azdelta/2)][int(altdelta*2)]/throw;
            except IndexError: sunsidelobe0=0;
            if sunaltitude<0: sunsidelobe0=0;
            
            # Plot them out in equatorial coordinate
            pix[alt][b]=sunsidelobe0;
            RA=(lst*15-LHA)%360;
            while RA<0: RA+=360;
            dec=math.degrees(dec1);
            pix3[int(dec+90)][int(RA)]+=1000/throw;
            #print(dec,RA)
            if a[1]=="wideAN" or a[1]=="wideAS" or a[1]=="wideB": pix6[int(dec+90)][int(RA)]=sunsidelobe0;
            if a[1]=="day_12h_n" or a[1]=="day_03h_s": pix7[int(dec+90)][int(RA)]=sunsidelobe0;
            if sunazimuth<=0: pix4[int(dec+90)][int(RA)]=sunsidelobe0;
            if sunazimuth>0: pix5[int(dec+90)][int(RA)]=sunsidelobe0;
            pix00[int(dec+90)][int(RA)]=sunsidelobe0;
            pix2[int(dec+90)][int(RA)]=sunsidelobe;

#pix: sidelobe map in horizontal; pix00: sidelobe map in equtorial, by original pato's map;
#pix0: pato's map originally; pix1: pato's map in polar; pix2: sidelobe map in equatorial, by polar pato's map;
#pix3: depth map; pix4: sidelobe map when sun on west, setting; pix5: when sun on east, rising; pix6: sidelobes of "Wide"
#fields; pix7: sidelobes of "Day" fields;

fig, ax = plt.subplots(3,2)
plt.subplot(321)
plt.imshow(pix3,origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.colorbar(ax=ax[0,0])
plt.subplot(322)
pix00[pix00==0]=0.0000000001;
plt.imshow(np.log(pix00),origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.colorbar(ax=ax[0,1])
plt.subplot(323)
pix4[pix4==0]=0.0000000001;
plt.imshow(np.log(pix4),origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.colorbar(ax=ax[1,0])
plt.subplot(324)
pix5[pix5==0]=0.0000000001;
plt.imshow(np.log(pix5),origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.colorbar(ax=ax[1,1])
plt.subplot(325)
pix6[pix6==0]=0.0000000001;
plt.imshow(np.log(pix6),origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.colorbar(ax=ax[2,0])
plt.subplot(326)
pix7[pix7==0]=0.0000000001;
plt.imshow(np.log(pix7),origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.colorbar(ax=ax[2,1])
'''fig, ax=plt.subplots()
plt.imshow(pix00,origin='lower',cmap=mpl.cm.get_cmap('jet'));
plt.colorbar(ax=ax)'''
plt.show()  
#savefig('suncontri.png')

