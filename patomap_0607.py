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

f=open('PA2_Center_Sidelobes_06012017.txt');
pix=np.zeros((200,200));
rad=0.01745329252;


#line=-1;
az=0;
for l in f:
    #line+=1
    a=l.split();
    #print(l)
    if len(a) == 0: break;
    cosmin=min(math.cos(az),math.cos(az+2*rad));
    cosmax=max(math.cos(az),math.cos(az+2*rad));
    #print(az)
    for radial in range(1,182):
        for r in np.arange((radial-1)*0.5,(radial)*0.5,0.3):
            for x in range(int(r*cosmin),int(r*cosmax)):
                if (r*r-x*x)<0:
                    y=0;
                    #print(line)
                if az>=180*rad:y=-int(math.sqrt(r*r-x*x));
                else: y=int(math.sqrt(r*r-x*x));
                pix[int(y+90)][int(x+90)]=float(a[radial]);
    az+=2*rad;
        



fig, ax = plt.subplots()
plt.imshow(pix,origin='lower',interpolation='gaussian', cmap=mpl.cm.get_cmap('jet'))
plt.colorbar(ax=ax)
plt.show()    
savefig('sunmap.png')    
    
    
    