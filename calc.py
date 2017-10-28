#!/usr/bin/env python3
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.image as mpimg
from pylab import *
import math
import mpld3
from mpl_toolkits.mplot3d import axes3d

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
r=0.22;
omega=math.pi*120.0;
l=0.9;
mag=0.005;
u0=1.2566370614*0.000001;
A=math.pi*r*r;

    
    
N=  float(raw_input("Enter turns, required: "))
Cin=  float(raw_input("Enter capacitance, if this is what you want to calc, enter -1: "))
density=  float(raw_input("Enter density, required: "))
voltin=  float(raw_input("Enter volt, if this is what you want to calc, enter -1: "))


L=u0*N*N*A*0.85/l;
R=N*2*math.pi*r*density;
narr = np.linspace(N*0.8, N*1.2, 10);
I=(2*mag*np.sqrt(l*l+r*r))/(u0*narr);
print(I)


if (voltin!=-1):
    varr= np.linspace(voltin*0.8, voltin*1.2, 10);
    narr, varr = np.meshgrid(narr, varr);
    x=omega*L-np.sqrt(np.absolute(np.multiply(varr,varr)/(I*I)-R*R));
    C=1/(x*omega);
    print(C)
    ax.plot_wireframe(narr, varr, C, rstride=int(N*0.05), cstride=int(voltin*0.05));
    plt.show()
    
    
if (Cin!=-1):
    carr = np.linspace(Cin*0.8, Cin*1.2, 10);
    narr, carr = np.meshgrid(narr, carr);
    V=I*np.sqrt(R*R+np.pow((omega*L-(1/(omega*carr))),2));
    print(V)
    ax.plot_wireframe(narr, carr, V, rstride=int(N*0.05), cstride=int(Cin*0.05));
    plt.show()
    



