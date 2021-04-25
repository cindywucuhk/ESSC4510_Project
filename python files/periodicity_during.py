#!/usr/bin/env python
# coding: utf-8

# In[5]:


import os
from datetime import date
from math import exp, pi, sin, sqrt, floor, ceil
import numpy as np
import scipy.linalg as la
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap

def eda_draw(*argv):
    bw = np.zeros((256,4));
    v = 0.9*(256 - np.linspace( 0, 255, 256 ))/255;
    bw[:,0] = v;
    bw[:,1] = v;
    bw[:,2] = v;
    bw[:,3] = np.ones(256);
    bwcmap = ListedColormap(bw);
    # size of plot
    W = 24;
    H = 6;
    fig1 = plt.figure(1);
    # figsize width and height in inches
    fig1.set_size_inches(W,H);
    ax1 = plt.subplot(1,1,1);
    plt.axis([0, W, -H/2, H/2]);
    plt.axis('off');
    LM = W/6;    # matrix width and heoght
    LV = W/40;   # vector width
    FS = 0.12;    # character width
    TO = 0.4;    # title vertical offset
    SP = 0.2;    # space between objects
    LS = 0.2;    # leading space
    p = LS; # starting x-position
    istitle=0; # flags presence of a title
    for a in argv:
        if isinstance(a,np.ndarray):
            sh = np.shape(a);
            if len(sh) == 1:  # conversion to nx1 array
                n = sh[0];
                m = 1;
                ap = a;
                a = np.zeros((n,1));
                a[:,0] = ap;
            else:
                n = sh[0];
                m = sh[1];
            if m==1:
                pold=p;
                left=p;
                right=p+LV;
                bottom=-LM/2;
                top=LM/2;
                plt.imshow( a, cmap=bwcmap, vmin=np.min(a), vmax=np.max(a), extent=(left,right,bottom,top) );
                p = p+LV;
                pm = (p+pold)/2;
                if istitle:
                    plt.text(pm,-(LM/2)-TO,titlestr,horizontalalignment='center');
                    istitle=0;
                p = p+SP;
            else:
                pold=p;
                left=p;
                right=p+LM;
                bottom=-LM/2;
                top=LM/2;
                plt.imshow( a, cmap=bwcmap, vmin=np.min(a), vmax=np.max(a), extent=(left,right,bottom,top) );
                p = p+LM;
                pm = (p+pold)/2;
                if istitle:
                    plt.text(pm,-(LM/2)-TO,titlestr,horizontalalignment='center');
                    istitle=0;
                p = p+SP;
        elif isinstance(a,str):
            ns = len(a);
            istitle=0;
            if( ns>=6 ):
                if 'title ' in a[0:6]:
                    istitle=1;
                    titlestr=a[6:];
            if( istitle != 1):
                plt.text(p,0,a);
                p = p + ns*FS + SP;
    #plt.show();
    return 1;

def plot_spectrum(f, s, Nf):
    fig2 = plt.figure(2,figsize=(10,4));
    plt.rcParams['xtick.bottom'] = True;
    plt.rcParams['xtick.labelbottom'] = True;
    plt.rcParams['xtick.top'] = False;
    plt.rcParams['xtick.labeltop'] = False;

    ax1 = plt.subplot(1,3,1);
    plt.axis([0, f[-1,0], 0, 1.1*np.max(s[1:Nf,0]) ]);
    plt.plot(f[1:Nf,0],s[1:Nf,0],'k-');
    plt.xlabel('frequency f (cycles per hour)');
    plt.ylabel('amplitude spectrum');
    plt.title('linear in frequency');

    ax1 = plt.subplot(1,3,2);
    Tmax = 1000;
    plt.axis([0, Tmax, 0, 1.1*np.max(s[1:Nf,0]) ]);
    plt.plot(np.reciprocal(f[1:Nf,0]),s[1:Nf,0],'k-');
    plt.xlabel('period T (hour)');
    plt.ylabel('amplitude spectrum');
    plt.title('linear in period');

    ax1 = plt.subplot(1,3,3);
    plt.semilogx(f[1:Nf,0], s[1:Nf,0], 'k-', basex=10);
    plt.axis([0.001, 1, 0, 1.1*np.max(s[1:Nf,0])])
    plt.xlabel('frequency f (cycles per hour)');
    plt.ylabel('amplitude spectrum');
    plt.title('logarithmic in frequency');

    fig2.tight_layout()
    #plt.show();
    return fig2
    
def time_series(N, dt, tmin):
    """ Create auxillary variables t, f and w
        :param N: even integer 
            No. of samples.
        :param dt: float
            Sampling spacing of the data
        :param tmin: float
            Start time of the series
        
        :return: t (np.array,(N,1))
            time series
        :return: f (np.array,(N/2,1))
            frequency series
        :return: w (np.array, (N/2,1))
            angular frequency sereis
        :return: Nf (int)
            floor(N/2+1)
    """
    # making the auxillary variable t       
    tmax = tmin + dt*(N-1)
    t = np.zeros((N,1))
    t[:,0] = np.linspace(tmin, tmax, N)

    # Nyquist frequencies and frequency
    Nf = floor(N/2+1)
    fmax = 1/(2*dt)
    df = fmax/(N/2)
    f = np.zeros((Nf,1))
    f[:,0] = np.linspace(0,fmax,Nf)
    w = f * 2*np.pi
    
    return t, f, w, Nf

def Fouier_kernal(N,M,t,w):
    G = np.zeros((N,M))

    # What is the value of zero frequency?
    G[:,0] = np.ones((N,1)).ravel()
    
    Mo2 = floor(M/2)
    for i in range(1,Mo2):
        j = 2*i-1
        k = j+1
        G[:,j]= np.cos(w[i,0]*t).ravel()
        G[:,k]= np.sin(w[i,0]*t).ravel()

    # nyquist column
    G[:,M-1] = np.cos(w[-1,0]*t).ravel()
    
    return G


# In[6]:


import pandas as pd

df = pd.read_csv("during.csv")
print(df)

df['datetime'] = pd.to_datetime( df['datetime'], format = '%Y%m%d%H' )
print(df['datetime'])


# ### No. of events/day against Time (year)

# In[11]:


D = df

[Nraw, K]=D.shape
traw = np.zeros((Nraw,1))
traw[:,0] = D.iloc[:,0]
Dt = traw[1,0]-traw[0,0]
draw = np.zeros((Nraw,1))
draw[:,0] = D.iloc[:,1]

# round off to even number of points
No2 = floor(Nraw/2)
N=2*No2
d = np.zeros((N,1))
d[:,0]=draw[0:N,0]
t = np.zeros((N,1))
t[:,0] = np.linspace(0,N+1,N)


# obtain the parameters from the dataset
Dt = t[1,0] - t[0,0]
tmin = t[0,0]
t, f, w, Nf = time_series(N, Dt, tmin)
Nw = Nf

# plot data
fig1 = plt.figure(1,figsize=(10,4));
plt.rcParams['xtick.bottom'] = True; plt.rcParams['xtick.labelbottom'] = True; plt.rcParams['xtick.top'] = False;plt.rcParams['xtick.labeltop'] = False;
ax1 = plt.subplot(1,1,1);
plt.axis([df['datetime'][0], df['datetime'][Nraw-1], 0, 1.1*np.max(d) ]);
# plot attributes
plt.plot(df['datetime'][0:Nraw-1],d,'k-'); plt.xlabel('time'); plt.ylabel('earthquake rate (no. of events/hr)'); plt.title('Number of earthquakes per hour');
#plt.show();
fig1.savefig("during.png")

# Obatain some parameters and data kernal
Mo2 = No2
M = N
G = Fouier_kernal(N,M,t,w)

# solve least squares problem using
# analytic formula for inv(GT*G)
import timeit
start = timeit.default_timer()  # Compute runtime for the code section

###analytic formula for inv(GT*G)

gtgi = 2*np.ones((M,1))/N
gtgi[0,0]=1/N
gtgi[M-1,0]=1/N
mest = np.multiply(gtgi, np.matmul(G.T, d))


# In[9]:


#  compute fourier coefficients
mest = np.zeros( (Nf,1), dtype=np.cdouble);
mest[:,0] = np.fft.rfft(d,axis=0).ravel()

# frequency
f = np.zeros( (Nf,1))
f[:,0] = np.fft.rfftfreq(N,Dt)

# ampltude spectrum
s = np.zeros( (Nf,1), dtype=np.double);
s[:,0] = np.abs(mest).ravel()

# compute power spectral density. note the normalization of (2/T).
s2 = np.zeros( (Nf,1), dtype=np.double);
s2[:,0] = (np.real(np.multiply(mest.conj(), mest))).ravel()

### plot spectrum
# Require, Amplitude Spectra - s and frequency -f
#fig2 = plot_spectrum(f,s,Nf)
fig3 = plot_spectrum(f,s2,Nf)
#fig2.savefig("during_amplitude.png")
fig3.savefig("during_power.png")

## export s2 as during_power.txt
np.savetxt("during_power.txt", s2)