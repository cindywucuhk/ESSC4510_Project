########################################
### This python file contain some modified
### useful functions. 
### Name: Wu Hei Tung
### SID: 1155109536
### Reference: 
###   - Detecting Correlations among data-Tutorial.ipynb
###   - eda_06.ipynb
########################################
import numpy as np
import matplotlib.pyplot as plt
from math import floor
import scipy.signal as sg

# 1) plot spectrum
def plot_spectrum(f, s, Nf):
    """ This function plot amplitude/power spectral desity plot.
    ## Tmax is the x-axis scale.
    Argument:
      f [types] - frequency series. [np array]
      s [types] - power spectral density/ amplitude spectral density. [np array]
      Nf[types] - number of data. [int]
    Returns:
    fig2 - matplotlib figure plotted. 
    """
    
    fig2 = plt.figure(figsize=(10, 4))
    #plt.cla()
    #plt.rcParams['xtick.bottom'] = True
    #plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = False
    #plt.rcParams['xtick.labeltop'] = False

    ax1 = fig2.add_subplot(1, 3, 1)
    plt.axis([0, 0.5, 0, 1.1 * np.max(s[1:Nf, 0])])
    plt.plot(f[1:Nf, 0], s[1:Nf, 0], 'k-')
    plt.xlabel('frequency f')
    plt.ylabel('amplitude spectrum')
    plt.title('linear in frequency')

    ax1 = fig2.add_subplot(1, 3, 2)
    Tmax = 10000
    plt.axis([0, Tmax, 0, 1.1*np.max(abs(s[1:Nf, 0].ravel()))])
    print(np.max(s[1:Nf, 0].ravel()))
    plt.plot(np.reciprocal(f[1:Nf, 0]), s[1:Nf, 0], 'k-')
    plt.xlabel('period T (hour)')
    plt.ylabel('amplitude spectrum')
    plt.title('linear in period')

    ax1 = fig2.add_subplot(1, 3, 3)
    plt.semilogx(f[1:Nf, 0], s[1:Nf, 0], 'k-', basex=10)
    plt.axis([0.001, 100, 0, 1.1 * np.max(s[1:Nf, 0])])
    plt.xlabel('frequency f')
    plt.ylabel('amplitude spectrum')
    plt.title('logarithmic in frequency')

    fig2.tight_layout()
    return fig2

# 2) IIR bandpass filter
def chebyshevfilt(d, Dt, flow, fhigh):
    """chebyshevfilt
    # chebyshev IIR bandpass filter
    # d - input array of data
    # Dt - sampling interval
    # flow - low pass frequency, Hz
    # fhigh - high pass frequency, Hz
    # dout - output array of data
    # u - the numerator filter
    # v - the denominator filter
    # these filters can be used again using dout=filter(u,v,din)
    """

    # make sure input timeseries is a column vector
    s = np.shape(d);
    N = s[0];
    if(N==1):
        dd = np.zeros((N,1));
        dd[:,0] = d;
    else:
        dd=d;
        
    # sampling rate
    rate=1/Dt;

    # ripple parameter, set to ten percent
    ripple=0.1;  

    # normalise frequency
    fl=2.0*flow/rate;
    fh=2.0*fhigh/rate;

    # center frequency 
    cf = 4 * np.tan( (fl*np.pi/2) ) * np.tan( (fh*np.pi/2) );

    # bandwidth
    bw = 2 * ( np.tan( (fh*np.pi/2) ) - np.tan( (fl*np.pi/2) ) );

    # ripple parameter factor
    rpf = np.sqrt((np.sqrt((1.0+1.0/(ripple*ripple))) + 1.0/ripple));
    a = 0.5*(rpf-1.0/rpf);
    b = 0.5*(rpf+1.0/rpf);

    u=np.zeros((5,1));
    v=np.zeros((5,1));
    theta = 3*np.pi/4;
    sr = a * np.cos(theta);
    si = b * np.sin(theta);
    es = np.sqrt(sr*sr+si*si);
    tmp= 16.0 - 16.0*bw*sr + 8.0*cf + 4.0*es*es*bw*bw - 4.0*bw*cf*sr + cf*cf;
    v[0,0] = 1.0;
    v[1,0] = 4.0*(-16.0 + 8.0*bw*sr - 2.0*bw*cf*sr + cf*cf)/tmp;
    v[2,0] = (96.0 - 16.0*cf - 8.0*es*es*bw*bw + 6.0*cf*cf)/tmp;
    v[3,0] = (-64.0 - 32.0*bw*sr + 8.0*bw*cf*sr + 4.0*cf*cf)/tmp;
    v[4,0] = (16.0 + 16.0*bw*sr + 8.0*cf + 4.0*es*es*bw*bw + 4.0*bw*cf*sr + cf*cf)/tmp;
    tmp = 4.0*es*es*bw*bw/tmp;
    u[0,0] = tmp;
    u[1,0] = 0.0;
    u[2,0] = -2.0*tmp;
    u[3,0] = 0.0;
    u[4,0] = tmp;

    dout = sg.lfilter(u.ravel(),v.ravel(),dd.ravel());
    return (dout,u,v);

#======================================#
### End ###
#======================================#
