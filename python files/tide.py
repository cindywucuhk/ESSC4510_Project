########################################
### This python file is used to process
### and plot the tide data.
### Name: Wu Hei Tung
### SID : 1155109536
### Reference: 
###   - Detecting Correlations among data-Tutorial.ipynb
###   - eda_06.ipynb
########################################
# import needed packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from math import floor

# import useful function from functions.py
from functions import chebyshevfilt, plot_spectrum

#======================================#
# tide data preprocessing
#======================================#
# import data
tide = pd.read_csv('tidal_data.txt', sep=" ")

# set columns name
tide.columns = ['YY', 'MM', 'DD', 'hh', 'mm', 'ss', 'T', 'Height (m)']

# create datetime columns
col_name = ['YY', 'MM', 'DD', 'hh', 'mm']

for name in col_name:
  tide[name] = tide[name].astype(str)
  tide[name].loc[tide[name].astype(int)<10] = '0' + tide[name].loc[tide[name].astype(int)<10]

tide['datetime'] =  tide['YY'] +  tide['MM'] +  tide['DD'] +  tide['hh'] +  tide['mm']
print(tide['datetime'])

tide['datetime'] = pd.to_datetime(tide['datetime'], format = "%Y%m%d%H%M")
print(tide['datetime'])

#======================================#
# Coping with missing data 1
# sort by date
tide = tide.sort_values(by="datetime")
print(tide)

# drop duplicate tide data
tide_dup = tide[tide.duplicated(subset = ['datetime'], keep = False)]
print(tide_dup)
tide.drop_duplicates(subset=['datetime'], keep = 'first', inplace=True)

# create time Dataframe to add missing data
time_range = pd.date_range(start = datetime(2015,3,1), end = datetime(2021,3,31), freq = '15min')
print(time_range)
time_df = pd.DataFrame(time_range, columns = ['datetime'])
print(time_df)
# merging the tide and time_df
tide = time_df.merge(tide, how='left', on =['datetime'] )
print(tide)

# show nan value
print(tide[(tide['Height (m)'].isna())& (tide['datetime'] < datetime(2021,3,1))& (tide['datetime'] > datetime(2020,12,1))])

#======================================#
# Plot raw tide data
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(tide['datetime'], tide['Height (m)'], 'k-')
ax.set_xlabel('time')
ax.set_ylabel('SSH (m)')
fig.savefig("tide_raw.png")

#======================================#
# get data in np
D = tide['Height (m)'].to_numpy()
print(D)
N = D.shape[0]
t = np.linspace(0,N+1,N)
Dt = t[1] - t[0]
tmax = t[-1]
tmin = t[0]
Nf = floor(N/2+1)
print(N, t)

#======================================#
# Coping with missing data 2
# filter out the 9999 data
D[D == 9999] = np.nan
print(D)

# get the missing data series
Missing = np.zeros((N))
Missing[np.isnan(D)] = 1
print(Missing)

#======================================#
# fourier transform of missing data
Msfft = np.zeros((N,1), dtype=np.complex )
Msfft[:,0] = np.fft.fft(Missing.ravel())

#======================================#
# power spectral density of missing data
norm_ms = (Dt**2)/tmax
Mss = np.zeros((Nf,1))
Mss[:,0]=norm_ms*np.power(abs(Msfft[0:Nf,0]),2)

#======================================#
# Demean
D = D - np.nanmean(D)
# sub 0 into nan data
D[np.isnan(D)] = 0
print(D)

#======================================#
# Demean plot
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
ax1.plot(tide['datetime'],D,'k-')
ax1.set_xlabel('time')
ax1.set_ylabel('Demeaned Height (m)')
fig1.savefig("Demeaned_tide.png")

#======================================#
# get needed variables such as frequencies
# angular frequencies and number of data
fmax = 1/(2*Dt)
No2 = floor(N/2)
Df=fmax/No2
f = np.zeros((Nf,1))
f[:,0] = np.linspace(0,fmax,Nf)
Nw=Nf
wmax = 2*np.pi*fmax
Dw = wmax/(N/2)
w = np.zeros((Nw,1))
w[:,0] = np.linspace(0,wmax,Nw)
fpos=np.zeros((Nf,1))
print((Df * np.linspace(0,No2,Nf)).shape)
fpos[:,0] = Df * np.linspace(0,No2,Nf)
wpos=np.zeros((Nf,1))
wpos[:,0] = Dw * np.linspace(0,No2,Nf)

#======================================#
# Hamming window function, W
W = np.zeros((N,1))
Nw = N
W[0:Nw,0]=0.54-0.46*np.cos(2*np.pi*np.linspace(0,Nw-1,Nw)/(Nw-1))
print(W)
# windowed timeseries
Wd = np.zeros((N,1))
Wd[:,0] = np.multiply( W.T, D ).ravel()
print(Wd)

#======================================#
# plot hammed data
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
ax2.plot(tide['datetime'],Wd,'k-')
ax2.set_xlabel('time')
ax2.set_ylabel('Demeaned Height (m)')
fig2.savefig("Hammed_tide.png")

#======================================#
# high pass filter
fcenter = 1/(24*60/15)
flow = 0.75*fcenter
fhigh = 1.25*fcenter

d = np.zeros((N,1))
Df = np.zeros((N,2))
Df[:,0] = t

d[:,0] = Wd[:,0]
[z, u, v] = chebyshevfilt(d[:,0], Dt, flow, fhigh);
Df[:,1] = z.ravel()

#======================================#
# plot Filtered data
fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.plot(tide['datetime'],Df[:,1].ravel(),'k-')
ax3.set_xlabel('time')
ax3.set_ylabel('Demeaned Height (m)')
fig3.savefig("Filtered_tide.png")

#======================================#
# fourier transform
Wdfft = np.zeros((N,1), dtype=np.complex )
Wdfft[:,0] = np.fft.fft(Df[:,1].ravel())

#======================================#
# power spectral density
norm = (Dt**2)/tmax
Wds = np.zeros((Nf,1))
Wds[:,0]=norm*np.power(abs(Wdfft[0:Nf,0]),2)

#======================================#
# plot power spectral density
# Tide
fig4 = plot_spectrum(f*60/15, Wds, Nf)
fig4.savefig('tidal_power_spectrum_filtered.png')

# plot spectrum
fig5 = plot_spectrum(f*60/15, Msfft, Nf)
fig5.savefig('missing_tidal_power_spectrum_overall.png')

#======================================#
### END ###
#======================================#
