########################################
### This python file downloaded, seperate
### and create csv file for earthquake data.
### Calculate earthquake rate
### Name: Wu Hei Tung
### SID: 1155109536
########################################
# import needed packages
import pandas as pd
import numpy as np
from urllib.request import urlretrieve 
from urllib import error 
from datetime import datetime

#======================================#
# download earthquake data from website
# set path
PATH = "earthquake_data.txt"

# set url
url = "http://axial.ocean.washington.edu/hypo71.dat"

# retrieved data
try:
  urlretrieve(url, PATH)
except error.HTTPError:
  pass

#======================================#
# save the raw data into csv
# open using np
D = pd.read_fwf('earthquake_data.txt')
D.to_csv('earthquake_data.csv', index = False)

#======================================#
# Earthquake data processing 1
# Data format processing
df = D
print(df, df.columns)

# get date
df['datetime'] = df['yyyymmdd'].astype(str) + df['HHMMSSS.SS']
print(df['datetime'])
df['datetime'] = pd.to_datetime(df['datetime'], format = '%Y%m%d%H%M %S.%f')
target_dates = (df['datetime'] > datetime(2015,4,23)) & (df['datetime'] < datetime(2015,5,22))
t = df[target_dates]['datetime']
print(t)

# seperate magnitude NWR
MW = df['MW NWR'].to_numpy()
MW_new = np.zeros((MW.shape[0], 1))
NWR_new = np.zeros((MW.shape[0], 1))

for i in range(0, MW.shape[0]):
  if MW[i][0:3] != 'NaN':
    MW_new[i,0] = MW[i][0:5]
  else:
    MW_new[i,0] = np.nan
  NWR_new[i,0] = MW[i][-2:]
print(MW_new, NWR_new)

# create new column that store MW and NWR
df['MW'] = MW_new.ravel()
df['NWR'] = NWR_new.ravel()

#======================================#
## create a better csv_file
# change the datetime back to string
df['datetime'] = df['datetime'].dt.strftime('%Y%m%d%H%M%S%f')
print(df)
# drop unnecessary column
df = df.drop(
  ['yyyymmdd', 'HHMMSSS.SS', 'MW NWR'], 
  axis = 1
  )
# create a new csv
df.to_csv(
  'earthquake_data_new.csv', 
  index = False
  )

#======================================#
# Earthquake data processing 2
# Calculate earthquake rate
# create a new column to store 
# the time from year to hour
df['hour'] = df['datetime'].str[0:10]
print(df['hour'])
dftime = df['hour'].drop_duplicates()
dftime = dftime.reset_index(drop = True)
print(dftime)

# count the number of earthquakes
earthquake_num = df.groupby(['hour']).count()['NWR'].to_numpy()
print(earthquake_num)

# create a datetime array
daterange = pd.date_range(start=datetime(2015,1,22), end=datetime(2021,3,31), freq="1H")
datestr = daterange.strftime('%Y%m%d%H')
print(datestr)

# create new dataframe to store the data
ratedf = pd.DataFrame(earthquake_num, columns = ['rate'])
ratedf['datetime'] = dftime
print(ratedf)
datedf = pd.DataFrame(datestr, columns = ['datetime'])
print(datedf)

# merge the datedf and ratedf
newdf = datedf.merge(ratedf, how = 'left', on = ['datetime'])
print(newdf)

# change all NaN to 0
newdf['rate'] = newdf['rate'].fillna(0)
print(newdf)
newdf.to_csv("overall_earthquake_rate.csv", index=False)

#======================================#
# seperate the new Dataframe in 3 data for analysis
# Before, During and After
newdf['datetime'] = pd.to_datetime(newdf['datetime'], format = '%Y%m%d%H')

# before volcano eruption
before = newdf[newdf['datetime'] < datetime(2015,4,24)]
before['datetime'] = before['datetime'].dt.strftime('%Y%m%d%H')
before.to_csv('before.csv', index = False)
# during volcano eruption
during = newdf[(newdf['datetime'] > datetime(2015,4,23)) & (newdf['datetime'] < datetime(2015,5,22))]
during['datetime'] = during['datetime'].dt.strftime('%Y%m%d%H')
during.to_csv('during.csv', index = False)
# after volcano eruption
after = newdf[newdf['datetime'] > datetime(2015,5,22)]
after['datetime'] = after['datetime'].dt.strftime('%Y%m%d%H')
after.to_csv('after.csv', index = False)

#======================================#
### End ###
#======================================#
