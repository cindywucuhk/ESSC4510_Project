# ESSC4510_Project
ESSC 4510 group project
Question 7
Download the catalog of earthquakes happening at a submarine volcano offshore the
west coast of the USA. Are there any
temporal trends (e.g. periodicity) in the earthquake rate (number/time)? The volcano
erupted from April 24th to May 21st 2015, so you might want to compare your results
from analyzing data from before, during vs after the eruption. How does the
earthquake rate temporal periodicity compare with the periodicity in ocean height
variation at this location due to tides (data will be provided)?

# Powerpoint
link: https://docs.google.com/presentation/d/1cRwa86Hu6ojwzb0phgGCeJlrYw7PBkfAfeG-m0L8qUY/edit?usp=sharing

# Data
Submarine volcano offshore the west coast of the USA

Axial Seamount Earthquake Catalog

http://axial.ocean.washington.edu/#name6

Tidal data of the west coast of the USA

https://www.ndbc.noaa.gov/station_page.php?station=46404

### Data lists
- earthquake_data_new.csv
  - New formatted earthquake data
- earthquake_data.csv
  - Retreived directly from the URL above
- before.csv
  - Earthquake rate before the volcanic eruption
- during.csv
  - Earthquake rate during the volcanic eruption
- after.csv
  - Earthquake rate after the volcanic eruption

# Python files
Here shows all the python files and their usage.
- plot_raw.py
  - plotting earthquake data
  - make the new formatted earthquake data
- seperate_data.py
  - calculate the earthquake rate per hour
  - seperate the data into before, during and after
  - plot earthquake rate
- periodicity_before.py
  - time series, power/amplitude spectral density for before the volcanic eruption
- periodicity_during.py
  - time series, power/amplitude spectral density for during the volcanic eruption
- periodicity_after.py
  - time series, power/amplitude spectral density for after the volcanic eruption
- periodicity_all.py
  - time series, power/amplitude spectral density for the whole set of data
- tide.py
  - Tide data processing, missing data processing, Fourier transform and power spectrum plotting

# Reference
Wilcock, W. S. D., M. Tolstoy, F. Waldhauser, C. Garcia, Y. J. Tan, D. R. Bohnenstiehl, J. Caplan-Auerbach, R. P. Dziak, A. Arnulf, & M. E. Mann (2016). Seismic constraints on caldera dynamics from the 2015 Axial Seamount eruption, Science, 354, 1395-1399. and

Wilcock, W. S. D., F. Waldhauser, & M. Tolstoy (2017). Catalogs of earthquake recorded on Axial Seamount from January, 2015 through November, 2015 (investigators William Wilcock, Maya Tolstoy, Felix Waldhauser). Interdisciplinary Earth Data Alliance. https://doi.org/10.1594/IEDA/323843.
