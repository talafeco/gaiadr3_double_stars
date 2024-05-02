#!/usr/local/bin/python3

# Original file can be found:
# https://gist.github.com/dokeeffe/18ae5c208aaac7aad1de42a411701919

import sys
import requests, math
import pandas as pd
import numpy as np
from photutils import DAOStarFinder
from astropy.stats import mad_std
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture
#%matplotlib inline
#plt.style.use('seaborn')

def get_comp_stars(ra,dec,filter_band='V',field_of_view=7.5):
    result = []
    vsp_template = 'https://www.aavso.org/apps/vsp/api/chart/?format=json&fov=' + str(field_of_view) + '&maglimit=18.5&ra=' + str(ra) + '&dec=' + str(dec)
    print(vsp_template)
    r = requests.get('https://app.aavso.org/vsp/api/chart/?format=json&fov=7.5&maglimit=20.5&ra=21%3A49%3A38.89&dec=09%3A14%3A17.4')
    #print('Downloaded Comparison Star Chart ID: '.format(r.json()['chartid']))
    for star in r.json()['photometry']:
        comparison = {}
        comparison['auid'] = star['auid']
        comparison['ra'] = star['ra']
        comparison['dec'] = star['dec']
        for band in star['bands']:
            if band['band'] == filter_band:
                comparison['vmag'] = band['mag']
                comparison['error'] = band['error']
        result.append(comparison)
    return result

TARGET_RA = '21:49:38.8897'
TARGET_DEC = '+9:14:17.410'
FITS_DATA_FILE = sys.argv[1]

# Source detection and photometry settings
FWHM = 7.0
SOURCE_SNR = 20
APERTURE_RADIUS = 7.0

comp_stars = get_comp_stars(TARGET_RA, TARGET_DEC)

hdulist = fits.open(FITS_DATA_FILE)
data = hdulist[0].data.astype(float)
wcs = WCS(hdulist[0].header)
bkg_sigma = mad_std(data)    
daofind = DAOStarFinder(fwhm=FWHM, threshold=SOURCE_SNR*bkg_sigma)    
sources = daofind(data)
sources

positions = (sources['xcentroid'], sources['ycentroid'])    
apertures = CircularAperture(positions, r=APERTURE_RADIUS)    
phot_table = aperture_photometry(data, apertures)    
print(phot_table)  

target = SkyCoord('08 10 56.65','28 08 33.2', unit=(u.hourangle, u.deg))
target_xy = SkyCoord.to_pixel(target, wcs=wcs, origin=1)
target_xy

comp_stars.append({'auid': 'target', 'ra': TARGET_RA, 'dec': TARGET_DEC})

for comp_star in comp_stars:
    comp_coord = SkyCoord(comp_star['ra'],comp_star['dec'], unit=(u.hourangle, u.deg))
    xy = SkyCoord.to_pixel(comp_coord, wcs=wcs, origin=1)
    x = xy[0].item(0)
    y = xy[1].item(0)
    for phot_measurement in phot_table:
        if(phot_measurement['xcenter'].value-4 < x < phot_measurement['xcenter'].value+4) and phot_measurement['ycenter'].value-4 < y < phot_measurement['ycenter'].value+4:
            comp_star['x'] = x
            comp_star['y'] = y
            comp_star['measured_flux'] = phot_measurement['aperture_sum']
            if not math.isnan(phot_measurement['aperture_sum']):
                comp_star['measured_instrumental_mag'] = -2.5 * math.log10(phot_measurement['aperture_sum'])
    
results = pd.DataFrame(comp_stars)
results

fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
positions = (results['x'], results['y'])    
position_apps = CircularAperture(positions, r=20.)    
target_app = CircularAperture(target_xy, r=20.)    
plt.imshow(data, cmap='gray_r', origin='lower', vmin=0, vmax=2500)
position_apps.plot(color='red', lw=1.5, alpha=0.5)
target_app.plot(color='blue', lw=1.5, alpha=0.5)

to_plot = results.query('vmag > -9999 and measured_flux > 0')
x = to_plot['measured_instrumental_mag'].values
y = to_plot['vmag'].values
fit = np.polyfit(x,y,1)
fit_fn = np.poly1d(fit) 
# fit_fn is now a function which takes in x and returns an estimate for y
fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
plt.plot(x,y, 'yo', x, fit_fn(x), '--k')
plt.plot(target_xy)
plt.xlim(np.min(x)-1, np.max(x)+1)
plt.ylim(np.min(y)-1, np.max(y)+1)

target_instrumental_magnitude = results[results.auid=='target']['measured_instrumental_mag'].values[0]
converted_magnitude = fit_fn(target_instrumental_magnitude)
converted_magnitude