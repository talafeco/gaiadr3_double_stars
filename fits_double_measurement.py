import sys
import csv
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.datasets import load_star_image
from photutils.detection import DAOStarFinder
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture
from astropy.wcs import WCS

import astropy.units as u
from astropy.coordinates import SkyCoord
#from astroquery.gaia import Gaia
#Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"


# Read fits file
fitsFileName = ('sao_127029.fits')
hdu = fits.open(fitsFileName)
mywcs = WCS(hdu[0].header)
#hdu.info()

# Estimate the background and background noise
data = hdu[0].data
mean, median, std = sigma_clipped_stats(data, sigma=5.0)  
#print((mean, median, std))  

daofind = DAOStarFinder(fwhm=10.0, threshold=26.0*std)  
sources = daofind(data - median)  
""" for col in sources.colnames:  
    sources[col].info.format = '%.8g'  # for consistent table output """
#print(type(sources))
# print(sources)

# Iterate trhough table elements
for row in sources:
    ra, dec = mywcs.all_pix2world([[row ['xcentroid'], row ['ycentroid']]], 0)[0]
    print(row['id'], ra, dec, row['mag'])

""" !Problem: only Gaia DR2 data is available..
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    radius = u.Quantity(0.001, u.deg)
    j = Gaia.cone_search_async(coord, radius)
    r = j.get_results()

    print(r['DESIGNATION'], row['id'], ra, dec, row['mag']) """

# Plot image sources

""" positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=14.)
# norm = ImageNormalize(stretch=SqrtStretch())
plt.figure(dpi=300)
plt.imshow(data, cmap='Greys',origin='lower', interpolation='nearest')
apertures.plot(color='blue', lw=0.5, alpha=0.5)
plt.savefig('sao_127029.png') """

# hdu.close()

# Convert pixel to coordinate
"""Example conversion from pixel coordinates to ra/dec."""
""" from astropy.io import fits
from astropy.wcs import WCS

x, y = 250, 250

f = fits.open('f242_e0_c98.fits')
mywcs = WCS(f[0].header)
ra, dec = mywcs.all_pix2world([[x, y]], 0)[0] """
