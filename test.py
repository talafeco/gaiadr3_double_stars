import os
import sys
import numpy as np
import math
import sys
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack, hstack
from astropy.table import Column, MaskedColumn
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
import warnings
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta
import ligo.skymap.plot

warnings.filterwarnings("ignore")

# Set working directory to read Double Star images
workingDirectory = "C:/Astro/testing/aligned"
#file = 'AVG_Aligned_9x10s.fits'
file = 'Gaia_DR3_4437183082038102528_Light_2023-07-27T22-11-00_010.new'
image_limit = 2000

#star_a_ra, star_a_dec = 327.41204013537777, 9.238097636554508
#star_b_ra, star_b_dec = 327.4110385552747, 9.225067236213862

star_a_ra, star_a_dec = 243.7383333, 3.8469444
star_b_ra, star_b_dec = 243.6287500, 3.8703889






ax = plt.axes(projection='astro zoom',
              center='5h -32d', radius='5 deg', rotate='20 deg')
ax.grid()
plt.imshow(ax) # , cmap='cividis'
plt.show()

"""# Create Image plot of the double stars
def imagePlot(filename, pairname, raa, deca, rab, decb):
    image_data = fits.open(workingDirectory + '/' + filename)
    header = image_data[0].header
    wcs_helix = WCS(image_data[0].header)
    image = image_data[0].data
    image_height = header['NAXIS2']

    star_a = SkyCoord(raa * u.deg, deca * u.deg, frame='icrs')
    star_b = SkyCoord(rab * u.deg, decb * u.deg, frame='icrs')
    star_a_pix = utils.skycoord_to_pixel(star_a, wcs_helix)
    star_b_pix = utils.skycoord_to_pixel(star_b, wcs_helix)
    plt.subplot(projection=wcs_helix)

    plt.scatter(star_a_pix[0] + (image_height / 30), star_a_pix[1], marker="_", s=(image_height / 10), color="grey")
    plt.scatter(star_a_pix[0], star_a_pix[1] + (image_height / 30), marker="|", s=(image_height / 10), color="grey")
    plt.scatter(star_b_pix[0] + (image_height / 30), star_b_pix[1], marker="_", s=(image_height / 10), color="grey")
    plt.scatter(star_b_pix[0], star_b_pix[1] + (image_height / 30), marker="|", s=(image_height / 10), color="grey")
    '''overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='grey', ls='dotted')
    overlay['ra'].set_ticklabel_position('rb')
    overlay['ra'].set_auto_axislabel(True)
    overlay['ra'].set_ticks(number=6)
    overlay['dec'].set_ticklabel_position('rb')
    overlay['dec'].set_auto_axislabel(True)
    overlay['dec'].set_ticks(number=6)'''
    plt.title(pairname)
    plt.imshow(image, origin='lower',cmap='grey', aspect='equal', vmax=image_limit, vmin=0) # , cmap='cividis'
    plt.savefig(workingDirectory + '/' + pairname + '_img.jpg',dpi=150.0, bbox_inches='tight', pad_inches=0.2)
    #plt.show()

imagePlot(file, 'Test_double_123AB', star_a_ra, star_a_dec, star_b_ra, star_b_dec)"""