from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture
from astroquery.mpc import MPC
from astropy.stats import sigma_clipped_stats
import astropy.units as u
import numpy as np

# Step 1: Load the FITS image
fits_file = 'your_image.fits'
with fits.open(fits_file) as hdul:
    data = hdul[0].data
    header = hdul[0].header

# Step 2: Perform source detection
mean, median, std = sigma_clipped_stats(data, sigma=3.0)
daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
sources = daofind(data - median)

# Get pixel positions of sources
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

# Step 3: Convert pixel to world coordinates
wcs = WCS(header)
world_coords = wcs.pixel_to_world(sources['xcentroid'], sources['ycentroid'])

# Step 4: Query asteroid positions
observation_time = Time(header['DATE-OBS'])  # Ensure the header has the observation time
mpc_result = MPC.get_ephemeris('asteroid_name_or_designation', location='your_observatory_code',
                               start=observation_time, number=1)

# Parse asteroid positions
asteroid_coord = SkyCoord(ra=mpc_result['RA'], dec=mpc_result['DEC'], unit='deg')

# Step 5: Cross-match sources with asteroid positions
source_coords = SkyCoord(ra=world_coords.ra, dec=world_coords.dec)
idx, d2d, _ = source_coords.match_to_catalog_sky(asteroid_coord)

# Set a matching tolerance (e.g., 5 arcseconds)
match_tolerance = 5.0  # arcseconds
matches = d2d < match_tolerance * u.arcsec

# Print identified asteroids
if matches.any():
    print("Asteroid(s) detected:")
    print(source_coords[matches])
else:
    print("No asteroids detected.")