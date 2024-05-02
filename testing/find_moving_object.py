'''
This function find_moving_objects takes two FITS files, representing earlier and later images, along with their acquisition times as input. It then performs the subtraction of the later-generated FITS file from the earlier one, identifies sources above a certain threshold in the later image, calculates their coordinates, measures the angular distance traveled by the sources in the sky, and finally computes the speed of movement in arcseconds per day.
'''

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.time import Time

def find_moving_objects(file_earlier, file_later, time_earlier, time_later, threshold_sigma=3):
    # Open FITS files and obtain data
    data_earlier, header_earlier = fits.getdata(file_earlier, header=True)
    data_later, header_later = fits.getdata(file_later, header=True)

    # Get WCS objects
    wcs_earlier = WCS(header_earlier)
    wcs_later = WCS(header_later)

    # Calculate image statistics for later image
    mean, median, std = sigma_clipped_stats(data_later, sigma=3.0)

    # Find sources in the later image above a certain threshold
    sources_later = data_later > threshold_sigma * std

    # Get coordinates of sources in the later image
    y_sources_later, x_sources_later = np.where(sources_later)
    world_coords_later = wcs_later.all_pix2world(np.column_stack((x_sources_later, y_sources_later)), 1)
    sky_coords_later = SkyCoord(world_coords_later, unit='deg')

    # Transform time strings into Astropy Time objects
    time_earlier_astropy = Time(time_earlier)
    time_later_astropy = Time(time_later)

    # Calculate time difference in days
    time_difference_days = (time_later_astropy - time_earlier_astropy).jd

    # Calculate the angular distance traveled by sources in the sky
    angular_distance_arcsec = sky_coords_later.separation(sky_coords_earlier).arcsec

    # Calculate speed in arcseconds per day
    speed_arcsec_per_day = angular_distance_arcsec / time_difference_days

    return speed_arcsec_per_day

# Example usage:
file_earlier = 'earlier_image.fits'
file_later = 'later_image.fits'
time_earlier = '2024-04-22T12:00:00'  # Time of earlier image acquisition
time_later = '2024-04-23T12:00:00'    # Time of later image acquisition

speed_arcsec_per_day = find_moving_objects(file_earlier, file_later, time_earlier, time_later)
print("Speed of movement of objects in arcseconds per day:", speed_arcsec_per_day)
