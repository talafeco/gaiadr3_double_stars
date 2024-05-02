#! /usr/bin/python3

import os
import sys
import csv
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils import DAOStarFinder, aperture_photometry, CircularAperture, CircularAnnulus
from astroquery.gaia import Gaia



def find_stars(fits_file):
    hdul = fits.open(fits_file)
    wcs = WCS(hdul[0].header)
    data = hdul[0].data
    hdul.close()
    
    # Find stars using DAOStarFinder
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*np.std(data))
    sources = daofind(data)
    
    # Get sky coordinates of stars
    sky_coords = wcs.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
    
    return sky_coords

def retrieve_gaia_info(sky_coords):
    gaia_info = Gaia.query_object_async(coordinate=sky_coords, width=1*u.arcsecond, height=1*u.arcsecond)
    return gaia_info

def perform_aperture_photometry(fits_file, sky_coords):
    hdul = fits.open(fits_file)
    data = hdul[0].data
    hdul.close()
    
    positions = [(coord.ra.deg, coord.dec.deg) for coord in sky_coords]
    apertures = CircularAperture(positions, r=4.)
    annulus_apertures = CircularAnnulus(positions, r_in=10., r_out=15.)
    annulus_masks = annulus_apertures.to_mask(method='center')
    phot_table = aperture_photometry(data, apertures)
    for aperture, annulus_mask in zip(apertures, annulus_masks):
        annulus_data = annulus_mask.multiply(data)
        annulus_data_1d = annulus_data[annulus_mask.data > 0]
        _, median, _ = sigma_clipped_stats(annulus_data_1d)
        phot_table['annulus_median'] = median
        bkg_mean = median
        bkg_sum = bkg_mean * aperture.area
        phot_table['aper_bkg'] = bkg_sum
        phot_table['aper_sum_bkgsub'] = phot_table['aperture_sum'] - phot_table['aper_bkg']
    
    return phot_table

def main(folder_path, output_csv):
    # Find all FITS files in the folder
    fits_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.new')]
    
    # Open CSV file for writing
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Gaia dr3 designation', 'Gaia dr3 ra coordinate', 'Gaia dr3 deg coordinate',
                      'Measured star ra coordinate', 'Measured star deg coordinate',
                      'Measured coordinate difference (arcsec)', 'Gaia dr3 g magnitude',
                      'Measured star magnitude', 'Magnitude difference']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        # Process each FITS file
        for fits_file in fits_files:
            # Find stars in the FITS file
            sky_coords = find_stars(fits_file)
            
            # Retrieve Gaia DR3 information
            gaia_info = retrieve_gaia_info(sky_coords)
            
            # Perform aperture photometry
            phot_table = perform_aperture_photometry(fits_file, sky_coords)
            
            # Match Gaia DR3 and measured stars
            for i, gaia_row in enumerate(gaia_info):
                gaia_coord = SkyCoord(ra=gaia_row['ra'], dec=gaia_row['dec'], unit=(u.deg, u.deg))
                separation = sky_coords[i].separation(gaia_coord).arcsec
                closest_match_idx = np.argmin(separation)
                matched_star = sky_coords[closest_match_idx]
                
                # Write data to CSV
                writer.writerow({'Gaia dr3 designation': gaia_row['designation'],
                                 'Gaia dr3 ra coordinate': gaia_row['ra'],
                                 'Gaia dr3 deg coordinate': gaia_row['dec'],
                                 'Measured star ra coordinate': matched_star.ra.deg,
                                 'Measured star deg coordinate': matched_star.dec.deg,
                                 'Measured coordinate difference (arcsec)': separation[closest_match_idx],
                                 'Gaia dr3 g magnitude': gaia_row['phot_g_mean_mag'],
                                 'Measured star magnitude': phot_table['aper_sum_bkgsub'][closest_match_idx],
                                 'Magnitude difference': phot_table['aper_sum_bkgsub'][closest_match_idx] - gaia_row['phot_g_mean_mag']})

# Example usage
folder_path = sys.argv[1]
output_csv = 'stars_magnitudes.csv'

main(folder_path, output_csv)
