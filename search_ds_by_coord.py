'''
Search double stars by coordinates in: WDS, WDSS, A million binaries from Gaia catalogs
'''

import sys
import numpy as np
import math
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack
import astropy.units as u
from astropy.coordinates import Angle
import warnings
from astropy.table import QTable
from astropy.table import Table, vstack, hstack
from datetime import datetime

# Configuration
# Insert the downloaded wds file path here
wds_file = "C:\Astro\catalogs\WDS\wdsweb_summ.txt"

# Get coordinates
star_ra = float(sys.argv[1].replace(",","."))*u.degree
star_dec = float(sys.argv[2].replace(",","."))*u.degree

print('### Script start! ###')
print('# Timestamp: ' + str(datetime.now()))

def calculate_wds_ra_hourangle(wds_ra_array):
    #print(wds_ra_array)
    wds_ra_hms = []
    for star in wds_ra_array:
        wds_ra_hms.append(str(star[0:2]) + 'h' + str(star[2:4]) + 'm' + str(star[4:9]) + 's')
    return wds_ra_hms

def calculate_wds_dec_hourangle(wds_dec_array):
    #print(wds_dec_array)
    wds_dec_dms = []
    for star in wds_dec_array:
        wds_dec_dms.append(str(star[0:3]) + 'd' + str(star[3:5]) + 'm' + str(star[5:9]) + 's')
    return wds_dec_dms

def create_unique_id(wds_id, discov):
    id_array = [i + j for i, j in zip(wds_id, discov)]
    return id_array

#def delete_invalid_ines_wds(wds_catalog):
    

# Function to search coordinates in the WDS catalog file
def search_in_wds(ra_source, dec_source):
    coordinates = SkyCoord(ra=ra_source, dec=dec_source) 
    catalog = SkyCoord(ra=wds_catalog['Coord (RA) hms'], dec=Angle(wds_catalog['Coord (DEC) dms']), unit='hour, degree')
    idx, d2d, d3d = coordinates.match_to_catalog_sky(catalog)
    star_coord = SkyCoord(ra=ra_source, dec=dec_source)
    #sep = coordinates.separation(d2d)*u.degree
    print(wds_data[idx])
    print('Coordinates: ' + str(ra_source) + '(Ra), ' + str(dec_source) + '(Dec), Separation: ' + str(d2d))

def filter_catalog(catalog_name, key_colnames):
    print(str(Angle(star_ra)) + ", " + str(Angle(star_dec)))
    distance_range = '1.0d'
    coord_ra_lower = star_ra - Angle(distance_range)
    coord_ra_upper = star_ra + Angle(distance_range)
    coord_dec_lower = star_dec - Angle(distance_range)
    coord_dec_upper = star_dec + Angle(distance_range)
    print(str(coord_ra_lower) + ", " + str(coord_ra_upper) + ", " + str(coord_dec_lower) + ", " + str(coord_dec_upper))
    if ((coord_ra_lower < Angle(catalog_name['Coord (RA) hms']) < coord_ra_upper) and (coord_dec_lower < Angle(catalog_name['Coord (DEC) dms']) < coord_dec_upper)).any():
        return True
    return False
    '''
    data = catalog_name
    distance_range = '1.0d'
    coord_ra_lower = star_ra - Angle(distance_range)
    coord_ra_upper = star_ra + Angle(distance_range)
    coord_dec_lower = star_dec - Angle(distance_range)
    coord_dec_upper = star_dec + Angle(distance_range)
    #mask = (coord_ra_lower <= Angle(data['Coord (RA) hms']) <= coord_ra_upper) and (coord_dec_lower <= Angle(data['Coord (DEC) dms']) <= coord_dec_upper)
    #mask = coord_ra_lower < Angle(data['Coord (RA) hms'])
    mask = int(data['Obs']) > 20
    new_data = data[mask]
    return new_data
    '''
'''
# Function to search coordinates in the WDS catalog file
def search_in_wdss(ra_source, dec_source):
    coordinates = SkyCoord(ra=ra_source, dec=dec_source)
    temp_catalog = SkyCoord(ra=wdss_catalog['Coord (RA) hms'], dec=Angle(wdss_catalog['Coord (DEC) dms']), unit='hour, degree')
    catalog = temp_catalog.group_by('2000 Coord')
    
    idx, d2d, d3d = coordinates.match_to_catalog_sky(catalog)
    star_coord = SkyCoord(ra=ra_source, dec=dec_source)
    sep = Angle(coordinates.separation(star_coord))
    print(wdss_file[idx])
    print('Separation: ' + str(sep))

# Function to search coordinates in the WDS catalog file
def search_in_million(ra_source, dec_source):
    coordinates = SkyCoord(ra=ra_source*u.degree, dec=dec_source*u.degree)  
    catalog = SkyCoord(ra=a_millio_binaries_file['ra1']*u.degree, dec=a_millio_binaries_file['dec1']*u.degree)
    idx, d2d, d3d = coordinates.match_to_catalog_sky(catalog)
    star_coord = SkyCoord(ra=ra_source*u.degree, dec=dec_source*u.degree)
    sep = Angle(coordinates.separation(star_coord))
    print(a_millio_binaries_file[idx], sep)
'''

## Read catalogs
# WDS
#https://docs.astropy.org/en/stable/io/ascii/index.html#supported-formats
wds_converters = {  '2000 Coord': np.str_,
                    'Discov': np.str_,
                    'Comp': np.str_,
                    'Date (first)': np.str_,
                    'Date (last)': np.str_,
                    'Obs': np.str_,
                    'PA_f': np.str_,
                    'PA_l': np.str_,
                    'Sep_f': np.str_,
                    'Sep_l': np.str_,
                    'Mag_A': np.str_,
                    'Mag_B': np.str_,
                    'Spectral A/B': np.str_,
                    'PM_A_ra': np.str_,
                    'PM_A_dec': np.str_,
                    'PM_B_ra': np.str_,
                    'PM_B_dec': np.str_,
                    'D. Number': np.str_,
                    'Notes': np.str_,
                    'Coord (RA)': np.str_,
                    'Coord (DEC)': np.str_
                    }
wds_data = Table.read(wds_file,
                      names=('2000 Coord', 'Discov', 'Comp', 'Date (first)', 'Date (last)', 'Obs',
                             'PA_f', 'PA_l', 'Sep_f', 'Sep_l', 'Mag_A',
                             'Mag_B', 'Spectral A/B', 'PM_A_ra', 'PM_A_dec',
                             'PM_B_ra', 'PM_B_dec', 'D. Number', 'Notes', 'Coord (RA)',
                             'Coord (DEC)'
                        ),
                      converters=wds_converters,
                      format='ascii.fixed_width',
                      header_start=2, data_start=5,
                      col_starts=(0, 10, 17, 23, 28, 33, 38, 42, 46, 52, 58, 64, 70, 80, 84, 89, 93, 98, 107, 112, 121),
                      col_ends=(9, 16, 21, 26, 31, 36, 40, 44, 50, 56, 61, 68, 78, 83, 87, 92, 96, 105, 110, 120, 129),
                      fill_values=[('.', 'N/A', 'Coord (RA)', ), ('.', 'N/A', 'Coord (DEC)')],
                      )
#print(wds_data)
#print(wds_data.info)
print('# Timestamp, creating wds_catalog: ' + str(datetime.now()))
wds_catalog = hstack([wds_data, calculate_wds_ra_hourangle(wds_data['Coord (RA)'])])
wds_catalog.rename_column('col0', 'Coord (RA) hms')
wds_catalog = hstack([wds_catalog, calculate_wds_dec_hourangle(wds_data['Coord (DEC)'])])
wds_catalog.rename_column('col0', 'Coord (DEC) dms')
wds_catalog = hstack([wds_catalog, create_unique_id(wds_data['2000 Coord'], wds_data['Discov'])])
wds_catalog.rename_column('col0', 'Unique ID')
print(wds_catalog)
print(wds_catalog.info)
print('# Timestamp, start wds catalog filter: ' + str(datetime.now()))
#wds_catalog.group_by('Unique ID')
#wds_catalog.groups.filter(filter_catalog)
print('# Timestamp, wds catalog filter complete: ' + str(datetime.now()))

#wds_catalog['Coord (RA) hms'].mask = [True]
print('# Timestamp, start coordinate search in wds catalog: ' + str(datetime.now()))
search_in_wds(star_ra, star_dec)
print('# Timestamp, search in wds catalog finished: ' + str(datetime.now()))


'''
# WDSS
# https://docs.astropy.org/en/stable/io/ascii/index.html#supported-formats
wdss_converters = {'WDSS id': np.str_, 
                   'Comp': np.str_, 
                   'Date': np.str_, 
                   'Obs': np.str_, 
                   'PA': np.str_, 
                   'Sep': np.str_,
                    'Sep unit': np.str_, 
                    'Mag vis': np.str_, 
                    'Filter vis': np.str_, 
                    'Mag inf': np.str_, 
                    'Filter inf': np.str_,
                    'Spectral': np.str_, 
                    'PM RA': np.str_, 
                    'PM DEC': np.str_, 
                    'Parallax': np.str_, 
                    'Alternate name': np.str_,
                    'Notes': np.str_, 
                    'Coord (RA)': np.str_, 
                    'Coord (DEC)': np.str_, 
                    'WDS designation': np.str_, 
                    'Discoverer': np.str_, 
                    'Component': np.str_    
                    }
print('# Timestamp, reading wdss file: ' + str(datetime.now()))
wdss_file = Table.read(f"~/Library/astro/catalogs/wds/wdss/wdss_summ.txt",
                      guess=False,
                      converters=wdss_converters,
                      names=('WDSS id', 'Comp', 'Date', 'Obs', 'PA', 'Sep',
                             'Sep unit', 'Mag vis', 'Filter vis', 'Mag inf', 'Filter inf',
                             'Spectral', 'PM RA', 'PM DEC', 'Parallax', 'Alternate name',
                             'Notes', 'Coord (RA)', 'Coord (DEC)', 'WDS designation', 'Discoverer', 'Component'
                            ), 
                      format='ascii.fixed_width_no_header',
                      data_start=0,
                      col_starts=(0, 15, 24, 29, 33, 37,
                                  43, 45, 50, 52, 57,
                                  59, 65, 73, 82, 90,
                                  115, 118, 127, 137, 148, 155),
                      col_ends=(13, 17, 27, 31, 35, 42,
                                43, 49, 50, 56, 57,
                                63, 72, 80, 88, 113,
                                116, 126, 136, 146, 154, 159),
                      )
#print(wdss_file)
#print(wdss_file.info)
print('# Timestamp, creating wdss_catalog: ' + str(datetime.now()))
wdss_catalog= hstack([wdss_file, calculate_wds_ra_hourangle(wdss_file['Coord (RA)'])])
wdss_catalog.rename_column('col0', 'Coord (RA) hms')
wdss_catalog= hstack([wdss_catalog, calculate_wds_dec_hourangle(wdss_file['Coord (DEC)'])])
wdss_catalog.rename_column('col0', 'Coord (DEC) dms')
print(wdss_catalog)
print(wdss_catalog.info)
print('# Timestamp, start coordinate search in wdss catalog: ' + str(datetime.now()))
#search_in_wdss(star_ra, star_dec)

print('# Timestamp, search in wdss catalog finished: ' + str(datetime.now()))
s
# A million binaries from Gaia

a_million_converters = {'source_id1': int64,
                        'source_id2': int64,
                        'ra1': float64,
                        'ra2': float64,
                        'dec1': float64,
                        'dec2': float64,
                        'parallax1': float64,
                        'parallax2': float64,
                        'parallax_error1': float64,
                        'parallax_error2': float64,
                        'pm1':float64,
                        'pm2': float64,
                        'pmra1': float64,
                        'pmra2': float64,
                        'pmra_error1': float64,
                        'pmra_error2': float64,
                        'pmdec1': float64,
                        'pmdec2': float64,
                        'pmdec_error1': float64,
                        'pmdec_error2': float64,
                        'phot_g_mean_mag1': float64,
                        'phot_g_mean_mag2': float64,
                        'phot_bp_mean_mag1': float64,
                        'phot_bp_mean_mag2': float64,
                        'phot_rp_mean_mag1': float64,
                        'phot_rp_mean_mag2': float64,
                        'pairdistance': float64,
                        'sep_AU': float64,
                        'binary_type': str4,
                        'dr2_source_id1': int64,
                        'dr2_source_id1_1': int64,
                    }

a_millio_binaries_file = Table.read(f"~/Library/astro/catalogs/wds/2023-05-12/a_million_binaries_from_gaia_erf3.csv",
                                    format='csv',
                                    header_start=0,
                                    data_start=1)

print(a_millio_binaries_file)
print(a_millio_binaries_file.info)
search_in_million(star_ra, star_dec)







for star in sources:
        ra2, dec2 = mywcs.all_pix2world([[star ['xcentroid'], star ['ycentroid']]], 0)[0]   
        c = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
        #catalog = SkyCoord(ra=gaiaStars[1:, 5]*u.degree, dec=gaiaStars[1:, 7]*u.degree)  
        catalog = SkyCoord(ra=gaiaStars['ra']*u.degree, dec=gaiaStars['dec']*u.degree)
        idx, d2d, d3d = c.match_to_catalog_sky(catalog)
        catalogstar = SkyCoord(ra=gaiaStars[idx]['ra']*u.degree, dec=gaiaStars[idx]['dec']*u.degree)
        sep = c.separation(catalogstar)
        if sep < Angle('00d00m02s'):
            sourceTable.add_row([fitsFile, gaiaStars[idx]['source_id'], gaiaStars[idx]['designation'], convertStringToNan(gaiaStars[idx]['ra']), convertStringToNan(gaiaStars[idx]['dec']), convertStringToNan(gaiaStars[idx]['parallax']), convertStringToNan(gaiaStars[idx]['parallax_error']), convertStringToNan(gaiaStars[idx]['pmra']), convertStringToNan(gaiaStars[idx]['pmdec']), convertStringToNan(gaiaStars[idx]['phot_g_mean_mag']), convertStringToNan(gaiaStars[idx]['phot_bp_mean_mag']), convertStringToNan(gaiaStars[idx]['phot_rp_mean_mag']), convertStringToNan(gaiaStars[idx]['radial_velocity']), convertStringToNan(gaiaStars[idx]['radial_velocity_error']), convertStringToNan(gaiaStars[idx]['teff_gspphot']), star['id'], ra2, dec2, star['mag']])
'''