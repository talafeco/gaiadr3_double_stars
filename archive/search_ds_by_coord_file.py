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
from astropy.table import Table, hstack
from datetime import datetime

# Configuration
# Insert the downloaded wds file path here
wds_file = "/usr/share/dr3map/wds/wdsweb_summ2.txt"

# Get coordinates
# star_ra = float(sys.argv[1].replace(",","."))*u.degree
# star_dec = float(sys.argv[2].replace(",","."))*u.degree

# Get coordinates
input_file = sys.argv[1]

results = []

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


def delete_invalid_lines_wds(catalog):
    rows_to_delete = np.where(catalog['Coord (RA) hms'] == '.hms')
    catalog.remove_rows(rows_to_delete)
    return catalog


# Function to search coordinates in the WDS catalog file
'''def search_in_wds(source_id, ra_source, dec_source):
    coordinates = SkyCoord(ra=ra_source, dec=dec_source, unit='degree, degree')
    catalog = SkyCoord(ra=wds_catalog['Coord (RA) hms'], dec=Angle(wds_catalog['Coord (DEC) dms']), unit='hour, degree', frame="icrs")
    idx, d2d, d3d = coordinates.match_to_catalog_sky(catalog)
    star_coord = SkyCoord(ra=ra_source, dec=dec_source, unit='degree, degree')
    print('Coordinates: ' + str(ra_source) + '(Ra), ' + str(dec_source) + '(Dec), Separation: ' + str(d2d))
    #sep = coordinates.separation(d2d)*u.degree
    print(wds_catalog[idx]['2000 Coord'])
    print(wds_catalog[np.where(wds_catalog['2000 Coord'] == wds_catalog[idx]['2000 Coord'])])'''

def search_all_in_wds(source_ids, ra_sources, dec_sources):
    print('### Starting batch match ###')
    coordinates = SkyCoord(ra=ra_sources, dec=dec_sources, unit='degree', frame='icrs')
    
    wds_coords = SkyCoord(
        ra=wds_catalog['Coord (RA) hms'],
        dec=Angle(wds_catalog['Coord (DEC) dms']),
        unit=('hour', 'degree'),
        frame='icrs'
    )

    idx, d2d, _ = coordinates.match_to_catalog_sky(wds_coords)

    matched_rows = wds_catalog[idx]

    # Build result table column-wise for performance
    results = Table()

    results['source_id'] = source_ids
    results['input_ra_deg'] = ra_sources
    results['input_dec_deg'] = dec_sources
    results['matched_unique_id'] = matched_rows['Unique ID']
    results['matched_2000_coord'] = matched_rows['2000 Coord']
    results['matched_discov'] = matched_rows['Discov']
    results['matched_comp'] = matched_rows['Comp']
    results['matched_date_last'] = matched_rows['Date (last)'].astype(str)
    results['matched_pa_last'] = matched_rows['PA_l'].astype(str)
    results['matched_sep_last'] = matched_rows['Sep_l'].astype(str)
    results['matched_mag_a'] = matched_rows['Mag_A'].astype(str)
    results['matched_mag_b'] = matched_rows['Mag_B'].astype(str)
    results['matched_coord_ra'] = matched_rows['Coord (RA)']
    results['matched_coord_dec'] = matched_rows['Coord (DEC)']
    results['matched_ra_hms'] = matched_rows['Coord (RA) hms']
    results['matched_dec_dms'] = matched_rows['Coord (DEC) dms']
    results['separation_arcsec'] = d2d.arcsecond

    return results

def write_wds_results_to_csv(results_table, output_filename='wds_search_results.csv'):
    # Optionally, convert any column types to string if needed here
    results_table.write(output_filename, format='csv', overwrite=True)
    print(f"# Saved results to {output_filename}")

 # Function to write results into a file

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
                      )

input_table_converters = {  'source_id': np.str_,
                    'ra_deg': np.float64,
                    'dec_deg': np.float64,
                    }

input_table = Table.read(input_file,
                         names=('source_id', 'ra_deg', 'dec_deg'),
                         #converters=input_table_converters,
                         format='ascii',
                         delimiter=','
                         )

print('# Timestamp, creating wds_catalog: ' + str(datetime.now()))
wds_catalog = hstack([wds_data, calculate_wds_ra_hourangle(wds_data['Coord (RA)'])])
wds_catalog.rename_column('col0', 'Coord (RA) hms')
wds_catalog = hstack([wds_catalog, calculate_wds_dec_hourangle(wds_data['Coord (DEC)'])])
wds_catalog.rename_column('col0', 'Coord (DEC) dms')
wds_catalog = hstack([wds_catalog, create_unique_id(wds_data['2000 Coord'], wds_data['Discov'])])
wds_catalog.rename_column('col0', 'Unique ID')
wds_catalog = delete_invalid_lines_wds(wds_catalog)
print('# Timestamp, start coordinate search in wds catalog: ' + str(datetime.now()))

'''for star in input_table:
    search_in_wds(star[0], star[1], star[2])'''

'''for star in input_table:
    result = search_in_wds(star[0], star[1], star[2])
    print(result)
    results.append(result)
results = search_in_wds(input_table['source_id'], input_table['ra_deg'], input_table['dec_deg'])'''

results_table = search_all_in_wds(input_table['source_id'], input_table['ra_deg'], input_table['dec_deg'])
write_wds_results_to_csv(results_table)

print('# Timestamp, search in wds catalog finished: ' + str(datetime.now()))