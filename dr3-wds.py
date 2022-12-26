# Search double stars from WDS in Gaia DR3
# Usage: python3 dr3-wds.py <WDS file>

# Import libraries

import sys
import csv
import numpy as np
from astropy.table import QTable
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import astropy.units as u
from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

# Read wds file (.csv)

filename = open(sys.argv[1], 'r')
file = csv.DictReader(filename)
fieldnames = ['wds_identifier', 'discoverer', 'component', 'epoch_first', 'epoch_last', 'observations', 'theta_first', 'theta_last', 'rho_first', 'rho_lat', 'magnitude_pri', 'magnitude_sec', 'spectral_type', 'pm_ra_pri', 'mp_dec_pri', 'pm_ra_sec', 'pm_dec_sec', 'note', 'ra_h', 'ra_m', 'ra_s', 'dec_d', 'dec_m', 'dec_s']

# Define Qtable for results
dr3Designation = np.array([], dtype=str)
dr3Ra = np.array([], dtype=np.float64)
dr3Dec = np.array([], dtype=np.float64)
dr3Parallax = np.array([], dtype=np.float64)
dr3ParallaxError = np.array([], dtype=np.float64)
dr3PmRa = np.array([], dtype=np.float64)
dr3PmDec = np.array([], dtype=np.float64)
imageId = np.array([], dtype=np.int32)
sourceId = np.array([], dtype=np.int64)
gMag = np.array([], dtype=np.float64)
wdsIdentifier = np.array([], dtype=str)
discoverer = np.array([], dtype=str)
magnitudePri = np.array([], dtype=np.float64)
sourceTable = QTable([wdsIdentifier, discoverer, dr3Designation, sourceId, dr3Ra, dr3Dec, dr3Parallax, dr3ParallaxError, dr3PmRa, dr3PmDec, magnitudePri, gMag], names=('wds_identifier', 'discoverer', 'designation', 'source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'magnitude_pri', 'phot_g_mean_mag'), meta={'name': 'first table'})

counter = int()
tableFileName = (str(sys.argv[1][:-4] + '_dr3.csv'))

for row in file:
    counter = counter + 1
    ra = row['ra_h'] + 'h' + row['ra_m'] + 'm' + row['ra_s'] + 's'
    dec = row['dec_d'] + 'd' + row['dec_m'] + 'm' + row['dec_s'] + 's'
    coord = SkyCoord(ra=Angle(ra), dec=Angle(dec), unit=(u.degree, u.degree), frame='icrs')
    radius = u.Quantity(0.0005, u.deg)
    j = Gaia.cone_search_async(coord, radius, table_name='gaiadr3.gaia_source')
    r = j.get_results()
    # Print findings
    print('\n#: ', counter, 'WDS Identifier: ' + row['wds_identifier'], 'Discoverer: ' + row['discoverer'], 'Mag Pri: ' + row['magnitude_pri'], 'RA: ' + ra, ' DEC: ', dec)
    print(r['DESIGNATION', 'source_id', 'phot_g_mean_mag', 'parallax', 'parallax_error', 'ra', 'dec', 'pmra', 'pmdec'])
    # Write findings to consolidated table
    if r:
        sourceTable.add_row([row['wds_identifier'], row['discoverer'], r[0]['DESIGNATION'], r[0]['source_id'], r[0]['ra'], r[0]['dec'], r[0]['parallax'], r[0]['parallax_error'], r[0]['pmra'], r[0]['pmdec'], row['magnitude_pri'], r[0]['phot_g_mean_mag']])
        sourceTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',')

#print(sourceTable)