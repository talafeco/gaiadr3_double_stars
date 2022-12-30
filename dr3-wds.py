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
magnitudePri = np.array([], dtype=str)
sourceTable = QTable([wdsIdentifier, discoverer, dr3Designation, sourceId, dr3Ra, dr3Dec, dr3Parallax, dr3ParallaxError, dr3PmRa, dr3PmDec, magnitudePri, gMag], names=('wds_identifier', 'discoverer', 'designation', 'source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'magnitude_pri', 'phot_g_mean_mag'), meta={'name': 'first table'})

counter = int()
importCounter = int()
tableFileName = (str(sys.argv[1][:-4] + '_dr3.csv'))

for row in file:
    counter = counter + 1
    importCounter = importCounter + 1
    ra = row['ra_h'] + 'h' + row['ra_m'] + 'm' + row['ra_s'] + 's'
    dec = row['dec_d'] + 'd' + row['dec_m'] + 'm' + row['dec_s'] + 's'
    coord = SkyCoord(ra=Angle(ra), dec=Angle(dec), unit=(u.degree, u.degree), frame='icrs', equinox='J2000.000')
    radius = u.Quantity(0.001, u.deg)
    j = Gaia.cone_search(coord, radius, table_name='gaiadr3.gaia_source')
    r = j.get_results()
    # Print findings
    print('\n#: ', counter, 'WDS Identifier: ' + row['wds_identifier'], 'Discoverer: ' + row['discoverer'], 'Mag Pri: ' + row['magnitude_pri'], 'RA: ' + ra, ' DEC: ', dec)
    print(r['DESIGNATION', 'source_id', 'phot_g_mean_mag', 'parallax', 'parallax_error', 'ra', 'dec', 'pmra', 'pmdec'])
    # Write findings to consolidated table
    if r:
        for line in r:
            sourceTable.add_row([row['wds_identifier'], row['discoverer'], r[line.index]['DESIGNATION'], r[line.index]['source_id'], r[line.index]['ra'], r[line.index]['dec'], r[line.index]['parallax'], r[line.index]['parallax_error'], r[line.index]['pmra'], r[line.index]['pmdec'], row['magnitude_pri'], r[line.index]['phot_g_mean_mag']])
            sourceTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',')
    elif not r:
        sourceTable.add_row([row['wds_identifier'], row['discoverer'], 'Not found in Gaia DR3', '0', '0', '0', '0', '0', '0', '0', row['magnitude_pri'], '0'])
        sourceTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',')
    

#print(sourceTable)

    #print(row['ra_h'], row['ra_m'], row['ra_s'], row['dec_d'], row['dec_m'], row['dec_s'])
    """ if row['ra_h']: #  and row['magnitude_pri']
        ra = row['ra_h'] + 'h' + row['ra_m'] + 'm' + row['ra_s'] + 's'
        dec = row['dec_d'] + 'd' + row['dec_m'] + 'm' + row['dec_s'] + 's'
        coord = SkyCoord(ra=Angle(ra), dec=Angle(dec), unit=(u.degree, u.degree), frame='icrs', equinox='J2000.000')
        radius = u.Quantity(0.001, u.deg)
    else:
        print('\n#: ', counter, 'WDS Identifier: ' + row['wds_identifier'], 'Discoverer: ' + row['discoverer'], 'Has no precise coodinates!')
        sourceTable.add_row([row['wds_identifier'], row['discoverer'], 'Precise coodinates or magnitude is missing!', '0', '0', '0', '0', '0', '0', '0', '0', '0'])
        sourceTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',') """