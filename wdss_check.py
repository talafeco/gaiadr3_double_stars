from dsmodules import dscalculation
from astropy.table import Table, vstack, hstack
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy import units as u
import numpy as np
import sys
import datetime


# Insert the downloaded wdss file path here
wdss_file = sys.argv[1]

coords_ra = float(sys.argv[2])
coords_dec = float(sys.argv[3])
target_coords = SkyCoord(ra=coords_ra * u.deg, dec=coords_dec * u.deg)

print('### Script start! ###')
print('# Timestamp: ' + str(datetime.datetime.now()))

# Create WDS table
wdss_converters = {  'WDSS identifier': np.str_,
                    'Component identifier': np.str_,
                    'First/last observation': np.str_,
                    'Number of astrometric observations': np.str_,
                    'Position angle': np.float64,
                    'Sepatarion': np.float64,
                    'Flag for separation units': np.str_,
                    'Magnitude': np.str_,
                    'G Filter': np.str_,
                    'Infrared magnitude': np.str_,
                    'Filter': np.str_,
                    'Spectral type': np.str_,
                    'Proper motion': np.str_,
                    'Parallax': np.str_,
                    'Alternate name': np.str_,
                    'Note flags': np.str_,
                    'Coord (RA)': np.str_,
                    'Coord (DEC)': np.str_,
                    'Designation in main WDS': np.str_,
                    'Discoverer designation': np.str_,
                    'Component designation': np.str_
                    }

wdss_data = Table.read(wdss_file,
                      names=('WDSS identifier',
                    'Component identifier',
                    'First/last observation',
                    'Number of astrometric observations',
                    'Position angle',
                    'Sepatarion',
                    'Flag for separation units',
                    'Magnitude',
                    'G Filter',
                    'Infrared magnitude',
                    'Filter',
                    'Spectral type',
                    'Proper motion',
                    'Parallax',
                    'Alternate name',
                    'Note flags',
                    'Coord (RA)',
                    'Coord (DEC)',
                    'Designation in main WDS',
                    'Discoverer designation',
                    'Component designation'
                        ),
                      converters=wdss_converters,
                      format='ascii.fixed_width',
                      header_start=0, data_start=0,
                      col_starts=(1, 15, 24, 29, 33, 37, 43, 45, 50, 52, 57, 59, 65, 82, 90, 115, 118, 127, 137, 148, 155),
                       col_ends=(14, 18, 28, 32, 36, 42, 44, 50, 51, 57, 58, 64, 81, 89, 114, 117, 126, 136, 147, 155, 160),
                      )


print('Creating wdss table.')
print('# Timestamp: ' + str(datetime.datetime.now()))
wdssTable = dscalculation.create_wdss_table(wdss_data)

print('Creating wdss catalog.')
print('# Timestamp: ' + str(datetime.datetime.now()))
wdss_catalog = SkyCoord(ra=wdssTable['Coord (RA) hms'], dec=Angle(wdssTable['Coord (DEC) dms']), unit='hour, degree', frame="icrs") #

# Folyt. köv. be kell olvasni a targeteket és megkeresni az adatbázisban a legközelebbi objektumokat, aztán kiíratni, hogy mi a helyzet

print('### Script end! ###')
print('# Timestamp: ' + str(datetime.datetime.now()))

#print(wdss_catalog[0])