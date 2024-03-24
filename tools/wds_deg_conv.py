# importing the module
import sys
import csv
from astropy.table import QTable
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import numpy as np

#filename = open(sys.argv[1], 'r')
#file = csv.DictReader(filename)
#fieldnames = ['wds_identifier', 'discovr', 'comp', 'theta', 'rho', 'mag_pri', 'mag_sec', 'spectra', 'pm_a_ra', 'pm_a_dec', 'pm_b_ra', 'pm_b_dec', 'ra_hrs', 'dec']

print('\n### Reading WDS database ###')
doubleStars = np.genfromtxt(sys.argv[1], delimiter=',', skip_header=1, dtype=str)

#starList = []
#for line in file:
#    starList.append(line)

#print(starList)
print('\n### Creating Qtable in the memory ###')
# Define Qtable for results
wdsIdentifier = np.array([], dtype=str)
discoverer = np.array([], dtype=str)
components = np.array([], dtype=str)
wdsTheta = np.array([], dtype=str)
wdsRho = np.array([], dtype=str)
magPri = np.array([], dtype=str)
magSec = np.array([], dtype=str)
wdsSpectra = np.array([], dtype=str)
pmARa = np.array([], dtype=str)
pmADec = np.array([], dtype=str)
pmBRa = np.array([], dtype=str)
pmBDec = np.array([], dtype=str)
raSysHrs = np.array([], dtype=str)
decSys = np.array([], dtype=str)
raDeg = np.array([], dtype=np.float64)
decDeg = np.array([], dtype=np.float64)
sourceTable = QTable([wdsIdentifier, discoverer, components, wdsTheta, wdsRho, magPri, magSec, wdsSpectra, pmARa, pmADec, pmBRa, pmBDec, raSysHrs, decSys, raDeg, decDeg], names=('wds_identifier', 'discovr', 'comp', 'theta', 'rho', 'mag_pri', 'mag_sec', 'spectra', 'pm_a_ra', 'pm_a_dec', 'pm_b_ra', 'pm_b_dec', 'ra_hms', 'dec_dms', 'ra_deg', 'dec_deg'), meta={'name': 'first table'})

print('\n### Start converting hourangles to degrees ###')
counter = ()

for line in doubleStars:
    #print(['ra_hms'])
    #print(line['dec'])
    raLine = (str(line[12][0:2]) + 'h' + str(line[12][2:4]) + 'm' + str(line[12][4:9]) + 's')
    decLine = line[13][0:3] + 'd' + line[13][3:5] + 'm' + line[13][5:9] + 's'
    #print(Angle(raLine), decLine)
    coord = SkyCoord(raLine, decLine, unit=(u.hourangle, u.deg), frame='icrs') # unit=(u.hourangle, u.deg)
    #print(line[0], line[1], coord)
    sourceTable.add_row([line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], coord.ra, coord.dec.degree])

#print(sourceTable)

tableFileName = (str(sys.argv[1][:-4] + '.dat'))
sourceTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',')