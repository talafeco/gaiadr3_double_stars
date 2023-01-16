# importing the module
import sys
import csv
from astropy.table import QTable
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import numpy as np

filename = open(sys.argv[1], 'r')
file = csv.DictReader(filename)
fieldnames = ['wds_identifier', 'discovr', 'comp', 'theta', 'rho', 'mag_pri', 'mag_sec', 'spectra', 'pm_a_ra', 'pm_a_dec', 'pm_b_ra', 'pm_b_dec', 'ra_hrs', 'dec']

starList = []
for line in file:
    starList.append(line)

#print(starList)

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
sourceTable = QTable([wdsIdentifier, discoverer, components, wdsTheta, wdsRho, magPri, magSec, wdsSpectra, pmARa, pmADec, pmBRa, pmBDec, raSysHrs, decSys, raDeg, decDeg], names=('wds_identifier', 'discovr', 'comp', 'theta', 'rho', 'mag_pri', 'mag_sec', 'spectra', 'pm_a_ra', 'pm_a_dec', 'pm_b_ra', 'pm_b_dec', 'ra_hrs', 'dec_deg', 'ra', 'dec'), meta={'name': 'first table'})

for line in starList:
    print(line['ra_hrs'])
    #print(line['dec'])
    raLine = (str(line['ra_hrs'][0:2]) + 'h' + str(line['ra_hrs'][2:4]) + 'm' + str(line['ra_hrs'][4:9]) + 's')
    decLine = line['dec'][0:3] + 'd' + line['dec'][3:5] + 'm' + line['dec'][5:9] + 's'
    #print(Angle(raLine), decLine)
    coord = SkyCoord(raLine, decLine, unit=(u.hourangle, u.deg), frame='icrs') # unit=(u.hourangle, u.deg)
    print(line['wds_identifier'], line['discovr'], coord)
    sourceTable.add_row([line['wds_identifier'], line['discovr'], line['comp'], line['theta'], line['rho'], line['mag_pri'], line['mag_sec'], line['spectra'], line['pm_a_ra'], line['pm_a_dec'], line['pm_b_ra'], line['pm_b_dec'], line['ra_hrs'], line['dec'], coord.ra, coord.dec.degree])

#print(sourceTable)

tableFileName = (str(sys.argv[1][:-4] + '.dat'))
sourceTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',')