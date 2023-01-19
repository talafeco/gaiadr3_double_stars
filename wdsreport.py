import csv
import os
import sys
import numpy as np
import datetime
import math
import sys
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import warnings
from io import StringIO
from astropy.io import ascii
warnings.filterwarnings("ignore")


def convertStringToNan(str):
    if str == '' or str == '.':
        str = np.nan
    return str

""" def convertRa(ra):
    raline = (str(ra[0:2]) + 'h' + str(ra[2:4]) + 'm' + str(ra[4:9]) + 's')
    radeg = Angle(raline)
    return radeg.degree

def convertDec(dec):
    decline = dec[0:3] + 'd' + dec[3:5] + 'm' + dec[5:9] + 's'
    decdeg = Angle(decline)
    return decdeg """

print('\n### Reading WDS database ###')
""" doubleStars = np.genfromtxt(sys.argv[1], delimiter=',', skip_header=1, dtype=str)
print(doubleStars) """
#wdsTable = np.genfromtxt(sys.argv[1], delimiter=',')
wdsTable = Table.read(sys.argv[1], delimiter=',', format='ascii')
print(wdsTable.info)


### Define Qtable for WDS

""" wds_identifier = np.array([], dtype=str)
wdsdiscovr = np.array([], dtype=str)
wdscomp = np.array([], dtype=str)
wdstheta = np.array([], dtype=np.float64)
wdsrho = np.array([], dtype=np.float64)
wdsmag_pri = np.array([], dtype=np.float64)
wdsmag_sec = np.array([], dtype=np.float64)
wdsspectra = np.array([], dtype=str)
wdspm_a_ra = np.array([], dtype=np.float64)
wdspm_a_dec = np.array([], dtype=np.float64)
wdspm_b_ra = np.array([], dtype=np.float64)
wdspm_b_dec = np.array([], dtype=np.float64)
wdsra = np.array([], dtype=str)
wdsdec = np.array([], dtype=str)
wdsdegra = np.array([], dtype=np.float64)
wdsdegdec = np.array([], dtype=np.float64)
wdsTable = QTable([wds_identifier, wdsdiscovr, wdscomp, wdstheta, wdsrho, wdsmag_pri, wdsmag_sec, wdsspectra, wdspm_a_ra, wdspm_a_dec, wdspm_b_ra, wdspm_b_dec, wdsra, wdsdec, wdsdegra, wdsdegdec], names=('wds_identifier', 'discovr', 'comp', 'theta', 'rho', 'mag_pri', 'mag_sec', 'spectra', 'pm_a_ra', 'pm_a_dec', 'pm_b_ra', 'pm_b_dec', 'ra_hms', 'dec_dms', 'ra_deg', 'dec_deg'), meta={'name': 'wds table'}) """

dsfilename = np.array([], dtype=str)
dswds_identifier = np.array([], dtype=str)
dsdiscovr = np.array([], dtype=str)
dscomp = np.array([], dtype=str)
dstheta = np.array([], dtype=np.float64)
dsrho = np.array([], dtype=np.float64)
dsmag_pri = np.array([], dtype=np.float64)
dsmag_sec = np.array([], dtype=np.float64)
dsmag_diff = np.array([], dtype=np.float32)
dsspectra = np.array([], dtype=str)
dspm_a_ra = np.array([], dtype=np.float64)
dspm_a_dec = np.array([], dtype=np.float64)
dspm_b_ra = np.array([], dtype=np.float64)
dspm_b_dec = np.array([], dtype=np.float64)
dsra = np.array([], dtype=str)
dsdec = np.array([], dtype=str)
dsdegra = np.array([], dtype=np.float64)
dsdegdec = np.array([], dtype=np.float64)
dsobjectid = np.array([], dtype=str)
dspaactual = np.array([], dtype=np.float64)
dssepactual = np.array([], dtype=np.float64)
dsmagdiff = np.array([], dtype=np.float64)
dsTable = QTable([dswds_identifier, dsdiscovr, dscomp, dstheta, dsrho, dsmag_pri, dsmag_sec, dsmag_diff, dsspectra, dspm_a_ra, dspm_a_dec, dspm_b_ra, dspm_b_dec, dsra, dsdec, dsdegra, dsdegdec, dsobjectid, dspaactual, dssepactual, dsmagdiff], names=('wds_identifier', 'discovr', 'comp', 'theta', 'rho', 'mag_pri', 'mag_sec', 'dsmag_diff', 'spectra', 'pm_a_ra', 'pm_a_dec', 'pm_b_ra', 'pm_b_dec', 'ra_hms', 'dec_dms', 'ra_deg', 'dec_deg', 'object_id', 'dspaactual', 'dssepactual', 'dsmagdiff'), meta={'name': 'ds table'})


#print('\n### Adding WDS Double stars to numpy array ###')
""" for line in doubleStars:
    wdsTable.add_row([line[0], line[1], line[2], convertStringToNan(line[3]), convertStringToNan(line[4]), convertStringToNan(line[5]), convertStringToNan(line[6]), line[7], convertStringToNan(line[8]), convertStringToNan(line[9]), convertStringToNan(line[10]), convertStringToNan(line[11]), line[12], line[13], line[14], line[15]]) """

""" print('\n### Print WDS table ###\n')
print(wdsTable) """

print('\n### Creating filelist ###')

workingDirectory = sys.argv[2]
directoryContent = os.listdir(workingDirectory)
print('Working directory: ', workingDirectory)

files = [f for f in directoryContent if os.path.isfile(workingDirectory+'/'+f) and f.endswith('.new')]
print('Files:', files)

wdscatalog = SkyCoord(ra=wdsTable['ra_deg']*u.degree, dec=wdsTable['dec_deg']*u.degree)
print('WDS Catalog:\n', wdscatalog)

### Run source detection, collect star data to Qtable
print('\n### Running source detection ###')
for fitsFile in files:
    # 1. Read the list of sources extracted from an image (fits) file
    print('\n\n### Processing file: ', fitsFile, '###')
    fitsFileName = workingDirectory + '/' + fitsFile
    hdu = fits.open(fitsFileName)
    mywcs = WCS(hdu[0].header)

    # Estimate the background and background noise
    data = hdu[0].data
    mean, median, std = sigma_clipped_stats(data, sigma=5.0)  

    daofind = DAOStarFinder(fwhm=10.0, threshold=18.0*std)  
    sources = daofind(data - median)

    for star in sources:
        ra2, dec2 = mywcs.all_pix2world([[star ['xcentroid'], star ['ycentroid']]], 0)[0]   
        mainstar = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
        #catalog = SkyCoord(ra=gaiaStars[1:, 5]*u.degree, dec=gaiaStars[1:, 7]*u.degree)  
        idx, d2d, d3d = mainstar.match_to_catalog_sky(wdscatalog)
        catalogstar = SkyCoord(ra=wdsTable[idx]['ra_deg']*u.degree, dec=wdsTable[idx]['dec_deg']*u.degree)
        sep = mainstar.separation(catalogstar)
        #print('Catalogstar:', catalogstar)
        #print('Separation:', sep)
        if sep < Angle(0.001 * u.deg):
            companion = mainstar.directional_offset_by(wdsTable[idx][3] * u.degree, (wdsTable[idx][4] / 3600) * u.deg)
            #print(companion)
            separation = mainstar.separation(companion)
            #print('Separation check:', separation)
            for star2 in sources:
                ra3, dec3 = mywcs.all_pix2world([[star2['xcentroid'], star2['ycentroid']]], 0)[0]   
                compstar = SkyCoord(ra=ra3*u.degree, dec=dec3*u.degree)  
                #catalog = SkyCoord(ra=gaiaStars[1:, 5]*u.degree, dec=gaiaStars[1:, 7]*u.degree)  
                sep2 = compstar.separation(companion)
                if sep2 <= Angle(0.002 * u.deg) and star != star2:
                    print('\nWDS Identifier:', wdsTable[idx][0], wdsTable[idx][1], wdsTable[idx][2])
                    print('Main star:', wdsTable[idx][0], wdsTable[idx][1], wdsTable[idx][3], wdsTable[idx][4])
                    print('Main star data:\n', star)
                    print('Separation:', Angle(sep))
                    print('Comapnion star data:\n', star2)
                    paactual = mainstar.position_angle(compstar).to(u.deg)
                    sepactual = mainstar.separation(compstar) * 3600
                    print('\nPA:', paactual.degree, 'Sep:', (sepactual.degree))
                    objectid = wdsTable[idx][0]
                    wdsmag_diff = wdsTable[idx]['mag_pri'] - wdsTable[idx]['mag_sec']
                    magdiff = star['mag'] - star2['mag']
                    dsTable.add_row([wdsTable[idx][0], wdsTable[idx][1], wdsTable[idx][2], wdsTable[idx][3], wdsTable[idx][4], wdsTable[idx][5], wdsTable[idx][6], wdsmag_diff, wdsTable[idx][7], wdsTable[idx][8], wdsTable[idx][9], wdsTable[idx][10], wdsTable[idx][11], str(wdsTable[idx][12]), str(wdsTable[idx][13]), str(wdsTable[idx][14]), str(wdsTable[idx][15]), objectid, paactual.degree, sepactual.degree, magdiff])

print('\n### Double stars ###')
print(dsTable)            
        
        