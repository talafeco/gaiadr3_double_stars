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
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import warnings
warnings.filterwarnings("ignore")

### Declare functions
# Function to calculate Delta RA
def deltaRa(raa, rab, decb):
    deltara = (rab - raa) * math.cos(math.radians(decb))
    return deltara

# Function to calculate Delta DEC
def deltaDec(decb, deca):
    deltadec = decb - deca
    return deltadec

# Function to calculate Theta (position angle)
def thetaCalc(deltara,deltadec):
    thetacalc = math.degrees(math.atan(deltara/deltadec))
    return thetacalc

# Function to calculate Rho (separation)
def rhoCalc(raa, deca, rab, decb):
    rhocalc = math.sqrt(((raa-rab) * math.cos(math.radians(deca))) ** 2 + (deca - decb) ** 2) * 3600
    return rhocalc

# Function to calculate the distance of the star based on the parallax
def calcDistance(par):
    dist = 1 / (par/1000)
    return dist

# Function to calculate the minimum distance of the star based on the parallax error
def calcDistanceMin(par, err):
    distmin = 1 / ((par + err) / 1000)
    return distmin

# Function to calculate the maximum distance of the star based on the parallax error
def calcDistanceMax(par, err):
    distmax = 1 / ((par - err) / 1000)
    return distmax

# Function to calculate the parallax factor of the pair
# EXCEL formula =1-ABS((raA-raB)/(0,5*(raA+raB)))
def calcParallaxFactor(para, parb):
    parfac = 1 - math.fabs((para-parb)/(0.5*(para+parb)))
    return parfac

# Function to calculate the proper motion factor of the stars and define if it is a CPM (common proper motion) pair
# EXCEL formula =ABS(1-((SQRT(pmraA-pmraB)^2+(pmdecA-pmdecB)^2)/(SQRT(pmraA^2+pmdecA^2)+(pmraB^2+pmdecB^2)))))
def calcPmFactor(pmraa, pmdeca, pmrab, pmdecb):
    pmfac = math.fabs(1-((math.sqrt(((pmraa-pmrab) ** 2) + ((pmdeca-pmdecb) ** 2))/(math.sqrt((pmraa ** 2) + (pmdeca ** 2))+((pmrab ** 2) + pmdecb ** 2)))))
    return pmfac

### Run source detection, collect star data to Qtable
workingDirectory = sys.argv[1]
directoryContent = os.listdir(workingDirectory)
print('Working directory: ', workingDirectory)

files = [f for f in directoryContent if os.path.isfile(workingDirectory+'/'+f) and f.endswith('.new')]
print(files)

# Define Qtable for results
fileName = np.array([], dtype=str)
sourceId = np.array([], dtype=np.int64)
dr3Designation = np.array([], dtype=str)
dr3Ra = np.array([], dtype=np.float64)
dr3Dec = np.array([], dtype=np.float64)
dr3Parallax = np.array([], dtype=np.float64)
dr3ParallaxError = np.array([], dtype=np.float64)
dr3PmRa = np.array([], dtype=np.float64)
dr3PmDec = np.array([], dtype=np.float64)
dr3gMag = np.array([], dtype=np.float64)
dr3bpMag = np.array([], dtype=np.float64) # phot_bp_mean_mag
dr3rpMag = np.array([], dtype=np.float64) # phot_rp_mean_mag
dr3RadVel = np.array([], dtype=np.float64) # radial_velocity
dr3RadVelErr = np.array([], dtype=np.float64)  # radial_velocity_error
dr3Temp = np.array([], dtype=np.float64)  # teff_gspphot
imageId = np.array([], dtype=np.int32)
sourceRa = np.array([], dtype=np.float64)
sourceDec = np.array([], dtype=np.float64)
sourceMag = np.array([], dtype=np.float64)
sourceTable = QTable([fileName, sourceId, dr3Designation, dr3Ra, dr3Dec, dr3Parallax, dr3ParallaxError, dr3PmRa, dr3PmDec, dr3gMag, dr3bpMag, dr3rpMag, dr3RadVel, dr3RadVelErr, dr3Temp, imageId, sourceRa, sourceDec, sourceMag], names=('filename', 'source_id', 'designation', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec','phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'radial_velocity', 'radial_velocity_error', 'teff_gspphot', 'image_id', 'source_ra', 'source_dec', 'source_mag'), meta={'name': 'first table'})

for fitsFile in files:
    # 1. Read the list of sources extracted from an image (fits) file
    print('Processing file: ', fitsFile)
    fitsFileName = workingDirectory + '/' + fitsFile
    hdu = fits.open(fitsFileName)
    mywcs = WCS(hdu[0].header)

    # Estimate the background and background noise
    data = hdu[0].data
    mean, median, std = sigma_clipped_stats(data, sigma=5.0)  

    daofind = DAOStarFinder(fwhm=10.0, threshold=18.0*std)  
    sources = daofind(data - median)

    # 2. Define the catalog file based on the source coordinate and read data from catalog file(s) to a catalog
    segments = []
    for star in sources:
        ra, dec = mywcs.all_pix2world([[star ['xcentroid'], star ['ycentroid']]], 0)[0]
        # Calculate the sky segment, which will indicate, which file needs to be populated with the star
        # print(star)
        segmentRaCalc = int((float(ra) // 5) + 1)
        segmentDecCalc = int((float(dec) // 5) + 1)
        segmentName = f"{segmentRaCalc}-{segmentDecCalc}.csv"
        if segmentName not in segments:
            segments.append(segmentName)

    # Read all segments into an array
    gaiaStars = np.empty((0, 152), float)

    # Add all segments to the numpy array
    for seg in segments:
        segmentpart = np.genfromtxt(f"/home/gergo/Documents/dr3_catalog/gaiadr3_15mag_catalog/{seg}", delimiter=",", skip_header=1)
        gaiaStars = np.append(gaiaStars, segmentpart, axis=0)

    # Search sources in the segment catalog
    for star in sources:
        ra2, dec2 = mywcs.all_pix2world([[star ['xcentroid'], star ['ycentroid']]], 0)[0]   
        c = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
        catalog = SkyCoord(ra=gaiaStars[1:, 5]*u.degree, dec=gaiaStars[1:, 7]*u.degree)  
        idx, d2d, d3d = c.match_to_catalog_sky(catalog)
        catalogstar = SkyCoord(ra=gaiaStars[idx + 1][5]*u.degree, dec=gaiaStars[idx + 1][7]*u.degree)
        sep = c.separation(catalogstar)
        if sep < Angle('00d00m02s'):
            sourceTable.add_row([fitsFile, gaiaStars[idx + 1][2], 'Gaia DR3 ' + str(int(gaiaStars[idx + 1][2])), gaiaStars[idx + 1][5], gaiaStars[idx + 1][7], gaiaStars[idx + 1][9], gaiaStars[idx + 1][10], gaiaStars[idx + 1][12], gaiaStars[idx + 1][14], gaiaStars[idx + 1][69], gaiaStars[idx + 1][74], gaiaStars[idx + 1][79], gaiaStars[idx + 1][89], gaiaStars[idx + 1][90], gaiaStars[idx + 1][130], star['id'], ra2, dec2, star['mag']])

# Write found sources into file
tableFileName = (workingDirectory + '/' + str(fitsFile[:-4] + '.csv'))
sourceTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',')

### Search double stars on the image sequence
sourceTable_by_file = sourceTable.group_by('filename')
print(sourceTable_by_file.groups.keys)

for key, group in zip(sourceTable_by_file.groups.keys, sourceTable_by_file.groups):
    # Creating empty arrays for Star related calculations
    StarA = []
    StarB = []
    for star in group: ## modify according to arrays instead of starlist
        StarA = (star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'], star['phot_g_mean_mag'], star['source_ra'], star['source_dec'])
        for star in group:
            StarB = (star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'], star['phot_g_mean_mag'], star['source_ra'], star['source_dec'])
            if StarA != StarB and float(StarA[7]) < float(StarB[7]) and float(StarA[3]) != 0 and float(StarB[3]) != 0:
                #Set input data
                starRa1 = float(StarA[1])
                starDec1 = float(StarA[2])
                starRa2 = float(StarB[1])
                starDec2 = float(StarB[2])
                starParallax1 = float(StarA[3])
                starParallaxError1 = float(StarA[4])
                            
                # Calculate the widest possible separation for StarA
                possSep1 = 10000 / calcDistanceMax(starParallax1, starParallaxError1)
                rhoStar = rhoCalc(starRa1, starDec1, starRa2, starDec2)
                if possSep1 > rhoStar:
                    starName1 = StarA[0]
                    starName2 = StarB[0]
                    starParallax2 = float(StarB[3])
                    starParallaxError2 = float(StarB[4])
                    starPmRa1 = float(StarA[5])
                    starPmDec1 = float(StarA[6])
                    starPmRa2 = float(StarB[5])
                    starPmDec2 = float(StarB[6])
                    starGMag1 = float(StarA[7])
                    starGMag2 = float(StarB[7])
                    starActualRa1 = float(StarA[8])
                    starActualDec1 = float(StarA[9])
                    starActualRa2 = float(StarB[8])
                    starActualDec2 = float(StarB[9])

                    # Value to modify Theta according to the appropriate quadrant
                    addThetaValue = ()
                    if deltaRa(starRa1, starRa2, starDec2) > 0 and deltaDec(starDec2, starDec1) > 0:
                        addThetaValue = 0
                    elif deltaRa(starRa1, starRa2, starDec2) > 0 and deltaDec(starDec2, starDec1) < 0:
                        addThetaValue = 180
                    elif deltaRa(starRa1, starRa2, starDec2) < 0 and deltaDec(starDec2, starDec1) < 0:
                        addThetaValue = 180
                    elif deltaRa(starRa1, starRa2, starDec2) < 0 and deltaDec(starDec2, starDec1) > 0:
                        addThetaValue = 360
                    
                    # Calculate actual data based on functions
                    thetaStar = thetaCalc(deltaRa(starRa1, starRa2, starDec2), deltaDec(starDec2, starDec1)) + addThetaValue
                    thetaActual = thetaCalc(deltaRa(starActualRa1, starActualRa2, starActualDec2), deltaDec(starActualDec2, starActualDec1)) + addThetaValue
                    rhoActual = rhoCalc(starActualRa1, starActualDec1, starActualRa2, starActualDec2)
                    starDistance1 = calcDistance(starParallax1)
                    starDistanceMax1 = calcDistanceMax(starParallax1, starParallaxError1)
                    starDistanceMin1 = calcDistanceMin(starParallax1, starParallaxError1)
                    starDistanceRange1 = starDistanceMax1 - starDistanceMin1
                    starDistance2 = calcDistance(starParallax2)
                    starDistanceMax2 = calcDistanceMax(starParallax2, starParallaxError2)
                    starDistanceMin2 = calcDistanceMin(starParallax2, starParallaxError2)
                    starDistanceRange2 = starDistanceMax2 - starDistanceMin2
                    starParallaxFactor = calcParallaxFactor(starParallax1, starParallax2)
                    starPmFactor = calcPmFactor(starPmRa1, starPmDec1, starPmRa2, starPmDec2)
                    
                    # Check if stars shares a common distance range
                    distanceCommon = ()
                    if starDistanceMin1 < starDistanceMin2 < starDistanceMax1 or starDistanceMin2 < starDistanceMin1 < starDistanceMax2:
                        distanceCommon = 'overlapping'
                    else:
                        distanceCommon = 'no'
                    
                    # Check if the pair is a Common Proper Motion pairs (CPM), Similar Proper Motion (SPM) or Different Proper Motion (DPM)
                    pmCommon = ()
                    if starPmFactor >= 0.8:
                        pmCommon = 'CPM'
                    elif 0.4 <= starPmFactor < 0.8:
                        pmCommon = 'SPM'
                    elif starPmFactor < 0.4:
                        pmCommon = 'DPM'
                    
                    #Print data, if stars are close and share a common distance range
                    if distanceCommon == 'overlapping':
                        print(star[0], '|', starName1,'|',starName2,'|',thetaStar,'|',rhoStar,'|',starGMag1,'|',starGMag2,'|',starDistance1,'|',starDistanceMax1,'|',starDistanceMin1,'|',starDistanceRange1,'|',starDistance2,'|',starDistanceMax2,'|',starDistanceMin2,'|',starDistanceRange2,'|',distanceCommon,'|',starParallaxFactor,'|',starPmFactor,'|',pmCommon,'|',thetaActual,'|',rhoActual)