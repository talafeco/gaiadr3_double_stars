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

# Function to calculate the Star's absolute magnitude
# Excel formula =phot_g_mean_mag-5*LOG10('distance from earth')+5
def calcAbsMag(gmag, par):
    dist = calcDistance(par)
    absmag = gmag - 5 * math.log(dist, 10) + 5
    return absmag

# Function to calculate the Star's luminosity

# Function to calculate the Star mass

# Function to calculate Harshaw probapility of duplicity based on the parallax and proper motion factors
def calcHarshaw(parallaxFactor, pmFactor):
    HarshawFactor = (parallaxFactor * 0.75) + (pmFactor * 0.15)
    return HarshawFactor

# Function to calculate the Relative velocity

# Function to calculate the Escape velocity of the system

# Function to calculate the Probability of binarity based on the Relative and the escape velocity

# Function to calculate the Standard deviation in RA/DEC measurements

# Function to calculate the Standard error in RA/DEC measurements


### Run source detection, collect star data to Qtable
workingDirectory = sys.argv[1]
directoryContent = os.listdir(workingDirectory)
print('Working directory: ', workingDirectory)

files = [f for f in directoryContent if os.path.isfile(workingDirectory+'/'+f) and f.endswith('.new')]
print(files)

# Define Qtable for sources
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
# Create source table
sourceTable = QTable([fileName, sourceId, dr3Designation, dr3Ra, dr3Dec, dr3Parallax, dr3ParallaxError, dr3PmRa, dr3PmDec, dr3gMag, dr3bpMag, dr3rpMag, dr3RadVel, dr3RadVelErr, dr3Temp, imageId, sourceRa, sourceDec, sourceMag], names=('filename', 'source_id', 'designation', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec','phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'radial_velocity', 'radial_velocity_error', 'teff_gspphot', 'image_id', 'source_ra', 'source_dec', 'source_mag'), meta={'name': 'first table'})

# Define Qtable for results
reportFileName = np.array([], dtype=str)
reportSourceIdA = np.array([], dtype=np.int64)
reportDr3DesignationA = np.array([], dtype=str)
reportDr3RaA = np.array([], dtype=np.float64)
reportDr3DecA = np.array([], dtype=np.float64)
reportDr3ParallaxA = np.array([], dtype=np.float64)
reportDr3ParallaxErrorA = np.array([], dtype=np.float64)
reportDr3PmRaA = np.array([], dtype=np.float64)
reportDr3PmDecA = np.array([], dtype=np.float64)
reportDr3gMagA = np.array([], dtype=np.float64)
reportDr3bpMagA = np.array([], dtype=np.float64) # phot_bp_mean_mag
reportDr3rpMagA = np.array([], dtype=np.float64) # phot_rp_mean_mag
reportDr3RadVelA = np.array([], dtype=np.float64) # radial_velocity
reportDr3RadVelErrA = np.array([], dtype=np.float64)  # radial_velocity_error
reportDr3TempA = np.array([], dtype=np.float64)  # teff_gspphot
reportImageIdA = np.array([], dtype=np.int32)
reportSourceIdB = np.array([], dtype=np.int64)
reportDr3DesignationB = np.array([], dtype=str)
reportDr3RaB = np.array([], dtype=np.float64)
reportDr3DecB = np.array([], dtype=np.float64)
reportDr3ParallaxB = np.array([], dtype=np.float64)
reportDr3ParallaxErrorB = np.array([], dtype=np.float64)
reportDr3PmRaB = np.array([], dtype=np.float64)
reportDr3PmDecB = np.array([], dtype=np.float64)
reportDr3gMagB = np.array([], dtype=np.float64)
reportDr3bpMagB = np.array([], dtype=np.float64) # phot_bp_mean_mag
reportDr3rpMagB = np.array([], dtype=np.float64) # phot_rp_mean_mag
reportDr3RadVelB = np.array([], dtype=np.float64) # radial_velocity
reportDr3RadVelErrB = np.array([], dtype=np.float64)  # radial_velocity_error
reportDr3TempB = np.array([], dtype=np.float64)  # teff_gspphot
reportImageIdB = np.array([], dtype=np.int32)
reportMeanRa = np.array([], dtype=np.float64)
reportMeanDec = np.array([], dtype=np.float64)
reportMeanMag = np.array([], dtype=np.float64)
reportMassA = np.array([], dtype=np.float64)
reportMassB = np.array([], dtype=np.float64)
reportAbsMagA = np.array([], dtype=np.float64)
reportAbsMagB = np.array([], dtype=np.float64)
reportLumA = np.array([], dtype=np.float64)
reportLumB = np.array([], dtype=np.float64)
reportThetaAll = np.array([], dtype=np.float64)
reportMeanTheta = np.array([], dtype=np.float64)
reportMeanThetaErr = np.array([], dtype=np.float64)
reportRhoAll = np.array([], dtype=np.float64)
reportMeanRho = np.array([], dtype=np.float64)
reportMeanRhoErr = np.array([], dtype=np.float64)
reportEscapeVelocity = np.array([], dtype=np.float64)
reportRelativeVelocity = np.array([], dtype=np.float64)
reportHarshawPhysicality = np.array([], dtype=str)
reportBinarity = np.array([], dtype=str)

# Create report table
reportTable = QTable([reportFileName, reportSourceIdA, reportDr3DesignationA, reportDr3RaA, reportDr3DecA, reportDr3ParallaxA, reportDr3ParallaxErrorA, reportDr3PmRaA, reportDr3PmDecA, reportDr3gMagA, reportDr3bpMagA, reportDr3rpMagA, reportDr3RadVelA, reportDr3RadVelErrA, reportDr3TempA, reportSourceIdB, reportDr3DesignationB, reportDr3RaB, reportDr3DecB, reportDr3ParallaxB, reportDr3ParallaxErrorB, reportDr3PmRaB, reportDr3PmDecB, reportDr3gMagB, reportDr3bpMagB, reportDr3rpMagB, reportDr3RadVelB, reportDr3RadVelErrB, reportDr3TempB, reportMeanRa, reportMeanDec, reportMeanMag, reportMassA, reportMeanRa, reportMeanDec, reportMeanMag, reportMassB, reportAbsMagA, reportAbsMagB, reportLumA, reportLumB, reportMeanTheta, reportMeanThetaErr, reportMeanRho, reportMeanRhoErr, reportEscapeVelocity, reportRelativeVelocity, reportHarshawPhysicality, reportBinarity], names=('filename', 'source_id_a', 'designation_a', 'ra_a', 'dec_a', 'parallax_a', 'parallax_error_a', 'pmra_a', 'pmdec_a', 'phot_g_mean_mag_a', 'phot_bp_mean_mag_a', 'phot_rp_mean_mag_a', 'radial_velocity_a', 'radial_velocity_error_a', 'teff_gspphot_a', 'source_id_b', 'designation_b', 'ra_b', 'dec_b', 'parallax_b', 'parallax_error_b', 'pmra_b', 'pmdec_b', 'phot_g_mean_mag_b', 'phot_bp_mean_mag_b', 'phot_rp_mean_mag_b', 'radial_velocity_b', 'radial_velocity_error_b', 'teff_gspphot_b', 'mean_ra_a', 'mean_dec_a', 'mean_mag_a', 'mass_a', 'mean_ra_b', 'mean_dec_b', 'mean_mag_b', 'mass_b', 'absmag_a', 'absmag_b', 'lum_a', 'lum_b', 'mean_theta', 'mean_theta_err', 'mean_rho', 'mean_rho_err', 'sys_esc_vel', 'sys_rel_vel', 'harshaw_physicality', 'binarity'), meta={'name': 'first table'})

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
        #StarA = (star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'], star['phot_g_mean_mag'], star['source_ra'], star['source_dec'])
        StarA = (star['filename'], star['source_id'], star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'],star['phot_g_mean_mag'], star['phot_bp_mean_mag'], star['phot_rp_mean_mag'], star['radial_velocity'], star['radial_velocity_error'], star['teff_gspphot'], star['image_id'], star['source_ra'], star['source_dec'], star['source_mag'])
        for star in group:
            #StarB = (star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'], star['phot_g_mean_mag'], star['source_ra'], star['source_dec'])
            StarB = (star['filename'], star['source_id'], star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'],star['phot_g_mean_mag'], star['phot_bp_mean_mag'], star['phot_rp_mean_mag'], star['radial_velocity'], star['radial_velocity_error'], star['teff_gspphot'], star['image_id'], star['source_ra'], star['source_dec'], star['source_mag'])
            if StarA != StarB and float(StarA[9]) < float(StarB[9]) and float(StarA[5]) != 0 and float(StarB[5]) != 0:
                #Set input data
                starRa1 = float(StarA[3])
                starDec1 = float(StarA[4])
                starRa2 = float(StarB[3])
                starDec2 = float(StarB[4])
                starParallax1 = float(StarA[5])
                starParallaxError1 = float(StarA[6])
                            
                # Calculate the widest possible separation for StarA
                possSep1 = 10000 / calcDistanceMax(starParallax1, starParallaxError1)
                rhoStar = rhoCalc(starRa1, starDec1, starRa2, starDec2)
                if possSep1 > rhoStar:
                    starName1 = StarA[2]
                    starName2 = StarB[2]
                    starParallax2 = float(StarB[5])
                    starParallaxError2 = float(StarB[6])
                    starPmRa1 = float(StarA[7])
                    starPmDec1 = float(StarA[8])
                    starPmRa2 = float(StarB[7])
                    starPmDec2 = float(StarB[8])
                    starGMag1 = float(StarA[9])
                    starGMag2 = float(StarB[9])
                    starActualRa1 = float(StarA[16])
                    starActualDec1 = float(StarA[17])
                    starActualRa2 = float(StarB[16])
                    starActualDec2 = float(StarB[17])

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

                    # Calculate Harshaw double star factor
                    HarshawPhisicality = ()
                    if distanceCommon == 'overlapping':
                        HarshawFactor = calcHarshaw(starParallaxFactor, starPmFactor)
                        if HarshawFactor > 0.85:
                            HarshawPhisicality = 'yes'
                        elif 0.65 < HarshawFactor < 0.85:
                            HarshawPhisicality = '?'
                        elif 0.5 < HarshawFactor < 0.65:
                            HarshawPhisicality = 'Maybe'
                        elif 0.35 < HarshawFactor < 0.5:
                            HarshawPhisicality = '??'
                        elif 0.0 < HarshawFactor < 0.35:
                            HarshawPhysicality = 'No'
                    
                    # Calculate Absolute magnitudes of the stars
                    absMagA = calcAbsMag(starGMag1, starParallax1)
                    absMagB = calcAbsMag(starGMag2, starParallax2)
                    
                    #Print data, if stars are close and share a common distance range
                    if distanceCommon == 'overlapping':
                        print(star[0], '|', starName1,'|',starName2,'|',thetaStar,'|',rhoStar,'|',starGMag1,'|',starGMag2,'|',starDistance1,'|',starDistanceMax1,'|',starDistanceMin1,'|',starDistanceRange1,'|',starDistance2,'|',starDistanceMax2,'|',starDistanceMin2,'|',starDistanceRange2,'|',distanceCommon,'|',starParallaxFactor,'|',starPmFactor,'|',pmCommon, '|', HarshawPhysicality, '|',thetaActual,'|',rhoActual)