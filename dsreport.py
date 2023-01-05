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

# Function to calculate the separation of the two stars in parsecs
# Excel formula =IF('min distance A'>'min distance b','min distance A'*'Rho','min distance b'*'Rho')
auToParsec = 0.0000048481368111358
def sepCalc(dist_a, dist_b, rho):
    if dist_a > dist_b:
        sep = (dist_a * rho) * auToParsec
    else:
        sep = (dist_b * rho) * auToParsec
    return sep

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
# Excel formula =2.52^(4.83-'Absolute magnitude')
def calcLuminosity(absmag):
    lum = 2.52 ** (4.83 - absmag)
    return lum

# Function to calculate the Star mass
# Excel formula M <0.43M =('luminosity'/0.23)^(1/2.3), M <2M ='luminosity'^(1/4), M < 20M =('luminosity'/1.4)^(1/3.5), M > 55M ='luminosity'/3200
def calcMass(lum):
    mass = ()
    mass_small = (lum / 0.23) ** (1 / 2.3)
    mass_med = lum ** (1 / 4)
    mass_lar = (lum / 1.4) ** (1 / 3.5)
    mass_ex = lum / 3200
    if mass_small <= 0.43:
        mass = mass_small
    elif 0.43 < mass_med < 2:
        mass = mass_med
    elif 2 < mass_lar < 55:
        mass = mass_lar
    elif mass_ex > 55:
        mass = mass_ex
    return mass

# Function to calculate Harshaw probapility of duplicity based on the parallax and proper motion factors
def calcHarshaw(parallaxFactor, pmFactor):
    HarshawFactor = (parallaxFactor * 0.75) + (pmFactor * 0.15)
    return HarshawFactor

# Function to calculate Harshaw physicality of the system based on the parallax and pm factors
def calcHarshawPhysicality(harfac):
    if harfac >= 0.85:
        HarshawPhysicality = 'yes'
    elif 0.65 <= harfac < 0.85:
        HarshawPhysicality = '?'
    elif 0.5 <= harfac < 0.65:
        HarshawPhysicality = 'Maybe'
    elif 0.35 <= harfac < 0.5:
        HarshawPhysicality = '??'
    elif 0.0 <= harfac < 0.35:
        HarshawPhysicality = 'No'
    elif harfac <= 0:
        HarshawPhysicality = 'Something went wrong...'
    return HarshawPhysicality

# Function to calculate the Tangential speed components from proper motin in km/s
# Excel formula to calculate the Proper motion in km/s =pm_ra(dec)/1000*distance from earth*4.74
def calcTangentialSpeedComponent(dist, pm):
    tanspeed = pm/1000*dist*4.74
    return tanspeed

# Function to calculate the Relative velocity
# Excel formula to calculate the difference to the Tangential speeds in km/s =SQRT((pm_ra_a-pm_ra_b)^2+(pm_dec_a-pm_dec_b)^2)
# Excel formula to calculate the difference to the Radial speeds in km/s =ABS(rad_vel_a-rad_vel_b)
# Excel formula to calculate the relative velocity =SQRT(L8^2+L7^2)
def calcRelativeVelocity(pmraa, pmdeca, pmrab, pmdecb, radvela, radvelb, dista, distb):
    tanraa = calcTangentialSpeedComponent(dista, pmraa)
    tandeca = calcTangentialSpeedComponent(dista, pmdeca)
    tanrab = calcTangentialSpeedComponent(distb, pmrab)
    tandecb = calcTangentialSpeedComponent(distb, pmdecb)
    tanspeeddiff = math.sqrt((tanraa - tanrab) ** 2 + (tandeca - tandecb) ** 2)
    radspeeddif = math.fabs(radvela - radvelb)
    sumspeeddiff = math.sqrt(tanspeeddiff ** 2 + radspeeddif ** 2)
    return sumspeeddiff


# Function to calculate the Escape velocity of the system, separation should be calculated in parsec!
gravConst = 0.0043009 # Gravitational constant is convenient if measure distances in parsecs (pc), velocities in kilometres per second (km/s) and masses in solar units M
def calcEscapevelocity(mass_a, mass_b, separation, gravconst):
    escvel = math.sqrt((2 * gravconst * (mass_a + mass_b)) / separation)
    return escvel

# Function to calculate the Probability of binarity based on the Relative and the escape velocity
# Excel formula =IF(relative speed<=escape velocity,"Y","No"))
def calcBinarity(relsped, escsped):
    binarity = ()
    if relsped == 'nan' or escsped == 'nan':
        binarity = 'missing data'
    else:
        if relsped < escsped:
            binarity = 'yes'
        else:
            binarity = 'no'
    return binarity

# Function to calculate the Standard error in RA/DEC measurements
def calcStandardError(arr):
    stderr = np.std(arr)
    return stderr

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
sourceTable = QTable([fileName, sourceId, dr3Designation, dr3Ra, dr3Dec, dr3Parallax, dr3ParallaxError, dr3PmRa, dr3PmDec, dr3gMag, dr3bpMag, dr3rpMag, dr3RadVel, dr3RadVelErr, dr3Temp, imageId, sourceRa, sourceDec, sourceMag], names=('filename', 'source_id', 'designation', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec','phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'radial_velocity', 'radial_velocity_error', 'teff_gspphot', 'image_id', 'source_ra', 'source_dec', 'source_mag'), meta={'name': 'source table'})

# Define Qtable to record Gaia data for each image about the doubles
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
reportRaMeasuredA = np.array([], dtype=np.float64)
reportDecMeasuredA = np.array([], dtype=np.float64)
reportMagMeasuredA = np.array([], dtype=np.float64)
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
reportRaMeasuredB = np.array([], dtype=np.float64)
reportDecMeasuredB = np.array([], dtype=np.float64)
reportMagMeasuredB = np.array([], dtype=np.float64)
reportThetaDr3 = np.array([], dtype=np.float64)
reportThetaMeasured = np.array([], dtype=np.float64)
reportRhoDr3= np.array([], dtype=np.float64)
reportRhoMeasured = np.array([], dtype=np.float64)
reportObjectId = np.array([], dtype=str)
reportTable = QTable([reportFileName, reportSourceIdA, reportDr3DesignationA, reportDr3RaA, reportDr3DecA, reportDr3ParallaxA, reportDr3ParallaxErrorA, reportDr3PmRaA, reportDr3PmDecA, reportDr3gMagA, reportDr3bpMagA, reportDr3rpMagA, reportDr3RadVelA, reportDr3RadVelErrA, reportDr3TempA, reportRaMeasuredA, reportDecMeasuredA, reportMagMeasuredA, reportSourceIdB, reportDr3DesignationB, reportDr3RaB, reportDr3DecB, reportDr3ParallaxB, reportDr3ParallaxErrorB, reportDr3PmRaB, reportDr3PmDecB, reportDr3gMagB, reportDr3bpMagB, reportDr3rpMagB, reportDr3RadVelB, reportDr3RadVelErrB, reportDr3TempB, reportRaMeasuredB, reportDecMeasuredB, reportMagMeasuredB, reportThetaDr3, reportThetaMeasured, reportRhoDr3, reportRhoMeasured, reportObjectId], names=('filename', 'source_id_a', 'designation_a', 'ra_a', 'dec_a', 'parallax_a', 'parallax_error_a', 'pmra_a', 'pmdec_a', 'phot_g_mean_mag_a', 'phot_bp_mean_mag_a', 'phot_rp_mean_mag_a', 'radial_velocity_a', 'radial_velocity_error_a', 'teff_gspphot_a', 'rameasured_a', 'decmeasured_a', 'magmeasured_a', 'source_id_b', 'designation_b', 'ra_b', 'dec_b', 'parallax_b', 'parallax_error_b', 'pmra_b', 'pmdec_b', 'phot_g_mean_mag_b', 'phot_bp_mean_mag_b', 'phot_rp_mean_mag_b', 'radial_velocity_b', 'radial_velocity_error_b', 'teff_gspphot_b', 'rameasured_b', 'decmeasured_b', 'magmeasured_b', 'theta_dr3', "theta_measured", 'rho_dr3', 'rho_measured', 'object_id'), meta={'name': 'report table'})
# , reportMassA, reportMassB, reportAbsMagA, reportAbsMagB, reportLumA, reportLumB, reportEscapeVelocity, reportRelativeVelocity, reportHarshawPhysicality, reportBinarity
# , 'mass_a', 'mass_b', 'absmag_a', 'absmag_b', 'lum_a', 'lum_b', 'sys_esc_vel', 'sys_rel_vel', 'harshaw_physicality', 'binarity'

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
            sourceTable.add_row([fitsFile, gaiaStars[idx + 1][2], 'Gaia DR3 ' + str(int(gaiaStars[idx + 1][2])), gaiaStars[idx + 1][5], gaiaStars[idx + 1][7], gaiaStars[idx + 1][9], gaiaStars[idx + 1][10], gaiaStars[idx + 1][13], gaiaStars[idx + 1][15], gaiaStars[idx + 1][69], gaiaStars[idx + 1][74], gaiaStars[idx + 1][79], gaiaStars[idx + 1][89], gaiaStars[idx + 1][90], gaiaStars[idx + 1][130], star['id'], ra2, dec2, star['mag']])

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
        StarA = (star['filename'], star['source_id'], star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'],star['phot_g_mean_mag'], star['phot_bp_mean_mag'], star['phot_rp_mean_mag'], star['radial_velocity'], star['radial_velocity_error'], star['teff_gspphot'], star['image_id'], star['source_ra'], star['source_dec'], star['source_mag'])
        for star in group:
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
                    starId1 = StarA[1]
                    starName1 = StarA[2]
                    starId2 = StarB[1]
                    starName2 = StarB[2]
                    starParallax2 = float(StarB[5])
                    starParallaxError2 = float(StarB[6])
                    starPmRa1 = float(StarA[7])
                    starPmDec1 = float(StarA[8])
                    starPmRa2 = float(StarB[7])
                    starPmDec2 = float(StarB[8])
                    starGMag1 = float(StarA[9])
                    starGMag2 = float(StarB[9])
                    starBpMag1 = float(StarA[10])
                    starBpMag2 = float(StarB[10])
                    starRpMag1 = float(StarA[11])
                    starRpMag2 = float(StarB[11])
                    starRadVel1 = float(StarA[12])
                    starRadVelErr1 = float(StarA[13])
                    starRadVel2 = float(StarB[12])
                    starRadVelErr2 = float(StarB[13])
                    starTemp1 = float(StarA[14])
                    starTemp2 = float(StarB[14])
                    starImageIdA = float(StarB[15])
                    starImageIdB = float(StarB[15])
                    starActualRa1 = float(StarA[16])
                    starActualDec1 = float(StarA[17])
                    starActualMag1 = float(StarA[18])
                    starActualRa2 = float(StarB[16])
                    starActualDec2 = float(StarB[17])
                    starActualMag2 = float(StarA[18])
                    starObjectId = starId1 + '_' + starId2

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
                    #starParallaxFactor = calcParallaxFactor(starParallax1, starParallax2)
                    #starPmFactor = calcPmFactor(starPmRa1, starPmDec1, starPmRa2, starPmDec2)
                    #starAbsMag1 = calcAbsMag(starGMag1, starParallax1) # Calculate Absolute magnitude
                    #starAbsMag2 = calcAbsMag(starGMag2, starParallax2) # Calculate Absolute magnitude
                    #starLum1 = calcLuminosity(starAbsMag1)
                    #starLum2 = calcLuminosity(starAbsMag2)
                    #starMass1 = calcMass(starLum1)
                    #starMass2 = calcMass(starLum2)
                    #starSepPar = sepCalc(starDistanceMin1, starDistanceMin2, rhoStar) # Separation of the stars in parsecs
                    #starEscapeVelocity = calcEscapevelocity(starMass1, starMass2, starSepPar, gravConst)
                    #starRelativeVelocity = calcRelativeVelocity(starPmRa1, starPmDec1, starPmRa2, starPmDec2, starRadVel1, starRadVel2, starDistanceMin1, #starDistanceMin2)
                    #starHarshawFactor = calcHarshaw(starParallaxFactor, starPmFactor)
                    #starHarshawPhysicality = calcHarshawPhysicality(starHarshawFactor)
                    #starBinarity = calcBinarity(starRelativeVelocity, starEscapeVelocity)

                    # Check if stars shares a common distance range
                    distanceCommon = ()
                    if starDistanceMin1 < starDistanceMin2 < starDistanceMax1 or starDistanceMin2 < starDistanceMin1 < starDistanceMax2:
                        distanceCommon = 'overlapping'
                    else:
                        distanceCommon = 'no'
                    
                    # Check if the pair is a Common Proper Motion pairs (CPM), Similar Proper Motion (SPM) or Different Proper Motion (DPM)


                    #Print data, if stars are close and share a common distance range
                    if distanceCommon == 'overlapping':
                        print(star[0], '|', starName1,'|',starName2,'|',thetaStar,'|',rhoStar,'|',starGMag1,'|',starGMag2,'|',starDistance1,'|',starDistanceMax1,'|',starDistanceMin1,'|',starDistanceRange1,'|',starDistance2,'|',starDistanceMax2,'|',starDistanceMin2,'|',starDistanceRange2,'|',distanceCommon,'|', thetaActual,'|',rhoActual)
                        reportTable.add_row([star[0], starId1, starName1, starRa1, starDec1, starParallax1, starParallaxError1, starPmRa1, starPmDec1, starGMag1, starBpMag1, starRpMag1, starRadVel1, starRadVelErr1, starTemp1, starActualRa1, starActualDec1, starActualMag1, starId2, starName2, starRa2, starDec2, starParallax2, starParallaxError2, starPmRa2, starPmDec2, starGMag2, starBpMag2, starRpMag2, starRadVel2, starRadVelErr2, starTemp2, starActualRa2, starActualDec2, starActualMag2, thetaStar, thetaActual, rhoStar, rhoActual, starObjectId])

print(reportTable)
tableFileName = (str('testQTab.csv'))
reportTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',')

# Create Qtable to list the measurements and the standard error per groups (double star)
measuredObject = np.array([], dtype=str)
measuredStarA = np.array([], dtype=str)
measuredStarB = np.array([], dtype=str)
measuredArrayTheta = np.array([], dtype=np.float64)
measuredArrayRho = np.array([], dtype=np.float64)
measuredMeanTheta = np.array([], dtype=np.float64)
measuredMeanThetaErr = np.array([], dtype=np.float64)
measuredMeanRho = np.array([], dtype=np.float64)
measuredMeanRhoErr = np.array([], dtype=np.float64)
meanTable = QTable([measuredObject, measuredStarA, measuredStarB, measuredArrayTheta, measuredArrayRho, measuredMeanTheta, measuredMeanThetaErr, measuredMeanRho, measuredMeanRhoErr], names=('object_name', 'star_a', 'star_b', 'measurements_theta', 'measurements_rho','mean_theta', 'mean_theta_error', 'mean_rho', 'mean_rho_error'), meta={'name': 'measured table'})

### Search double stars on the image sequence
reportTable_by_object = reportTable.group_by('object_id')
print(reportTable_by_object.groups.keys)

for key, group in zip(reportTable_by_object.groups.keys, reportTable_by_object.groups):
    objectId = group
    objectStarA = 
    objectStarB = 
    objectArrayTheta = 
    objectArrayRho = 
    objectMeanTheta = reportTable_by_object['theta_measured'].groups.aggregate(np.mean)  
    objectMeanThetaErr = reportTable_by_object['theta_measured'].groups.aggregate(np.std)
    objectMeanRho = reportTable_by_object['rho_measured'].groups.aggregate(np.mean)
    objectMeanRhoErr = reportTable_by_object['rho_measured'].groups.aggregate(np.std)



# Create QTable to collect all data for the final record summary per object



pairParallaxFactor = calcParallaxFactor(starParallax1, starParallax2)
pairPmFactor = calcPmFactor(starPmRa1, starPmDec1, starPmRa2, starPmDec2)
pairAbsMag1 = calcAbsMag(starGMag1, starParallax1) # Calculate Absolute magnitude
pairAbsMag2 = calcAbsMag(starGMag2, starParallax2) # Calculate Absolute magnitude
pairLum1 = calcLuminosity(starAbsMag1)
pairLum2 = calcLuminosity(starAbsMag2)
pairMass1 = calcMass(starLum1)
pairMass2 = calcMass(starLum2)
pairSepPar = sepCalc(starDistanceMin1, starDistanceMin2, rhoStar) # Separation of the stars in parsecs
pairEscapeVelocity = calcEscapevelocity(starMass1, starMass2, starSepPar, gravConst)
pairRelativeVelocity = calcRelativeVelocity(starPmRa1, starPmDec1, starPmRa2, starPmDec2, starRadVel1, starRadVel2, starDistanceMin1, starDistanceMin2)
pairHarshawFactor = calcHarshaw(starParallaxFactor, starPmFactor)
pairHarshawPhysicality = calcHarshawPhysicality(starHarshawFactor)
pairBinarity = calcBinarity(starRelativeVelocity, starEscapeVelocity)
pairMassA = np.array([], dtype=np.float64)
pairMassB = np.array([], dtype=np.float64)
pairAbsMagA = np.array([], dtype=np.float64)
pairAbsMagB = np.array([], dtype=np.float64)
pairLumA = np.array([], dtype=np.float64)
pairLumB = np.array([], dtype=np.float64)
pairEscapeVelocity = np.array([], dtype=np.float64)
pairRelativeVelocity = np.array([], dtype=np.float64)
pairHarshawPhysicality = np.array([], dtype=str)
pairBinarity = np.array([], dtype=str)
pairpmCommon = ()
#, starMass1, starMass2, starAbsMag1, starAbsMag2, starLum1, starLum2, starEscapeVelocity, starRelativeVelocity, starHarshawPhysicality, starBinarity
if starPmFactor >= 0.8:
    pmCommon = 'CPM'
elif 0.4 <= starPmFactor < 0.8:
    pmCommon = 'SPM'
elif starPmFactor < 0.4:
    pmCommon = 'DPM'