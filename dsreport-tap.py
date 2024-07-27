#! /usr/bin/python3

# DSReport tool to find phisical double stars on astronomical images based on Gaia DR3 data
# Version: 1.0
# Install: copy file to /usr/local/bin folder
# Usage: dsreport <image_folder>

#tasks
#rhocalc to be fixed!
#show progress by x/y files processed - OK
#show the number of measurements - OK

import os
import sys
import numpy as np
import math
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, FK5, ICRS, match_coordinates_sky, search_around_sky, position_angle, angular_separation
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta
from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
Gaia.ROW_LIMIT = -1 # To return an unlimited number of rows

import warnings
warnings.filterwarnings("ignore")

# List of constances
#wdsFile = Table.read(f"/usr/share/dr3map/wds/dr3-wds.csv", format='ascii')
hipparcos_file = Table.read(f"/usr/share/dr3map/hipparcos/I_239_selection.csv", format='ascii')
segment_lib = "/usr/share/dr3map/gaiadr3_15mag_catalog/"
hipparcos_abs_mag = hipparcos_file['Abs_mag']
hipparcos_bv_index = hipparcos_file['B-V']

# Configuration for the ATIK camera

dao_sigma = 3.0
dao_fwhm = 9.0
dao_threshold = 9.0
possible_distance = 30000.0 # AU
search_cone = 0.002 # Decimal degree
image_limit = 2000

# Gravitational constant is convenient if measure distances in parsecs (pc), velocities in kilometres per second (km/s) and masses in solar units M
gravConst = 0.0043009 

# Constant to calculate star luminosity and mass
sun_luminosity = 3.0128 * (10 ** 28)
sun_absolute_luminosity = 3.828 * (10 ** 26)

# Configuration for the CANON camera
'''
dao_sigma = 2.0
dao_fwhm = 3.0
dao_threshold = 5.0
possible_distance = 10000.0 # AU
search_cone = 0.001 # Decimal degree
'''

#########################
### Declare functions ###
#########################

# Function to calculate Star positions based on Gaia DR3 coordinates and proper motion
def calcCurrentDR3Coord(date, star_ra, star_dec, star_pr_ra, star_pr_dec):
    date_time = Time(date, format='jyear')
    star = SkyCoord(ra = star_ra * u.degree,
                dec = star_dec * u.degree,
                pm_ra_cosdec = star_pr_ra * u.mas/u.yr,
                pm_dec = star_pr_dec * u.mas/u.yr,
                frame = 'icrs',
                obstime=Time('2016-01-01 00:00:00.0'))
    starUpdCoord = star.apply_space_motion(new_obstime=date_time)
    return starUpdCoord

# Function to check, if a number is 'nan'
def isNaN(num):
    return num != num

# Function to create precise coordinates for Epoch 2000
def getPreciseCoord(ra, dec, date):
    coord_object = {}
    coords = [str(ra) + ' ' + str(dec)]
    if isNaN(date):
        coord_object = SkyCoord(coords, frame='icrs', unit=(u.degree, u.degree))
    else:
        coord_object = SkyCoord(coords, frame='icrs', unit=(u.degree, u.degree), obstime=date)
    j2000_coord = coord_object.transform_to(FK5(equinox='J2000.0'))
    j2000_coords = j2000_coord.to_string(style='hmsdms', precision=2)[0]
    j2000_coord_formatted = str(j2000_coords).replace("d",":").replace("h",":").replace("m",":").replace("s","")
    return j2000_coord_formatted

# Function to calculate utc from cet
def getUTC(date_time):
    utc_date_time = ''
    if isNaN(date_time):
        utc_date_time = str(date_time)
    else:
        date_of_observation_time_cet = Time(date_time, precision=0)
        time_zone_delta = TimeDelta(-3600, format='sec')
        date_of_observation_time_utc = date_of_observation_time_cet + time_zone_delta
        utc_date_time = str(date_of_observation_time_utc.jyear)
    return utc_date_time


# Function to calculate Delta RA
'''def deltaRa(raa, rab, decb):
    deltara = (rab - raa) * math.cos(math.radians(decb))
    return deltara'''

# Function to calculate Delta DEC
'''def deltaDec(decb, deca):
    deltadec = decb - deca
    return deltadec'''

# Function to calculate Theta (position angle)
'''def thetaCalc(deltara,deltadec):
    thetacalc = math.degrees(math.atan(deltara/deltadec))
    return thetacalc'''

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

def calcLuminosityAlternate(absmag):
    lum = sun_luminosity * (10 ** (-absmag / 2.512))
    lum2 = lum / sun_absolute_luminosity
    return lum2

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
    return float(HarshawFactor)

# Function to calculate Harshaw physicality of the system based on the parallax and pm factors
def calcHarshawPhysicality(harfac):
    if harfac >= 85.0:
        HarshawPhysicality = 'yes'
    elif 65.0 <= harfac < 85.0:
        HarshawPhysicality = '?'
    elif 50.0 <= harfac < 65.0:
        HarshawPhysicality = 'Maybe'
    elif 35.0 <= harfac < 50.0:
        HarshawPhysicality = '??'
    elif 0.0 <= harfac < 35.0:
        HarshawPhysicality = 'No'
    elif harfac <= 0.0:
        HarshawPhysicality = 'Something went wrong... (so NO)'
    else:
        HarshawPhysicality = 'Something went wrong... (so NO)'
    return HarshawPhysicality

# Function to calculate the Tangential speed components from proper motin in km/s
# Excel formula to calculate the Proper motion in km/s =pm_ra(dec)/1000*distance from earth*4.74(au per year to km/s)
def calcTangentialSpeedComponent(dist, pm):
    tanspeed = pm/1000*dist*4.74372
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

# Calculate Common Proper Motion category
def calcPmCategory(pmfact):
    pmCommon = ()
    if pmfact >= 80.0:
        pmCommon = 'CPM'
    elif 40.0 <= pmfact < 80.0:
        pmCommon = 'SPM'
    elif pmfact < 40.0:
        pmCommon = 'DPM'
    return pmCommon

# Convert string to numpy 'nan'
def convertStringToNan(str):
    if str == 'null':
        str = np.nan
    return str

# Search pair in Washington double Star Catalog
'''def searchWds(pairadesig):
    wdsIdx = np.where(pairadesig == wdsFile['designation'])
    wdsRow = wdsFile[wdsIdx]
    if wdsRow:
        wdsPair = wdsRow[0]
    elif not wdsRow:
        wdsPair = 'Not found'
    return wdsPair'''

# Create HRD plot of the double stars based on Hipparcos
# Create HRD plot of the double stars based on Hipparcos
def hrdPlot(designation_a, designation_b, mag_abs_a, mag_abs_b, bv_a, bv_b):
    print(designation_a, designation_b, mag_abs_a, mag_abs_b, bv_a, bv_b)
    if designation_a and designation_b and mag_abs_a and mag_abs_b and bv_a and bv_b:
        hipparcos_abs_mag = hipparcos_file['Abs_mag']
        hipparcos_bv_index = hipparcos_file['B-V']
        colors = (hipparcos_bv_index)
        plt.scatter(hipparcos_bv_index, hipparcos_abs_mag, c=colors, s=0.5, alpha=0.1, cmap='RdYlBu_r', vmax=1.9, vmin=-0.4) #, 
        plt.scatter(bv_a, mag_abs_a, s=14, color="blue", label='Main star') # s= 1 / mag_abs_a
        plt.scatter(bv_b, mag_abs_b, s=7, color="red", label='Companion star') # s= 1 / mag_abs_a
        plt.legend(loc="upper left")
        plt.axis((-0.4,1.9,15,-10))
        plt.title('Double Star ' + designation_a + ' - ' + designation_b + ' H-R Diagram')
        plt.xlabel('B-V index')
        plt.ylabel('Absolute magnitude')
        plt.gca().set_aspect(0.1)
        savename = str(workingDirectory + '/' + designation_a + ' - ' + designation_b + '_hrd.jpg').replace(' ', '')
        plt.savefig(savename, bbox_inches='tight', dpi=300.0)
        plt.close()
    else:
        print('Data is missiong, HRD plot cannot be created!')

'''def hrdPlot(designation_a, designation_b, mag_abs_a, mag_abs_b, bv_a, bv_b):
    hipparcos_abs_mag = hipparcos_file['Abs_mag']
    hipparcos_bv_index = hipparcos_file['B-V']
    plt.scatter(hipparcos_bv_index, hipparcos_abs_mag, s=0.5, alpha=0.2, color="grey") #, 
    plt.scatter(bv_a, mag_abs_a, s= 1 / mag_abs_a * 40, color="blue", label='Main star')
    plt.scatter(bv_b, mag_abs_b, s= 1 / mag_abs_b * 40, color="red", label='Companion star')
    plt.legend(loc="upper left")
    plt.axis((-0.4,1.9,21,-16))
    plt.title('Double Star ' + designation_a + ' - ' + designation_b + ' H-R Diagram')
    plt.xlabel('B-V index')
    plt.ylabel('Absolute magnitude')
    plt.gca().set_aspect(0.07)
    plt.savefig(workingDirectory + '/' + designation_a + ' - ' + designation_b + '_hrd.jpg', dpi=150.0, bbox_inches='tight')
    plt.close()'''

# Create Image plot of the double stars
'''def imagePlot(filename, designation_a, designation_b, raa, deca, rab, decb):
    image_data = fits.open(workingDirectory + '/' + filename)
    wcs_helix = WCS(image_data[0].header)
    image = image_data[0].data
    star_a = SkyCoord(raa * u.deg, deca * u.deg, frame='icrs')
    star_b = SkyCoord(rab * u.deg, decb * u.deg, frame='icrs')
    star_a_pix = utils.skycoord_to_pixel(star_a, wcs_helix)
    star_b_pix = utils.skycoord_to_pixel(star_b, wcs_helix)
    plt.figure(figsize=(10, 10), frameon=False) # 
    ax = plt.subplot(projection=wcs_helix)
    plt.scatter(star_a_pix[0] + 30, star_a_pix[1], marker="_", s=50, color="grey")
    plt.scatter(star_a_pix[0], star_a_pix[1] + 30, marker="|", s=50, color="grey")
    plt.scatter(star_b_pix[0] + 30, star_b_pix[1], marker="_", s=50, color="grey")
    plt.scatter(star_b_pix[0], star_b_pix[1] + 30, marker="|", s=50, color="grey")
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='grey', ls='dotted')
    plt.imshow(image, origin='lower',cmap='grey', aspect='equal', vmax=2000, vmin=0) # , cmap='cividis'
    plt.savefig(workingDirectory + '/' + str(designation_a) + '-' + str(designation_b) + '_img.jpg', bbox_inches='tight')'''

# Create Image plot of the double stars
def imagePlot(filename, designation_a, designation_b, raa, deca, rab, decb):
    '''coord_meta = {}
    coord_meta['type'] = ('longitude', 'latitude')
    coord_meta['wrap'] = (None, None)
    coord_meta['unit'] = (u.degree, u.degree)
    coord_meta['format_unit'] = (u.hour, u.degree)
    coord_meta['name'] = 'ra', 'dec'
    '''
    
    image_data = fits.open(workingDirectory + '/' + filename)
    header = image_data[0].header
    wcs_helix = WCS(image_data[0].header, naxis=2)
    image = image_data[0].data[0]
    image_height = header['NAXIS2']

    star_a = SkyCoord(raa * u.deg, deca * u.deg, frame='icrs')
    star_b = SkyCoord(rab * u.deg, decb * u.deg, frame='icrs')
    star_a_pix = utils.skycoord_to_pixel(star_a, wcs_helix)
    star_b_pix = utils.skycoord_to_pixel(star_b, wcs_helix)
    '''
    ax = plt.subplot(projection=wcs_helix)
    
    plt.scatter(star_a_pix[0] + (image_height / 30), star_a_pix[1], marker="_", s=(image_height / 10), color="grey")
    plt.scatter(star_a_pix[0], star_a_pix[1] + (image_height / 30), marker="|", s=(image_height / 10), color="grey")
    plt.scatter(star_b_pix[0] + (image_height / 30), star_b_pix[1], marker="_", s=(image_height / 10), color="grey")
    plt.scatter(star_b_pix[0], star_b_pix[1] + (image_height / 30), marker="|", s=(image_height / 10), color="grey")
    overlay = ax.get_coords_overlay('icrs', coord_meta=coord_meta)
    overlay.grid(color='grey', ls='dotted')
    overlay['ra'].set_axislabel('Lon')
    overlay['dec'].set_axislabel('Lat')
    overlay['ra'].set_ticklabel_position('bt')
    overlay['ra'].set_ticks(number=6)
    overlay['dec'].set_ticklabel_position('lr')
    overlay['dec'].set_ticks(number=6)'''

    plt.scatter(star_a_pix[0] + 40, star_a_pix[1], marker="_", s=50, color="grey")
    plt.scatter(star_a_pix[0], star_a_pix[1] + 40, marker="|", s=50, color="grey")
    plt.scatter(star_b_pix[0] + 40, star_b_pix[1], marker="_", s=50, color="grey")
    plt.scatter(star_b_pix[0], star_b_pix[1] + 40, marker="|", s=50, color="grey")

    #plt.title(str(designation_a) + ' - ' + str(designation_b), pad=50.0)
    plt.title(str(designation_a) + ' - ' + str(designation_b))

    plt.imshow(image, origin='lower',cmap='Greys_r', aspect='equal', vmax=image_limit, vmin=0) # , cmap='cividis'
    plt.savefig(workingDirectory + '/' + str(designation_a) + '-' + str(designation_b) + '_img.jpg', dpi=150.0, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    #plt.show()

### Run source detection, collect star data to Qtable
workingDirectory = sys.argv[1]
directoryContent = os.listdir(workingDirectory)
print('Working directory: ', workingDirectory)

files = [f for f in directoryContent if os.path.isfile(workingDirectory+'/'+f) and (f.endswith('.new') or f.endswith('.fit') or f.endswith('.fits'))]
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
imageDate = np.array([], dtype=str)
# Create source table
sourceTable = QTable([fileName, sourceId, dr3Designation, dr3Ra, dr3Dec, dr3Parallax, dr3ParallaxError, dr3PmRa, dr3PmDec, dr3gMag, dr3bpMag, dr3rpMag, dr3RadVel, dr3RadVelErr, dr3Temp, imageId, sourceRa, sourceDec, sourceMag, imageDate], names=('filename', 'source_id', 'designation', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec','phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'radial_velocity', 'radial_velocity_error', 'teff_gspphot', 'image_id', 'source_ra', 'source_dec', 'source_mag', 'image_date'), meta={'name': 'source table'})

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
reportImageDate = np.array([], dtype=str)
reportTable = QTable([reportFileName, reportSourceIdA, reportDr3DesignationA, reportDr3RaA, reportDr3DecA, reportDr3ParallaxA, reportDr3ParallaxErrorA, reportDr3PmRaA, reportDr3PmDecA, reportDr3gMagA, reportDr3bpMagA, reportDr3rpMagA, reportDr3RadVelA, reportDr3RadVelErrA, reportDr3TempA, reportRaMeasuredA, reportDecMeasuredA, reportMagMeasuredA, reportSourceIdB, reportDr3DesignationB, reportDr3RaB, reportDr3DecB, reportDr3ParallaxB, reportDr3ParallaxErrorB, reportDr3PmRaB, reportDr3PmDecB, reportDr3gMagB, reportDr3bpMagB, reportDr3rpMagB, reportDr3RadVelB, reportDr3RadVelErrB, reportDr3TempB, reportRaMeasuredB, reportDecMeasuredB, reportMagMeasuredB, reportThetaDr3, reportThetaMeasured, reportRhoDr3, reportRhoMeasured, reportObjectId, reportImageDate], names=('filename', 'source_id_a', 'designation_a', 'ra_a', 'dec_a', 'parallax_a', 'parallax_error_a', 'pmra_a', 'pmdec_a', 'phot_g_mean_mag_a', 'phot_bp_mean_mag_a', 'phot_rp_mean_mag_a', 'radial_velocity_a', 'radial_velocity_error_a', 'teff_gspphot_a', 'rameasured_a', 'decmeasured_a', 'magmeasured_a', 'source_id_b', 'designation_b', 'ra_b', 'dec_b', 'parallax_b', 'parallax_error_b', 'pmra_b', 'pmdec_b', 'phot_g_mean_mag_b', 'phot_bp_mean_mag_b', 'phot_rp_mean_mag_b', 'radial_velocity_b', 'radial_velocity_error_b', 'teff_gspphot_b', 'rameasured_b', 'decmeasured_b', 'magmeasured_b', 'theta_dr3', 'theta_measured', 'rho_dr3', 'rho_measured', 'object_id', 'image_date'), meta={'name': 'report table'})

# Set observation date and time


file_counter = 0

for fitsFile in files:
    # 1. Read the list of sources extracted from an image (fits) file
    file_counter = file_counter + 1
    print('Processing file', file_counter, 'out of', len(files),': ', fitsFile)
    fitsFileName = workingDirectory + '/' + fitsFile
    hdu = fits.open(fitsFileName)
    mywcs = WCS(hdu[0].header, naxis=2)
    file_header = hdu[0].header

    fitsFileDate = ''
    key_to_lookup_a = 'DATE-OBS'
    key_to_lookup_b = 'DATE'
    if key_to_lookup_a in file_header:
        fitsFileDate = file_header['DATE-OBS']
    elif key_to_lookup_b in file_header:
        fitsFileDate = file_header['DATE']
    else:
        fitsFileDate = np.nan

    # Estimate the background and background noise
    data = hdu[0].data
    mean, median, std = sigma_clipped_stats(data, sigma=dao_sigma)  

    daofind = DAOStarFinder(fwhm=dao_fwhm, threshold=dao_threshold*std)  
    sources = daofind(data[0] - median)
    
    photo_left_upper = SkyCoord.from_pixel(0, 0, mywcs, origin=0, mode='all')
    photo_right_lower = SkyCoord.from_pixel(file_header['NAXIS2'], file_header['NAXIS1'], mywcs, origin=0, mode='all')
    photo_center = SkyCoord(file_header['CRVAL1'] * u.degree, file_header['CRVAL2'] * u.degree)
    photo_radius = photo_left_upper.separation(photo_right_lower) / 2
    print('Center of photo: ', photo_center.to_string('hmsdms'), '/', photo_center.to_string('decimal'),
      '\nRadius of photo: ', photo_radius)
    gaia_photo_catalog = Gaia.cone_search_async(photo_center, radius=u.Quantity(photo_radius))
    gaiaStars = gaia_photo_catalog.get_results()
 
    #dr3TableFileName = (str('dr3stars.csv'))
    #gaiaStars.write(dr3TableFileName, format='ascii.ecsv', overwrite=True, delimiter=',')
    # Search sources in the segment catalog
    for star in sources:
        ra2, dec2 = mywcs.all_pix2world([[star ['xcentroid'], star ['ycentroid']]], 0)[0]   
        c = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
        #catalog = SkyCoord(ra=gaiaStars[1:, 5]*u.degree, dec=gaiaStars[1:, 7]*u.degree)  
        catalog = SkyCoord(ra=gaiaStars['ra'], dec=gaiaStars['dec'])
        idx, d2d, d3d = c.match_to_catalog_sky(catalog)
        catalogstar = SkyCoord(ra=gaiaStars[idx]['ra']*u.degree, dec=gaiaStars[idx]['dec']*u.degree)
        sep = c.separation(catalogstar)
        if sep < Angle('00d00m02s'):
            sourceTable.add_row([fitsFile, gaiaStars[idx]['SOURCE_ID'], gaiaStars[idx]['DESIGNATION'], convertStringToNan(gaiaStars[idx]['ra']), convertStringToNan(gaiaStars[idx]['dec']), convertStringToNan(gaiaStars[idx]['parallax']), convertStringToNan(gaiaStars[idx]['parallax_error']), convertStringToNan(gaiaStars[idx]['pmra']), convertStringToNan(gaiaStars[idx]['pmdec']), convertStringToNan(gaiaStars[idx]['phot_g_mean_mag']), convertStringToNan(gaiaStars[idx]['phot_bp_mean_mag']), convertStringToNan(gaiaStars[idx]['phot_rp_mean_mag']), convertStringToNan(gaiaStars[idx]['radial_velocity']), convertStringToNan(gaiaStars[idx]['radial_velocity_error']), convertStringToNan(gaiaStars[idx]['teff_gspphot']), star['id'], ra2, dec2, star['mag'], fitsFileDate])

#gaiaStarsTableFileName = (str('gaiaStarsTab.csv'))

# Write found sources into file
tableFileName = (workingDirectory + '/' + str(fitsFile[:-4] + '.csv'))
sourceTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',')

### Search double stars on the image sequence
sourceTable_by_file = sourceTable.group_by('filename')
#print(sourceTable_by_file.info)
#print(sourceTable_by_file)

for group in sourceTable_by_file.groups:
    # Creating empty arrays for Star related calculations
    StarA = []
    StarB = []
    for star in group: ## modify according to arrays instead of starlist
        StarA = (star['filename'], star['source_id'], star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'],star['phot_g_mean_mag'], star['phot_bp_mean_mag'], star['phot_rp_mean_mag'], star['radial_velocity'], star['radial_velocity_error'], star['teff_gspphot'], star['image_id'], star['source_ra'], star['source_dec'], star['source_mag'], star['image_date'])
        for star in group:
            StarB = (star['filename'], star['source_id'], star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'],star['phot_g_mean_mag'], star['phot_bp_mean_mag'], star['phot_rp_mean_mag'], star['radial_velocity'], star['radial_velocity_error'], star['teff_gspphot'], star['image_id'], star['source_ra'], star['source_dec'], star['source_mag'], star['image_date'])
            if StarA != StarB and float(StarA[9]) < float(StarB[9]) and float(StarA[5]) != 0 and float(StarB[5]) != 0:
                #Set input data
                starRa1 = float(StarA[3])
                starDec1 = float(StarA[4])
                starRa2 = float(StarB[3])
                starDec2 = float(StarB[4])
                starCoord1 = SkyCoord(StarA[3], StarA[4], unit="deg")
                starCoord2 = SkyCoord(StarB[3], StarB[4], unit="deg")
                starParallax1 = float(StarA[5])
                starParallaxError1 = float(StarA[6])
                            
                # Calculate the widest possible separation for StarA
                possSep1 = possible_distance / calcDistanceMax(starParallax1, starParallaxError1)
                rhoStar = starCoord1.separation(starCoord2).arcsecond
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
                    starActualMag2 = float(StarB[18])
                    starActualCoord1 = SkyCoord(starActualRa1, starActualDec1, unit="deg")
                    starActualCoord2 = SkyCoord(starActualRa2, starActualDec2, unit="deg")
                    starObjectId = (str(starId1) + '_' + str(starId2))
                    starImageDate = StarA[19]

                    # Value to modify Theta according to the appropriate quadrant
                    '''addThetaValue = ()
                    if deltaRa(starRa1, starRa2, starDec2) > 0 and deltaDec(starDec2, starDec1) > 0:
                        addThetaValue = 0
                    elif deltaRa(starRa1, starRa2, starDec2) > 0 and deltaDec(starDec2, starDec1) < 0:
                        addThetaValue = 180
                    elif deltaRa(starRa1, starRa2, starDec2) < 0 and deltaDec(starDec2, starDec1) < 0:
                        addThetaValue = 180
                    elif deltaRa(starRa1, starRa2, starDec2) < 0 and deltaDec(starDec2, starDec1) > 0:
                        addThetaValue = 360'''
                    
                    # Calculate actual data based on functions
                    #thetaStar = thetaCalc(deltaRa(starRa1, starRa2, starDec2), deltaDec(starDec2, starDec1)) + addThetaValue
                    thetaStar = starCoord1.position_angle(starCoord2).degree
                    #thetaActual = thetaCalc(deltaRa(starActualRa1, starActualRa2, starActualDec2), deltaDec(starActualDec2, starActualDec1)) + addThetaValue
                    thetaActual = starActualCoord1.position_angle(starActualCoord2).degree
                    #rhoActual = rhoCalc(starActualRa1, starActualDec1, starActualRa2, starActualDec2)
                    rhoActual = starActualCoord1.separation(starActualCoord2).arcsecond
                    starDistance1 = calcDistance(starParallax1)
                    starDistanceMax1 = calcDistanceMax(starParallax1, starParallaxError1)
                    starDistanceMin1 = calcDistanceMin(starParallax1, starParallaxError1)
                    starDistanceRange1 = starDistanceMax1 - starDistanceMin1
                    starDistance2 = calcDistance(starParallax2)
                    starDistanceMax2 = calcDistanceMax(starParallax2, starParallaxError2)
                    starDistanceMin2 = calcDistanceMin(starParallax2, starParallaxError2)
                    starDistanceRange2 = starDistanceMax2 - starDistanceMin2

                    # Check if stars shares a common distance range
                    distanceCommon = ()
                    if starDistanceMin1 < starDistanceMin2 < starDistanceMax1 or starDistanceMin2 < starDistanceMin1 < starDistanceMax2:
                        distanceCommon = 'overlapping'
                    else:
                        distanceCommon = 'no'
                    
                    # Check if the pair is a Common Proper Motion pairs (CPM), Similar Proper Motion (SPM) or Different Proper Motion (DPM)


                    #Print data, if stars are close and share a common distance range
                    if distanceCommon == 'overlapping':
                        reportTable.add_row([star[0], starId1, starName1, starRa1, starDec1, starParallax1, starParallaxError1, starPmRa1, starPmDec1, starGMag1, starBpMag1, starRpMag1, starRadVel1, starRadVelErr1, starTemp1, starActualRa1, starActualDec1, starActualMag1, starId2, starName2, starRa2, starDec2, starParallax2, starParallaxError2, starPmRa2, starPmDec2, starGMag2, starBpMag2, starRpMag2, starRadVel2, starRadVelErr2, starTemp2, starActualRa2, starActualDec2, starActualMag2, thetaStar, thetaActual, rhoStar, rhoActual, starObjectId, starImageDate])

# print('### Report Table ###')
# print(reportTable)

# Create Qtable to list the measurements and the standard error per groups (double star)
measuredObject = np.array([], dtype=str)
measuredMeanTheta = np.array([], dtype=np.float64)
measuredMeanThetaErr = np.array([], dtype=np.float64)
measuredMeanRho = np.array([], dtype=np.float64)
measuredMeanRhoErr = np.array([], dtype=np.float64)
meanTable = QTable([measuredObject, measuredMeanTheta, measuredMeanThetaErr, measuredMeanRho, measuredMeanRhoErr], names=('object_name', 'mean_theta', 'mean_theta_error', 'mean_rho', 'mean_rho_error'), meta={'name': 'measured table'}) # measuredArrayTheta, measuredArrayRho, # 'measurements_theta', 'measurements_rho', measuredStarA, measuredStarB, 'star_a', 'star_b', 


### Search double stars on the image sequence
reportTable_by_object = reportTable.group_by('object_id')
print('\n### Report Table by object ###')
print(reportTable_by_object.info)

objectMean = reportTable_by_object.groups.aggregate(np.mean)

count = 1
for ds in reportTable_by_object.groups:
    print('\n### Group index:', count, '###')
    count = count + 1
    rhoPairDr3 = rhoCalc(ds[0][3], ds[0][4], ds[0][20], ds[0][21])
    pairFileName = ds[0][0]
    pairRaA = ds[0][3]
    pairDecA = ds[0][4]
    pairRaB = ds[0][20]
    pairDecB = ds[0][21]
    pairDistanceMinA = calcDistanceMin(ds[0][5], ds[0][6])
    pairDistanceMinB = calcDistanceMin(ds[0][22], ds[0][23])
    pairMeanTheta = ds['theta_measured'].groups.aggregate(np.mean)
    pairMeanThetaErr = ds['theta_measured'].groups.aggregate(np.std)
    pairMeanRho = ds['rho_measured'].groups.aggregate(np.mean)
    pairMeanRhoErr = ds['rho_measured'].groups.aggregate(np.std)
    pairMagMeasuredA = ds['magmeasured_a'].groups.aggregate(np.mean)
    pairMagMeasuredAErr = ds['magmeasured_a'].groups.aggregate(np.std)
    pairMagMeasuredB = ds['magmeasured_b'].groups.aggregate(np.mean)
    pairMagMeasuredBErr = ds['magmeasured_b'].groups.aggregate(np.std)
    pairDesignationA = ds[0][2]
    pairDesignationB = ds[0][19]
    pairGMagnitudeA = ds[0][9]
    pairGMagnitudeB = ds[0][26]
    pairMagDiff = math.fabs(pairMagMeasuredA - pairMagMeasuredB)
    pairMagDiffDr3 = math.fabs(pairGMagnitudeA - pairGMagnitudeB)
    pairMagDiffErr = math.fabs(((ds['magmeasured_a']) - (ds['magmeasured_b'])).groups.aggregate(np.std))
    pairRadVelRatioA = math.fabs(ds[0][13] / ds[0][12]) * 100
    pairRadVelRatioB = math.fabs(ds[0][30] / ds[0][29]) * 100
    pairRadVelA, pairRadVelAErr, pairRadVelB, pairRadVelBErr = ds[0][12], ds[0][13], ds[0][29], ds[0][30]
    pairParallaxFactor = (calcParallaxFactor(ds[0][5], ds[0][22])) * 100
    pairPmFactor = (calcPmFactor(ds[0][7], ds[0][8], ds[0][24], ds[0][25])) * 100
    pairPmCommon = calcPmCategory(pairPmFactor)
    pairAbsMag1 = calcAbsMag(pairGMagnitudeA, ds[0][5]) # Calculate Absolute magnitude
    pairAbsMag2 = calcAbsMag(pairGMagnitudeB, ds[0][22]) # Calculate Absolute magnitude
    pairLum1 = calcLuminosity(pairAbsMag1)
    pairLum2 = calcLuminosity(pairAbsMag2)
    pairAltLum1 = calcLuminosityAlternate(pairAbsMag1)
    pairAltLum2 = calcLuminosityAlternate(pairAbsMag2)
    pairMass1 = calcMass(pairAltLum1)
    pairMass2 = calcMass(pairAltLum2)
    pairBVIndexA = ds[0][10] - ds[0][11]
    pairBVIndexB = ds[0][27] - ds[0][28]
    pairSepPar = sepCalc(pairDistanceMinA, pairDistanceMinB, rhoPairDr3) # Separation of the pairs in parsecs
    pairEscapeVelocity = calcEscapevelocity(pairMass1, pairMass2, pairSepPar, gravConst)
    pairRelativeVelocity = calcRelativeVelocity(ds[0][7], ds[0][8], ds[0][24], ds[0][25], ds[0][12], ds[0][29], pairDistanceMinA, pairDistanceMinB)
    pairHarshawFactor = calcHarshaw(pairParallaxFactor, pairPmFactor)
    pairHarshawPhysicality = calcHarshawPhysicality(pairHarshawFactor)
    pairBinarity = calcBinarity(pairRelativeVelocity, pairEscapeVelocity)
    dateOfObservation = getUTC(Time(ds['image_date']).mean())

    pairACurrentCoord = calcCurrentDR3Coord(dateOfObservation, pairRaA, pairDecA, ds[0][7], ds[0][8])
    pairBCurrentCoord = calcCurrentDR3Coord(dateOfObservation, pairRaB, pairDecB, ds[0][24], ds[0][25])
    pairAMeasuredCoord = SkyCoord(ra=ds['rameasured_a'].groups.aggregate(np.mean) * u.deg, dec=ds['decmeasured_a'].groups.aggregate(np.mean) * u.deg)
    pairBMeasuredCoord = SkyCoord(ra=ds['rameasured_b'].groups.aggregate(np.mean) * u.deg, dec=ds['decmeasured_b'].groups.aggregate(np.mean) * u.deg)
    pairACoordErr = pairACurrentCoord.separation(pairAMeasuredCoord)
    pairBCoordErr = pairBCurrentCoord.separation(pairBMeasuredCoord)

    preciseCoord = str(getPreciseCoord(pairRaA, pairDecA, fitsFileDate))

    pairNumTheta = len(ds['theta_measured'])
    pairNumRho = len(ds['rho_measured'])
    pairNumMagMeasureadA = len(ds['magmeasured_a'])
    pairNumMagMeasureadB = len(ds['magmeasured_b'])

    reportName = (workingDirectory + '/' + ds[0][39] + '.txt')
    reportFile = open(reportName, "a")

    hrdPlot(pairDesignationA, pairDesignationB, pairAbsMag1, pairAbsMag2, pairBVIndexA, pairBVIndexB)
    imagePlot(pairFileName, pairDesignationA, pairDesignationB, pairRaA, pairDecA, pairRaB, pairDecB)

    print('### COMPONENTS ###')
    print('\nDate of observation: ' + dateOfObservation)
    print('\nPrecise coordinates (J2000): ' + preciseCoord)
    print('\nComponent A:', pairDesignationA)
    print('Component B:', pairDesignationB)
    print('\nCalculated coordinates')
    print('\nComponent A DR3 2016:', pairRaA, pairDecA)
    print('Component A DR3 on date:', pairACurrentCoord.ra.degree, pairACurrentCoord.dec.degree)
    print('Component A measured:', pairAMeasuredCoord.ra.degree, pairAMeasuredCoord.dec.degree)
    print('Component A error (on date - measured):', pairACoordErr.arcsecond)
    print('\nComponent B DR3 2016:', pairRaB, pairDecB)
    print('Component B DR3 on date:', pairBCurrentCoord.ra.degree, pairBCurrentCoord.dec.degree)
    print('Component B measured:', pairBMeasuredCoord.ra.degree, pairBMeasuredCoord.dec.degree)
    print('Component B error (on date - measured):', pairBCoordErr.arcsecond)
    print('2016 Calculared Position angle / Separation: ', pairACurrentCoord.position_angle(pairBCurrentCoord).degree, pairACurrentCoord.separation(pairBCurrentCoord).arcsecond)
    print('Current Calculared Position angle / Separation: ', SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).degree, SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).arcsecond)
    print('\nTheta measurements\n', ds['theta_measured'])
    print('Number of measurements: ', pairNumTheta)
    print('Mean:', pairMeanTheta[0])
    print('Error:', pairMeanThetaErr[0])
    print('\nRho measurements\n', ds['rho_measured'])
    print('Number of measurements: ', pairNumRho)
    print('Mean:', pairMeanRho[0])
    print('Error:', pairMeanRhoErr[0])
    print('\nMagnitude A DR3: \n', pairGMagnitudeA)
    print('\nMagnitude A measurements\n', ds['magmeasured_a'])
    print('Number of measurements: ', pairNumMagMeasureadA)
    print('Mean:', pairMagMeasuredA[0])
    print('Error:', pairMagMeasuredAErr[0])
    print('\nMagnitude B DR3: \n', pairGMagnitudeB)
    print('\nMagnitude B measuremets\n', ds['magmeasured_b'])
    print('Number of measurements: ', pairNumMagMeasureadB)
    print('Mean:', pairMagMeasuredB[0])
    print('Error:', pairMagMeasuredBErr[0])
    print('\nMagnitude difference (DR3):', pairMagDiffDr3)
    print('Magnitude difference (measured):', pairMagDiff)
    print('Parallax factor:', pairParallaxFactor, '%')
    print('Proper motion factor:', pairPmFactor, '%')
    print('Proper motion category:', pairPmCommon)
    print('Absolute magnitude A:', pairAbsMag1)
    print('Absolute magnitude B:', pairAbsMag2)
    print('Luminosity A:', pairLum1)
    print('Luminosity B:', pairLum2)
    print('Luminosity Alternate A:', pairAltLum1)
    print('Luminosity Alternate B:', pairAltLum2)
    print('Mass A:', pairMass1)
    print('Mass B:', pairMass2)
    print('BV index A:', pairBVIndexA, 'B:', pairBVIndexB)
    print('Radial velocity of the stars', 'A:', pairRadVelA, 'km/s (Err:', pairRadVelAErr, 'km/s)', 'B:', pairRadVelB, 'km/s (Err:', pairRadVelBErr, 'km/s)')
    print('Radial velocity ratio A:', pairRadVelRatioA, '%')
    print('Radial velocity ratio B:', pairRadVelRatioB, '%')
    print('Separation:', pairSepPar, 'parsec,', pairSepPar * 206265, 'AU')
    print('Pair Escape velocity:', pairEscapeVelocity, 'km/s')
    print('Pair Relative velocity:', pairRelativeVelocity, 'km/s')
    print('Pair Harshaw factor:', pairHarshawFactor)
    print('Pair Harshaw physicality:', pairHarshawPhysicality)
    print('Pair binarity:', pairBinarity)
    
    reportFile.write('### COMPONENTS ###')
    reportFile.write('\nDate of observation: ' + dateOfObservation)
    reportFile.write('\nPrecise coordinates (J2000): ' + preciseCoord)     
    reportFile.write('\nComponent A: ' + pairDesignationA)
    reportFile.write('\nComponent B: ' + pairDesignationB)
    reportFile.write('\n\nCalculated coordinates')
    reportFile.write('\nComponent A DR3 2016: ' + str(pairRaA) + ' ' + str(pairDecA))
    reportFile.write('\nComponent A DR3 on date: ' + str(pairACurrentCoord.ra.degree) + ' ' + str(pairACurrentCoord.dec.degree))
    reportFile.write('\nComponent A measured: ' + str(pairAMeasuredCoord.ra.degree) + ' ' + str(pairAMeasuredCoord.dec.degree))
    reportFile.write('\nComponent A error (on date - measured): ' + str(pairACoordErr.arcsecond))
    reportFile.write('\nComponent B DR3 2016: ' + str(pairRaB) + ' ' + str(pairDecB))
    reportFile.write('\nComponent B DR3 on date: ' + str(pairBCurrentCoord.ra.degree) + ' ' + str(pairBCurrentCoord.dec.degree))
    reportFile.write('\nComponent B measured: ' + str(pairBMeasuredCoord.ra.degree) + ' ' + str(pairBMeasuredCoord.dec.degree))
    reportFile.write('\nComponent B error (on date - measured): ' + str(pairBCoordErr.arcsecond))
    reportFile.write('\n\n2016 Calculared Position angle / Separation: '  + str(pairACurrentCoord.position_angle(pairBCurrentCoord).degree) + ' ' + str(pairACurrentCoord.separation(pairBCurrentCoord).arcsecond))
    reportFile.write('\nCurrent Calculared Position angle / Separation: ' + str(SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).degree) + ' ' + str(SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).arcsecond))
    reportFile.write('\n\nTheta measurements\n' + str(ds['theta_measured']))
    reportFile.write('\nMean: ' + str(pairMeanTheta[0]))
    reportFile.write('\nError: ' + str(pairMeanThetaErr[0]))
    reportFile.write('\n\nRho measurements\n' + str(ds['rho_measured']))
    reportFile.write('\nMean: ' + str(pairMeanRho[0]))
    reportFile.write('\nError: ' + str(pairMeanRhoErr[0]))
    reportFile.write('\n\nMagnitude A DR3:  \n' + str(pairGMagnitudeA))
    reportFile.write('\n\nMagnitude A measurements\n' + str(ds['magmeasured_a']))
    reportFile.write('\nMean: ' + str(pairMagMeasuredA[0]))
    reportFile.write('\nError: ' + str(pairMagMeasuredAErr[0]))
    reportFile.write('\n\nMagnitude B DR3:  \n' + str(pairGMagnitudeB))
    reportFile.write('\n\nMagnitude B measuremets\n' + str(ds['magmeasured_b']))
    reportFile.write('\nMean: ' + str(pairMagMeasuredB[0]))
    reportFile.write('\nError: ' + str(pairMagMeasuredBErr[0]))
    reportFile.write('\n\nMagnitude difference (DR3): ' + str(pairMagDiffDr3))
    reportFile.write('\nMagnitude difference (measured): ' + str(pairMagDiff))
    reportFile.write('\n\nParallax factor: ' + str(pairParallaxFactor) + '%')
    reportFile.write('\nProper motion factor: ' + str(pairPmFactor) + '%')
    reportFile.write('\nProper motion category: ' + str(pairPmCommon))
    reportFile.write('\nAbsolute magnitude A: ' + str(pairAbsMag1))
    reportFile.write('\nAbsolute magnitude B: ' + str(pairAbsMag2))
    reportFile.write('\nLuminosity A: ' + str(pairLum1))
    reportFile.write('\nLuminosity B: ' + str(pairLum2))
    reportFile.write('\nMass A: ' + str(pairMass1))
    reportFile.write('\nMass B: ' + str(pairMass2))
    reportFile.write('\nBV index A: ' + str(pairBVIndexA) + ' B: ' + str(pairBVIndexB))
    reportFile.write('\nRadial velocity of the stars A: ' + str(pairRadVelA) + ' km/s (Err: ' + str(pairRadVelAErr) + ' km/s) B: ' + str(pairRadVelB) + ' km/s (Err: ' + str(pairRadVelBErr) + ' km/s)')
    reportFile.write('\nRadial velocity ratio A: ' + str(pairRadVelRatioA) + ' %')
    reportFile.write('\nRadial velocity ratio B: ' + str(pairRadVelRatioB) + ' %')
    reportFile.write('\nSeparation: ' + str(pairSepPar) + ' parsec, ' + str(pairSepPar * 206265) + ' AU')
    reportFile.write('\nPair Escape velocity: ' + str(pairEscapeVelocity) + ' km/s')
    reportFile.write('\nPair Relative velocity: ' + str(pairRelativeVelocity) + ' km/s')
    reportFile.write('\nPair Harshaw factor: ' + str(pairHarshawFactor))
    reportFile.write('\nPair Harshaw physicality: ' + str(pairHarshawPhysicality))
    reportFile.write('\nPair binarity: ' + str(pairBinarity))

    reportFile.write('\n\n### Publication table 1. ###')
    reportFile.write('\n' + str(pairDesignationA) + ',' + str(pairMagMeasuredA[0]) + ',' + str(pairMagMeasuredAErr[0]) + ',' + str(pairGMagnitudeA) + ',' + str(pairAbsMag1) + ',' + str(pairLum1) + ',' + str(pairMass1) + ',' + dateOfObservation)
    reportFile.write('\n' + str(pairDesignationB) + ',' + str(pairMagMeasuredB[0]) + ',' + str(pairMagMeasuredBErr[0]) + ',' + str(pairGMagnitudeB) + ',' + str(pairAbsMag2) + ',' + str(pairLum2) + ',' + str(pairMass2) + ',' + dateOfObservation)
    
    reportFile.write('\n\n### Publication table 2. ###')
    reportFile.write('\n' + preciseCoord + ',' + str(pairMeanTheta[0]) + ',' + str(pairMeanThetaErr[0]) + ',' + str(pairMeanRho[0]) + ',' + str(pairMeanRhoErr[0]) + ',' + str(pairSepPar * 206265) + ' AU')
 
    reportFile.write('\n\n### Publication table 3. ###')
    reportFile.write('\n' + str(pairParallaxFactor) + '%' + ',' + str(pairPmFactor) + '%' + ',' + str(pairPmCommon) + ',' + str(pairEscapeVelocity) + ',' + str(pairRelativeVelocity) + ',' + str(pairHarshawFactor) + ',' + str(pairHarshawPhysicality) + ',' + str(pairBinarity))

    reportFile.write('\n\n### WDS table ###')
    c = str(SkyCoord(pairRaA, pairDecA, unit='deg', frame='icrs').to_string('hmsdms'))
    pair_wds_name = str(c[0:2]) + str(c[3:5]) + str(c[6:8]) + str(c[9:10]) + str(c[19:22]) + str(c[23:25]) + str(c[26:28]) + str(c[29:30])
    wdsform = pair_wds_name + ',' + dateOfObservation + ',' +  str(pairMeanTheta[0]) + ',' + str(pairMeanThetaErr[0]) + ',' + str(pairMeanRho[0]) + ',' + str(pairMeanRhoErr[0]) + ',' +  'nan' + ',' +  'nan' + ',' +  str(pairMagDiff) + ',' +  str(pairMagDiffErr) + ',' + 'Filter wawelenght' + ',' + 'filter FWHM' + ',' + '0.2' + ',' + '1' + ',' + 'TLB_2023' + ',' +  'C' + ',' + '7' + ',' + preciseCoord.replace(":","").replace(" ","")

    reportFile.write('\n\n### Publication table sum ###')
    reportFile.write('\n' + str(pairDesignationA) + ',' + str(pairMagMeasuredA[0]) + ',' + str(pairMagMeasuredAErr[0]) + ',' + str(pairGMagnitudeA) + ',' + str(pairAbsMag1) + ',' + str(pairLum1) + ',' + str(pairMass1) + ',' + dateOfObservation + ',' + str(pairDesignationB) + ',' + str(pairMagMeasuredB[0]) + ',' + str(pairMagMeasuredBErr[0]) + ',' + str(pairGMagnitudeB) + ',' + str(pairAbsMag2) + ',' + str(pairLum2) + ',' + str(pairMass2) + ',' + preciseCoord + ',' + str(pairMeanTheta[0]) + ',' + str(pairMeanThetaErr[0]) + ',' + str(pairMeanRho[0]) + ',' + str(pairMeanRhoErr[0]) + ',' + str(pairSepPar * 206265) + ' AU' + ',' + str(pairParallaxFactor) + '%' + ',' + str(pairPmFactor) + '%' + ',' + str(pairPmCommon) + ',' + str(pairEscapeVelocity) + ',' + str(pairRelativeVelocity) + ',' + str(pairHarshawFactor) + ',' + str(pairHarshawPhysicality) + ',' + str(pairBinarity) + ',' + pair_wds_name + ',' + dateOfObservation + ',' +  str(pairMeanTheta[0]) + ',' + str(pairMeanThetaErr[0]) + ',' + str(pairMeanRho[0]) + ',' + str(pairMeanRhoErr[0]) + ',' +  'nan' + ',' +  'nan' + ',' +  str(pairMagDiff) + ',' +  str(pairMagDiffErr) + ',' + 'Filter wawelenght' + ',' + 'filter FWHM' + ',' + '0.2' + ',' + '1' + ',' + 'TLB_2023' + ',' +  'C' + ',' + '7' + ',' + preciseCoord.replace(":","").replace(" ",""))

    reportFile.write('\n' + str(wdsform))
    reportFile.close()
