#! /usr/bin/python3
# WDS Report tool to measure double stars on astronomical images based on Gaia DR3 data
# Version: 1.0
# Usage: wdsreport <image_folder>

import os
import sys
import numpy as np
import math
import sys
import datetime
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy.coordinates import match_coordinates_sky, search_around_sky, position_angle, angular_separation
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack, hstack
from astropy.table import Column, MaskedColumn
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
import warnings
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta
warnings.filterwarnings("ignore")
from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
 
# Configuration for the ATIK camera

dao_sigma = 3.0
dao_fwhm = 7.0
dao_threshold = 12.0
possible_distance = 10000.0 # AU
search_cone = 0.001 # Decimal degree
dummyObservationDate = "2022-01-01T12:00:00"


# Configuration for the CANON camera

'''dao_sigma = 2.0
dao_fwhm = 3.0
dao_threshold = 5.0
possible_distance = 10000.0 # AU
search_cone = 0.001 # Decimal degree'''

# Gravitational constant is convenient if measure distances in parsecs (pc), velocities in kilometres per second (km/s) and masses in solar units M
gravConst = 0.0043009 
image_limit = 2000

# Constant to calculate star luminosity and mass
sun_luminosity = 3.0128 * (10 ** 28)
sun_absolute_luminosity = 3.828 * (10 ** 26)

# Constant variables
hipparcos_file = Table.read(f"/usr/share/dr3map/hipparcos/I_239_selection.csv", format='ascii')

# Insert the downloaded wds file path here
wds_file = "/usr/share/dr3map/wds/wdsweb_summ2.txt"

# WDS table to be used to identify double stars on the image
# wdsTable = Table.read(sys.argv[1], delimiter=',', format='ascii')

#wdsTable = Table.read(f"/usr/share/dr3map/dr3-wds/wdsweb_summ2.dat", format='ascii')
#print(wdsTable.info)


# Create WDS table
wds_converters = {  '2000 Coord': np.str_,
                    'Discov': np.str_,
                    'Comp': np.str_,
                    'Date (first)': np.str_,
                    'Date (last)': np.str_,
                    'Obs': np.str_,
                    'PA_f': np.float64,
                    'PA_l': np.float64,
                    'Sep_f': np.float64,
                    'Sep_l': np.float64,
                    'Mag_A': np.str_,
                    'Mag_B': np.str_,
                    'Spectral A/B': np.str_,
                    'PM_A_ra': np.str_,
                    'PM_A_dec': np.str_,
                    'PM_B_ra': np.str_,
                    'PM_B_dec': np.str_,
                    'D. Number': np.str_,
                    'Notes': np.str_,
                    'Coord (RA)': np.str_,
                    'Coord (DEC)': np.str_
                    }
wds_data = Table.read(wds_file,
                      names=('2000 Coord', 'Discov', 'Comp', 'Date (first)', 'Date (last)', 'Obs',
                             'PA_f', 'PA_l', 'Sep_f', 'Sep_l', 'Mag_A',
                             'Mag_B', 'Spectral A/B', 'PM_A_ra', 'PM_A_dec',
                             'PM_B_ra', 'PM_B_dec', 'D. Number', 'Notes', 'Coord (RA)',
                             'Coord (DEC)'
                        ),
                      converters=wds_converters,
                      format='ascii.fixed_width',
                      header_start=2, data_start=5,
                      col_starts=(0, 10, 17, 23, 28, 33, 38, 42, 46, 52, 58, 64, 70, 80, 84, 89, 93, 98, 107, 112, 121),
                      col_ends=(9, 16, 21, 26, 31, 36, 40, 44, 50, 56, 61, 68, 78, 83, 87, 92, 96, 105, 110, 120, 129),
                      )

# Set working directory to read Double Star images
workingDirectory = sys.argv[1]

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

def convertStringToNan(str):
    string = np.nan
    if str == 'null' or str == '' or str == '.':
        string = np.nan
    else:
        string = str
    return string

def roundNumber(num):
    if type(num) == np.float32 or type(num) == np.float64 or type(num) == float or type(num) == int:
        finalNumber = np.round(num, 3)
    else:
        finalNumber = num
    return finalNumber

'''def rhoCalc(raa, deca, rab, decb):
    rhocalc = math.sqrt(((raa-rab) * math.cos(math.radians(deca))) ** 2 + (deca - decb) ** 2) * 3600
    return rhocalc'''

# Function to calculate Delta RA
def deltaRa(raa, rab, decb):
    deltara = (rab - raa) * math.cos(math.radians(decb))
    return deltara

# Function to calculate Delta DEC
def deltaDec(decb, deca):
    deltadec = decb - deca
    return deltadec

# Function to calculate Theta (position angle)
'''def thetaCalc(deltara,deltadec):
    if (deltara != 0 and deltadec != 0):
        thetacalc = math.degrees(math.atan(deltara/deltadec))
    else:
        thetacalc = 0
    return thetacalc'''

# Function to calculate Rho (separation)
'''def rhoCalc(raa, deca, rab, decb):
    rhocalc = math.sqrt(((raa-rab) * math.cos(math.radians(deca))) ** 2 + (deca - decb) ** 2) * 3600
    return rhocalc'''

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
    dist = 1 / (math.fabs(par/1000))
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

# Calculate the radius of the star
# Excel formula =SQRT(Luminosity/(T eff/5778))
def calcRadius(luminosity, teffStar):
    if bool(luminosity) and bool(teffStar):
        teffSun = 5778
        radius = math.sqrt(luminosity / ((teffStar / teffSun) ** 4))
    else:
        radius = 'Cannot be calculated, T(eff) or Luminosity is missing'
    return radius

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
    if bool(mass_a) and bool(mass_b) and bool(separation) and bool(gravconst):
        print('calcEscapevelocity: ' + str(mass_a) + ' ' + str(mass_b) + ' ' + str(separation) + ' ' + str(gravconst))
        escvel = math.sqrt((2 * gravconst * (mass_a + mass_b)) / separation)
    else:
        escvel = 0.0
    return escvel

# Function to calculate the Probability of binarity based on the Relative and the escape velocity
# Excel formula =IF(relative speed<=escape velocity,"Y","No"))
def calcBinarity(relsped, escsped):
    binarity = ()
    #if relsped == 'nan' or escsped == 'nan':
    if math.isnan(relsped) or math.isnan(escsped):
        binarity = 'Missing data, binarity cannot be determined.'
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
    if pmfact >= 0.8:
        pmCommon = 'CPM'
    elif 0.4 <= pmfact < 0.8:
        pmCommon = 'SPM'
    elif pmfact < 0.4:
        pmCommon = 'DPM'
    return pmCommon

# Search pair in Washington double Star Catalog
'''
def searchWds(pairadesig):
    wdsIdx = np.where(pairadesig == wdsFile['designation'])
    wdsRow = wdsFile[wdsIdx]
    if wdsRow:
        wdsPair = wdsRow[0]
    elif not wdsRow:
        wdsPair = 'Not found'
    return wdsPair
'''
    
# Function to search coordinates in the WDS catalog file
def search_in_wds(ra_source, dec_source):
    coordinates = SkyCoord(ra=ra_source, dec=dec_source)
    idx, d2d, d3d = coordinates.match_to_catalog_sky(wds_catalog)
    star_coord = SkyCoord(ra=ra_source, dec=dec_source)
    print('Coordinates: ' + str(ra_source) + '(Ra), ' + str(dec_source) + '(Dec), Separation: ' + str(d2d))
    #sep = coordinates.separation(d2d)*u.degree
    print(wdsTable[idx]['2000 Coord'])
    print(wdsTable[np.where(wds_catalog['2000 Coord'] == wdsTable[idx]['2000 Coord'])])

def create_unique_id(wds_id, discov):
    id_array = [i + j for i, j in zip(wds_id, discov)]
    return id_array

def delete_invalid_lines_wds(catalog):
    rows_to_delete = np.where(catalog['Coord (RA) hms'] == '.hms')
    catalog.remove_rows(rows_to_delete)
    return catalog

def calculate_wds_ra_hourangle(wds_ra_array):
    #print(wds_ra_array)
    wds_ra_hms = []
    for star in wds_ra_array:
        wds_ra_hms.append(str(star[0:2]) + 'h' + str(star[2:4]) + 'm' + str(star[4:9]) + 's')
    return wds_ra_hms

def calculate_wds_dec_hourangle(wds_dec_array):
    #print(wds_dec_array)
    wds_dec_dms = []
    for star in wds_dec_array:
        wds_dec_dms.append(str(star[0:3]) + 'd' + str(star[3:5]) + 'm' + str(star[5:9]) + 's')
    return wds_dec_dms

# Create HRD plot of the double stars based on Hipparcos
def hrdPlot(pairname, mag_abs_a, mag_abs_b, bv_a, bv_b):
    print(pairname, mag_abs_a, mag_abs_b, bv_a, bv_b)
    if pairname and mag_abs_a and mag_abs_b and bv_a and bv_b:
        hipparcos_abs_mag = hipparcos_file['Abs_mag']
        hipparcos_bv_index = hipparcos_file['B-V']
        plt.scatter(hipparcos_bv_index, hipparcos_abs_mag, s=0.5, alpha=0.2, color="grey") #, 
        plt.scatter(bv_a, mag_abs_a, s=14, color="blue", label='Main star') # s= 1 / mag_abs_a
        plt.scatter(bv_b, mag_abs_b, s=7, color="red", label='Companion star') # s= 1 / mag_abs_a
        plt.legend(loc="upper left")
        plt.axis((-0.4,1.9,21,-16))
        plt.title('Double Star ' + pairname + ' H-R Diagram')
        plt.xlabel('B-V index')
        plt.ylabel('Absolute magnitude')
        plt.gca().set_aspect(0.07)
        savename = str(workingDirectory + '/' + pairname + '_hrd.jpg')
        plt.savefig(savename, bbox_inches='tight', dpi=150.0)
        plt.close()
    else:
        print('Data is missiong, HRD plot cannot be created!')


# Create Image plot of the double stars
def imagePlot(filename, pairname, raa, deca, rab, decb):
    '''coord_meta = {}
    coord_meta['type'] = ('longitude', 'latitude')
    coord_meta['wrap'] = (None, None)
    coord_meta['unit'] = (u.degree, u.degree)
    coord_meta['format_unit'] = (u.hour, u.degree)
    coord_meta['name'] = 'ra', 'dec'
    '''
    
    image_data = fits.open(workingDirectory + '/' + filename)
    header = image_data[0].header
    wcs_helix = WCS(image_data[0].header)
    image = image_data[0].data
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

    #plt.title(pairname, pad=50.0)
    plt.title(pairname)

    plt.imshow(image, origin='lower',cmap='grey', aspect='equal', vmax=image_limit, vmin=0) # , cmap='cividis'
    plt.savefig(workingDirectory + '/' + pairname + '_img.jpg',dpi=150.0, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    #plt.show()

def get_gaia_dr3_data(doublestars):
    pairACoord = SkyCoord(ra=doublestars[0]['ra_deg_1'], dec=doublestars[0]['dec_deg_1'], unit=(u.degree, u.degree), frame='icrs')
    pairBCoord = SkyCoord(ra=doublestars[0]['ra_deg_2'], dec=doublestars[0]['dec_deg_2'], unit=(u.degree, u.degree), frame='icrs')
    a = Gaia.cone_search_async(pairACoord, radius=u.Quantity(search_cone, u.deg))
    b = Gaia.cone_search_async(pairBCoord, radius=u.Quantity(search_cone, u.deg))
    a_query = a.get_results()
    b_query = b.get_results()
    return a_query, b_query

###################################################################################################################################

print('\n### Reading WDS database ###')

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

# Define Report table
reportw_identifier = np.array([], dtype=str)
reportdate = np.array([], dtype=str)
reporttheta = np.array([], dtype=np.float64)
reportthetaerr = np.array([], dtype=np.float64)
reportrho = np.array([], dtype=np.float64)
reportrhoerr = np.array([], dtype=np.float64)
reportmag_pri = np.array([], dtype=np.float64)
reportmag_prierr = np.array([], dtype=np.float64)
reportmag_sec = np.array([], dtype=np.float64)
reportmag_secerr = np.array([], dtype=np.float64)
reportfilter =  np.array([], dtype=str)
reportfilterfwhm = np.array([], dtype=str)
reporttelescopeap = np.array([], dtype=str)
reportnights = np.array([], dtype=str)
reportrefcode = np.array([], dtype=str)
reporttech = np.array([], dtype=str)
reportcat = np.array([], dtype=str)
preccoord = np.array([], dtype=str)
reportTable = QTable([reportw_identifier, reportdate, reporttheta, reportthetaerr, reportrho, reportrhoerr, reportmag_pri, reportmag_prierr, reportmag_sec, reportmag_secerr, reportfilter, reportfilterfwhm, reporttelescopeap, reportnights, reportrefcode, reporttech, reportcat, preccoord], names=('wds_identifier', 'date_of_obs', 'mean_theta', 'mean_theta_err', 'mean_rho', 'mean_rho_err', 'mag_pri', 'mag_pri_err', 'mag_sec', 'mag_sec_err', 'filter', 'filter_fwhm', 'telescope_ap', 'nights_of_obs', 'reference_code', 'tech_code', 'catalog_code', 'precise_cordinates'), meta={'name': 'report table'})

print('\n### Creating filelist ###')

directoryContent = os.listdir(workingDirectory)
print('Working directory: ', workingDirectory)

# Add fit, fits file extensions too
files = [f for f in directoryContent if os.path.isfile(workingDirectory+'/'+f) and (f.endswith('.new') or f.endswith('.fit') or f.endswith('.fits'))]
print('Files:', files)

# Create the WDS catalog table

print('Creating wds catalog - ' , datetime.datetime.now())

wdsTable = hstack([wds_data, calculate_wds_ra_hourangle(wds_data['Coord (RA)'])])
wdsTable.rename_column('col0', 'Coord (RA) hms')
wdsTable = hstack([wdsTable, calculate_wds_dec_hourangle(wds_data['Coord (DEC)'])])
wdsTable.rename_column('col0', 'Coord (DEC) dms')
wdsTable = hstack([wdsTable, create_unique_id(wds_data['2000 Coord'], wds_data['Discov'])])
wdsTable.rename_column('col0', 'Unique ID')
wdsTable = delete_invalid_lines_wds(wdsTable)

print('WDS angle formatting done, creating catalog... - ' , datetime.datetime.now())

wds_catalog = SkyCoord(ra=wdsTable['Coord (RA) hms'], dec=Angle(wdsTable['Coord (DEC) dms']), unit='hour, degree', frame="icrs")

print('WDS Catalog created successfully - ' , datetime.datetime.now())

sources_ds = Table()

# Set observation date and time
fitsFileDate = ''
fitsHeader = fits.open(workingDirectory + '/' + files[0])[0].header
key_to_lookup = 'DATE-OBS'
if key_to_lookup in fitsHeader:
    fitsFileDate = fits.open(workingDirectory + '/' + files[0])[0].header['DATE-OBS']
else:
    fitsFileDate = dummyObservationDate

### Run source detection, collect star data to Qtable
print('\n### Running source detection ###\n' , datetime.datetime.now())
for fitsFile in files:
    # 1. Read the list of sources extracted from an image (fits) file
    print('\n\n### Processing file: ', fitsFile, '###')
    fitsFileName = workingDirectory + '/' + fitsFile
    hdu = fits.open(fitsFileName)
    mywcs = WCS(hdu[0].header)

    # Estimate the background and background noise
    data = hdu[0].data
    mean, median, std = sigma_clipped_stats(data, sigma = dao_sigma)  

    daofind = DAOStarFinder(fwhm=dao_fwhm, threshold=dao_threshold*std)  
    sources = daofind(data - median)
    ra2, dec2 = mywcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 1)
    sources.add_column(ra2, name='ra_deg') 
    sources.add_column(dec2, name='dec_deg')

    sources_catalog = SkyCoord(ra=sources['ra_deg']*u.degree, dec=sources['dec_deg']*u.degree)
    idxw, idxs, wsd2d, wsd3d = search_around_sky(wds_catalog, sources_catalog, search_cone*u.deg)
    composit_catalog = hstack([wdsTable[idxw]['2000 Coord', 'Discov', 'Comp', 'PA_l', 'Sep_l', 'Mag_A', 'Mag_B'], sources[idxs]['id', 'mag', 'ra_deg', 'dec_deg']])
    companion_catalog = SkyCoord(ra=composit_catalog['ra_deg'] * u.degree, dec=composit_catalog['dec_deg'] * u.degree).directional_offset_by(composit_catalog['PA_l']*u.degree, composit_catalog['Sep_l']*u.arcsec)
    idxs2, d2ds2, d3ds2 = match_coordinates_sky(companion_catalog, sources_catalog)
    composit_catalog2 = hstack([composit_catalog, sources[idxs2]]) #['id', 'mag', 'ra_deg', 'dec_deg']

    sources_pa = SkyCoord(ra=composit_catalog2['ra_deg_1']*u.degree, dec=composit_catalog2['dec_deg_1']*u.degree).position_angle(SkyCoord(ra=composit_catalog2['ra_deg_2']*u.degree, dec=composit_catalog2['dec_deg_2']*u.degree)).to(u.deg)
    sources_sep = SkyCoord(ra=composit_catalog2['ra_deg_1']*u.degree, dec=composit_catalog2['dec_deg_1']*u.degree).separation(SkyCoord(ra=composit_catalog2['ra_deg_2']*u.degree, dec=composit_catalog2['dec_deg_2']*u.degree)).to(u.arcsec)
    sources_mag_diff = composit_catalog2['mag_2'] - composit_catalog2['mag_1']
    
    composit_catalog2.add_column(sources_pa, name='theta_measured')
    composit_catalog2.add_column(sources_sep, name='rho_measured')
    composit_catalog2.add_column(sources_mag_diff, name='mag_diff')
    sources_ds = vstack([sources_ds, composit_catalog2])

print('### Sources DS ###')
print(sources_ds)

upd_sources_ds = sources_ds[sources_ds['rho_measured'] != 0]
upd_sources_ds_by_object = upd_sources_ds.group_by(['2000 Coord', 'Discov', 'Comp'])

print(upd_sources_ds.info)
print('### Updated sources DS table grouped by WDS Identifier, Discoverer and Components ###')
print(upd_sources_ds_by_object)

upd_sources_ds_by_object.write(workingDirectory + '/double_stars.csv', format='ascii', overwrite=True, delimiter=',')


objectMean = upd_sources_ds_by_object.groups.aggregate(np.mean)

count = 1
for ds in upd_sources_ds_by_object.groups:
    print('\n#---------------------------------------------------------------------------------------------------------------------#')
    print('\n### Group index:', count, '###')
    print(ds)
    count = count + 1
    pairObjectId = ds[0]['2000 Coord'] + ds[0]['Discov'] + str(ds[0]['Comp'])
    
    # Search component in the Gaia DR3 database
    gaiaAStar, gaiaBStar = get_gaia_dr3_data(ds)

    if gaiaAStar and gaiaBStar:
        print('Gaia A star: ' + str(gaiaAStar[0]['DESIGNATION']))
        print('Gaia B star: ' + str(gaiaBStar[0]['DESIGNATION']))
        pairDistanceMinA = calcDistanceMin(float(gaiaAStar[0]['parallax']), float(gaiaAStar[0]['parallax_error']))
        pairDistanceMinB = calcDistanceMin(float(gaiaBStar[0]['parallax']), float(gaiaBStar[0]['parallax_error']))
        
        # Calculate physical binarity

        # Creating empty arrays for Star related calculations
        #Set input data
        starRa1 = float(gaiaAStar[0]['ra'])
        starDec1 = float(gaiaAStar[0]['dec'])
        starRa2 = float(gaiaBStar[0]['ra'])
        starDec2 = float(gaiaBStar[0]['dec'])
        starParallax1 = float(gaiaAStar[0]['parallax'])
        starParallaxError1 = float(gaiaAStar[0]['parallax_error'])

        # Value to modify Theta according to the appropriate quadrant

        '''addThetaValue = ()
        if deltaRa(starRa1, starRa2, starDec2) > 0 and deltaDec(starDec2, starDec1) > 0:
            addThetaValue = 0
        elif deltaRa(starRa1, starRa2, starDec2) > 0 and deltaDec(starDec2, starDec1) < 0:
            addThetaValue = 180
        elif deltaRa(starRa1, starRa2, starDec2) < 0 and deltaDec(starDec2, starDec1) < 0:
            addThetaValue = 180
        elif deltaRa(starRa1, starRa2, starDec2) < 0 and deltaDec(starDec2, starDec1) > 0:
            addThetaValue = 360
        elif deltaRa(starRa1, starRa2, starDec2) == 0 or deltaDec(starDec2, starDec1) == 0:
            addThetaValue = 0'''
                    
        # Calculate the widest possible separation for StarA
        possSep1 = possible_distance / calcDistanceMax(starParallax1, starParallaxError1)
        # rhoStar = rhoCalc(starRa1, starDec1, starRa2, starDec2)
        rhoStar = angular_separation(starRa1, starDec1, starRa2, starDec2)
        if possSep1 > rhoStar:
            starId1 = gaiaAStar[0]['solution_id']
            starName1 = gaiaAStar[0]['DESIGNATION']
            starId2 = gaiaBStar[0]['solution_id']
            starName2 = gaiaBStar[0]['DESIGNATION']
            starParallax2 = float(gaiaBStar[0]['parallax'])
            starParallaxError2 = float(gaiaBStar[0]['parallax_error'])
            starActualRa1 = float(ds['ra_deg_1'].mean())
            starActualDec1 = float(ds[0]['dec_deg_1'].mean())
            starActualRa2 = float(ds[0]['ra_deg_2'].mean())
            starActualDec2 = float(ds[0]['dec_deg_2'].mean())

            # Calculate actual data based on functions
            # thetaStar = thetaCalc(deltaRa(starRa1, starRa2, starDec2), deltaDec(starDec2, starDec1)) + addThetaValue
            # thetaActual = thetaCalc(deltaRa(starActualRa1, starActualRa2, starActualDec2), deltaDec(starActualDec2, starActualDec1)) + addThetaValue
            # rhoActual = rhoCalc(starActualRa1, starActualDec1, starActualRa2, starActualDec2)
            rhoActual = angular_separation(starActualRa1, starActualDec1, starActualRa2, starActualDec2)
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
        
        # Calculate attributes
        pairParallaxFactor, pairPmFactor, pairPmFactor, pairPmCommon, pairAbsMag1, pairAbsMag2, pairLum1, pairLum2, pairRad1, pairRad2, pairDR3Theta, pairDR3Rho, pairMass1, pairMass2, pairBVIndexA, pairBVIndexB, pairSepPar, pairEscapeVelocity, pairRelativeVelocity, pairHarshawFactor, pairHarshawPhysicality, pairBinarity = 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 

        pairParallaxFactor = (calcParallaxFactor(gaiaAStar[0]['parallax'], gaiaBStar[0]['parallax'])) * 100
        pairPmFactor = (calcPmFactor(gaiaAStar[0]['pmra'], gaiaAStar[0]['pmdec'], gaiaBStar[0]['pmra'], gaiaBStar[0]['pmdec']))
        pairPmCommon = calcPmCategory(pairPmFactor)
        pairAbsMag1 = calcAbsMag(gaiaAStar[0]['phot_g_mean_mag'], gaiaAStar[0]['parallax']) # Calculate Absolute magnitude
        pairAbsMag2 = calcAbsMag(gaiaBStar[0]['phot_g_mean_mag'], gaiaBStar[0]['parallax']) # Calculate Absolute magnitude
        pairLum1 = calcLuminosity(pairAbsMag1)
        pairLum2 = calcLuminosity(pairAbsMag2)
        pairAltLum1 = calcLuminosityAlternate(pairAbsMag1)
        pairAltLum2 = calcLuminosityAlternate(pairAbsMag2)
        pairRad1 = calcRadius(pairLum1, gaiaAStar[0]['teff_gspphot'])
        pairRad2 = calcRadius(pairLum2, gaiaBStar[0]['teff_gspphot'])
        # pairDR3Theta = thetaCalc(deltaRa(gaiaAStar[0]['ra'], gaiaBStar[0]['ra'], gaiaBStar[0]['dec']), deltaDec(gaiaBStar[0]['dec'], gaiaAStar[0]['dec'])) + addThetaValue
        pairDR3Theta = position_angle(gaiaAStar[0]['ra'], gaiaAStar[0]['dec'], gaiaBStar[0]['ra'], gaiaBStar[0]['dec'])
        # pairDR3Rho = rhoCalc(gaiaAStar[0]['ra'], gaiaAStar[0]['dec'], gaiaBStar[0]['ra'], gaiaBStar[0]['dec'])
        pairDR3Rho = angular_separation(gaiaAStar[0]['ra'], gaiaAStar[0]['dec'], gaiaBStar[0]['ra'], gaiaBStar[0]['dec'])
        pairMass1 = calcMass(pairLum1)
        pairMass2 = calcMass(pairLum2)
        pairBVIndexA = gaiaAStar[0]['phot_bp_mean_mag'] - gaiaAStar[0]['phot_g_mean_mag']
        pairBVIndexB = gaiaBStar[0]['phot_bp_mean_mag'] - gaiaBStar[0]['phot_g_mean_mag']
        pairSepPar = sepCalc(pairDistanceMinA, pairDistanceMinB, rhoStar) # Separation of the pairs in parsecs
        pairEscapeVelocity = calcEscapevelocity(pairMass1, pairMass2, pairSepPar, gravConst)
        pairRelativeVelocity = calcRelativeVelocity(gaiaAStar[0]['pmra'], gaiaAStar[0]['pmdec'], gaiaBStar[0]['pmra'], gaiaBStar[0]['pmdec'], gaiaAStar[0]['radial_velocity'], gaiaBStar[0]['radial_velocity'], pairDistanceMinA, pairDistanceMinB)
        pairHarshawFactor = calcHarshaw((pairParallaxFactor / 100), (pairPmFactor /100))
        pairHarshawPhysicality = calcHarshawPhysicality(pairHarshawFactor)
        pairBinarity = calcBinarity(pairRelativeVelocity, pairEscapeVelocity)
        
        # Calculate values for each pair based on the groups
        pairDesignationA = gaiaAStar[0]['DESIGNATION']
        pairDesignationB = gaiaBStar[0]['DESIGNATION']
        pairRaA = gaiaAStar[0]['ra']
        pairDecA = gaiaAStar[0]['dec']
        pairRaB = gaiaBStar[0]['ra']
        pairDecB = gaiaBStar[0]['dec']
        pairMagA = gaiaAStar[0]['phot_g_mean_mag']
        pairMagB = gaiaBStar[0]['phot_g_mean_mag']
        pairMeanTheta = float(ds['theta_measured'].degree.mean())
        pairMeanThetaErr = float(ds['theta_measured'].degree.std())
        pairMeanRho = float(ds['rho_measured'].arcsec.mean())
        pairMeanRhoErr = float(ds['rho_measured'].arcsec.std())
        pairMagnitudeA = ds[0]['mag_1']
        pairMagnitudeB = ds[0]['mag_2']
        pairGMagDiff = gaiaBStar[0]['phot_g_mean_mag'] - gaiaAStar[0]['phot_g_mean_mag']
        pairMagDiff = float((ds['mag_diff']).mean())
        pairMagDiffErr = (ds['mag_diff']).std()
        pairRadVelA = convertStringToNan(gaiaAStar[0]['radial_velocity'])
        pairRadVelErrA = convertStringToNan(gaiaAStar[0]['radial_velocity_error'])
        pairRadVelB = convertStringToNan(gaiaBStar[0]['radial_velocity'])
        pairRadVelErrB = convertStringToNan(gaiaBStar[0]['radial_velocity_error'])
        pairRadVelRatioA = convertStringToNan(math.fabs(gaiaAStar[0]['radial_velocity_error'] / gaiaAStar[0]['radial_velocity']) * 100)
        pairRadVelRatioB = convertStringToNan(math.fabs(gaiaBStar[0]['radial_velocity_error'] / gaiaBStar[0]['radial_velocity']) * 100)
        pairDesA = str(gaiaAStar[0]['DESIGNATION'])
        pairDesB = str(gaiaBStar[0]['DESIGNATION'])
        dateOfObservation = getUTC(fitsFileDate)
        #
        print('dateOfObservation: ', dateOfObservation)
        pairACurrentCoord = calcCurrentDR3Coord(dateOfObservation, pairRaA, pairDecA, gaiaAStar[0]['pmra'], gaiaAStar[0]['pmdec'])
        pairBCurrentCoord = calcCurrentDR3Coord(dateOfObservation, pairRaB, pairDecB, gaiaBStar[0]['pmra'], gaiaBStar[0]['pmdec'])
        pairAMeasuredCoord = SkyCoord(ra=ds['ra_deg_1'].groups.aggregate(np.mean) * u.deg, dec=ds['dec_deg_1'].groups.aggregate(np.mean) * u.deg)
        pairBMeasuredCoord = SkyCoord(ra=ds['ra_deg_2'].groups.aggregate(np.mean) * u.deg, dec=ds['dec_deg_2'].groups.aggregate(np.mean) * u.deg)
        pairACoordErr = pairACurrentCoord.separation(pairAMeasuredCoord)
        pairBCoordErr = pairBCurrentCoord.separation(pairBMeasuredCoord)
        #
        preciseCoord = str(getPreciseCoord(pairRaA, pairDecA, fitsFileDate))
        reportName = (workingDirectory + '/' + pairObjectId + '.txt')
        reportFile = open(reportName, "a")
        gaiaData = str(ds[0]['2000 Coord']) + ',' + str(ds[0]['Discov']) + ',' + str(gaiaAStar[0]['pmra']) + ',' + str(gaiaAStar[0]['pmdec']) + ',' + str(gaiaBStar[0]['pmra']) + ',' + str(gaiaBStar[0]['pmdec']) + ',' + str(gaiaAStar[0]['parallax']) + ',' + str(gaiaBStar[0]['parallax']) + ',' + str(calcDistance(gaiaAStar[0]['parallax'])) + ',' + str(calcDistance(gaiaBStar[0]['parallax'])) + ',' + str(gaiaAStar[0]['radial_velocity']) + ',' + str(gaiaBStar[0]['radial_velocity']) + ',' + 'pairRad1' + ',' + 'pairRad2' + ',' + str(pairLum1) + ',' + str(pairLum2) + ',' + str(gaiaAStar[0]['teff_gspphot']) + ',' + str(gaiaBStar[0]['teff_gspphot']) + ',' + str(gaiaAStar[0]['phot_g_mean_mag']) + ',' + str(gaiaBStar[0]['phot_g_mean_mag']) + ',' + str(gaiaAStar[0]['phot_bp_mean_mag']) + ',' + str(gaiaBStar[0]['phot_bp_mean_mag']) + ',' + str(gaiaAStar[0]['phot_rp_mean_mag']) + ',' + str(gaiaBStar[0]['phot_rp_mean_mag']) + ',' + str(pairDR3Theta) + ',' + str(pairDR3Rho) + ',' + str(gaiaAStar[0]['ra']) + ',' + str(gaiaAStar[0]['dec']) + ',' + str(gaiaBStar[0]['ra']) + ',' + str(gaiaBStar[0]['dec']) + ',' + str(gaiaAStar[0]['parallax_error']) + ',' + str(gaiaBStar[0]['parallax_error'])
        hrdPlot(pairObjectId, pairAbsMag1, pairAbsMag2, pairBVIndexA, pairBVIndexB)
        
    
    imagePlot(files[0], pairObjectId, pairRaA, pairDecA, pairRaB, pairDecB)
    
    # Print temp data
    print('\n### COMPONENTS ###')
    print('\nWDS Identifier:', ds[0]['2000 Coord'], ds[0]['Discov'], ds[0]['Comp'])
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
    print('\nTheta measurements\n') # , ds['dspaactual']
    print('Mean:', pairMeanTheta)
    print('Error:', pairMeanThetaErr)
    print('\nRho measurements\n') # , ds['dssepactual']
    print('Mean:', pairMeanRho)
    print('Error:', pairMeanRhoErr)
    print('\nMagnitude measurements\n') # , ds['dsmagdiff']
    print('Mean:', pairMagDiff)
    print('Error:', pairMagDiffErr)
    print('\n\nParallax factor:', pairParallaxFactor, '%')
    print('Proper motion factor:', pairPmFactor * 100, '%')
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
    print('Radial velocity of the stars', 'A:', pairRadVelA, 'km/s (Err:', pairRadVelErrA, 'km/s)', 'B:', pairRadVelB, 'km/s (Err:', pairRadVelErrB, 'km/s)')
    print('Radial velocity ratio A:', pairRadVelRatioA, '%')
    print('Radial velocity ratio B:', pairRadVelRatioB, '%')
    print('Separation:', pairSepPar, 'parsec,', pairSepPar * 206265, 'AU')
    print('Pair Escape velocity:', pairEscapeVelocity, 'km/s')
    print('Pair Relative velocity:', pairRelativeVelocity, 'km/s')
    print('Pair Harshaw factor:', pairHarshawFactor)
    print('Pair Harshaw physicality:', pairHarshawPhysicality)
    print('Pair binarity:', pairBinarity)
    print('Analysis finished: ' , datetime.datetime.now())
    
    # Write results to file
    reportTable.add_row([ds[0]['2000 Coord'] + ds[0]['Discov'] + str(ds[0]['Comp']), dateOfObservation, pairMeanTheta, pairMeanThetaErr, pairMeanRho, pairMeanRhoErr, np.nan, np.nan, pairMagDiff, pairMagDiffErr, 'Filter wawelenght', 'filter FWHM', '0.2', '1', 'TAL_2022', 'C', '7', preciseCoord])
    reportFile.write('### WDS Data ###')
    reportFile.write('\nWDS Identifier: ' + ds[0]['2000 Coord'])
    reportFile.write('\nDiscoverer and components: ' + str(ds[0]['Discov']) + ' ' + str(ds[0]['Comp']))
    reportFile.write('\nMagnitude Pri: ' + str(ds[0]['Mag_A']))
    reportFile.write('\nMagnitude Sec: ' + str(ds[0]['Mag_B']))
    reportFile.write('\nPA last: ' + str(ds[0]['PA_l']))
    reportFile.write('\nSep last: ' +  str(ds[0]['Sep_l']))
    reportFile.write('\n\n### Gaia DR3 Data ###')
    reportFile.write('\nMain star: ' + pairDesA)
    reportFile.write('\nCompanion: ' + pairDesB)
    reportFile.write('\nPair G magnitudes A: ' + str(roundNumber(pairMagA)) + ' B: ' + str(roundNumber(pairMagB)))
    reportFile.write('\nPosition angle: ' + str(roundNumber(pairDR3Theta)))
    reportFile.write('\nSeparation: ' + str(roundNumber(pairDR3Rho)))
    reportFile.write('\nMagnitude difference: ' + str(roundNumber(pairGMagDiff)))
    reportFile.write('\nPrecise coordinates (J2000): ' + preciseCoord)
    reportFile.write('\nDate of observation: ' + dateOfObservation)
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
    reportFile.write('\n\n### Measurements ###')
    reportFile.write('\nPosition angle:')
    reportFile.write('\nTheta measurements' + str(ds['theta_measured'].degree))
    reportFile.write('\nMean: ' + str(roundNumber(pairMeanTheta)))
    reportFile.write('\nError: ' + str(roundNumber(pairMeanThetaErr)))
    reportFile.write('\n\nSeparation:')
    reportFile.write('\nRho measurements\n' + str(ds['rho_measured'].arcsec))
    reportFile.write('\nMean: ' + str(roundNumber(pairMeanRho)))
    reportFile.write('\nError: ' + str(roundNumber(pairMeanRhoErr)))
    reportFile.write('\n\nMagnitude measurements\n'  + str(ds['mag_diff']))
    reportFile.write('\nMean: ' + str(roundNumber(pairMagDiff)))
    reportFile.write('\nError: ' + str(roundNumber(pairMagDiffErr)))
    reportFile.write('\n\n### Calculated attributes ###')
    reportFile.write('\nSeparation (Measured): ' + str(roundNumber(pairMeanRho)))
    reportFile.write('\nPosition angle (Measured): ' + str(roundNumber(pairMeanTheta)))
    reportFile.write('\nMagnitude difference (Measured): ' + str(roundNumber(pairMagDiff)) + ' (Err: ' + str(roundNumber(pairMagDiffErr)) + ')')
    reportFile.write('\nAbsolute magnitude A: ' + str(roundNumber(pairAbsMag1)))
    reportFile.write('\nAbsolute magnitude B: ' + str(roundNumber(pairAbsMag2)))
    reportFile.write('\nLuminosity A: ' + str(roundNumber(pairLum1)))
    reportFile.write('\nLuminosity B: ' + str(roundNumber(pairLum2)))
    reportFile.write('\nRad A: ' + str(roundNumber(pairRad1)))
    reportFile.write('\nRad B: ' + str(roundNumber(pairRad2)))
    reportFile.write('\nMass A: ' + str(roundNumber(pairMass1)))
    reportFile.write('\nMass B: ' + str(roundNumber(pairMass2)))
    reportFile.write('\nBV index A: ' + str(roundNumber(pairBVIndexA)) + ' B: ' + str(roundNumber(pairBVIndexB)))
    reportFile.write('\nRadial velocity of the stars ' + 'A:' + str(roundNumber(pairRadVelA)) + 'km/s (Err:' + str(roundNumber(pairRadVelErrA)) + 'km/s)' + ' B:' + str(roundNumber(pairRadVelB)) + 'km/s (Err:' + str(roundNumber(pairRadVelErrB)) + 'km/s)')
    reportFile.write('\nRadial velocity ratio A: ' + str(roundNumber(pairRadVelRatioA)) + ' %')
    reportFile.write('\nRadial velocity ratio B: ' + str(roundNumber(pairRadVelRatioB)) + ' %')
    reportFile.write('\nSeparation: ' + str(roundNumber(pairSepPar)) + ' parsec, ' + str(roundNumber((pairSepPar * 206265))) + ' AU')
    reportFile.write('\nPair Escape velocity: ' + str(roundNumber(pairEscapeVelocity)) + ' km/s')
    reportFile.write('\nPair Relative velocity: ' + str(roundNumber(pairRelativeVelocity)) + ' km/s')
    reportFile.write('\n\n### Analysis ###')
    reportFile.write('\nParallax factor: ' + str(roundNumber(pairParallaxFactor)) + ' %')
    reportFile.write('\nProper motion factor: ' + str(roundNumber(pairPmFactor) * 100) + ' %')
    reportFile.write('\nProper motion category: '+ str(pairPmCommon))
    reportFile.write('\nPair Harshaw factor: ' + str(roundNumber(pairHarshawFactor)))
    reportFile.write('\nPair Harshaw physicality: ' + str(pairHarshawPhysicality))
    reportFile.write('\nPair binarity: ' + str(pairBinarity))
    reportFile.write('\n\n### WDS form:\n')
    wdsform = str(ds[0]['2000 Coord']) + ',' + dateOfObservation + ',' +  str(roundNumber(pairMeanTheta)) + ',' +  str(roundNumber(pairMeanThetaErr)) + ',' +  str(roundNumber(pairMeanRho)) + ',' +  str(roundNumber(pairMeanRhoErr)) + ',' +  'nan' + ',' +  'nan' + ',' +  str(roundNumber(pairMagDiff)) + ',' +  str(roundNumber(pairMagDiffErr)) + ',' + 'Filter wawelenght' + ',' + 'filter FWHM' + ',' + '0.2' + ',' + '1' + ',' + 'TLB_2023' + ',' +  'C' + ',' + '7'+ ',' + str(getPreciseCoord(pairRaA, pairDecA, fitsFileDate))
    reportFile.write(str(wdsform))
    reportFile.write('\n\n### Gaia data:\n')
    reportFile.write(str(gaiaData))

reportTable.write(workingDirectory + '/double_stars_wds_format.txt', format='ascii', overwrite=True, delimiter=',')