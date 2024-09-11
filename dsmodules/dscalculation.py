# This file contains functions needed for the neccessary calcualtions of double star image processing.

import pprint

import os
import sys
import numpy as np
import numpy.ma as ma
import math
import sys
import datetime
from time import gmtime, strftime
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy.coordinates import match_coordinates_sky, search_around_sky, position_angle, angular_separation
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack, hstack
from astropy.table import Column, MaskedColumn
from astropy.utils.masked import MaskedNDArray
from photutils.detection import DAOStarFinder
from astropy.wcs import WCS
import astropy.units as u
import warnings
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta
warnings.filterwarnings("ignore")
from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default

# configurations
possible_distance = 30000.0 # AU
search_cone = 0.001 # Decimal degree
dummyObservationDate = "2022-01-01T12:00:00"
gaia_dr3_epoch = 2016.0

# Gravitational constant is convenient if measure distances in parsecs (pc), velocities in kilometres per second (km/s) and masses in solar units M
gravConst = 0.0043009 
auToParsec = 0.0000048481368111358
image_limit = 10000

# Constant to calculate star luminosity and mass
sun_luminosity = 3.0128 * (10 ** 28)
sun_absolute_luminosity = 3.828 * (10 ** 26)

# Constant variables
hipparcos_file = Table.read(f"C:/Users/gerge/Documents/Catalogs/Hipparcos/I_239_selection.csv", format='ascii')

# Insert the downloaded wds file path here
wds_file = "C:/Users/gerge/Documents/Catalogs/WDS/wdsweb_summ2.txt"



# function to test the package
def print_hello_world():
    print('Hello World!')

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

# Function to calculate Rho (separation)
def rhoCalc(raa, deca, rab, decb):
    rhocalc = math.sqrt(((raa-rab) * math.cos(math.radians(deca))) ** 2 + (deca - decb) ** 2) * 3600
    return rhocalc

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
    if (deltara != 0 and deltadec != 0):
        thetacalc = math.degrees(math.atan(deltara/deltadec))
    else:
        thetacalc = 0
    return thetacalc

# Function to calculate the separation of the two stars in parsecs
# Excel formula =IF('min distance A'>'min distance b','min distance A'*'Rho','min distance b'*'Rho')
# Excel formula (error percent): K7, L7=Parallax error/Parallax
# Excel formula: =(((1-(parallax error a / parallax a))*K5)+((1-L7)*L5))/((1-K7)+(1-L7))
def sepCalc(dist_a, dist_b, rho):
    if dist_a > dist_b:
        sep = (dist_b * rho) * auToParsec
    else:
        sep = (dist_a * rho) * auToParsec
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
    # =HA(VAGY(A plx<0;B plx<0);"N/R";HA(VAGY(Dist A=0;Dist B=0);0;ABS(HA(ÉS(A plx=0;B plx=0);0;HA((Dist A-Dist B)=0;0.75;((1-ABS((Dist A-Dist B)/((Dist A+Dist B)/2)))*0.75))))))
    parfac = 1 - math.fabs((para-parb)/(0.5*(para+parb)))
    return parfac

# Function to calculate the proper motion factor of the stars and define if it is a CPM (common proper motion) pair
# EXCEL formula =ABS(1-((SQRT(pmraA-pmraB)^2+(pmdecA-pmdecB)^2)/(SQRT(pmraA^2+pmdecA^2)+(pmraB^2+pmdecB^2)))))
def calcPmFactor(pmraa, pmdeca, pmrab, pmdecb):
    pmfac = math.fabs(1-((math.sqrt(((pmraa-pmrab) ** 2) + ((pmdeca-pmdecb) ** 2))/(math.sqrt((pmraa ** 2) + (pmdeca ** 2))+((pmrab ** 2) + pmdecb ** 2)))))
    ### 2 Átlag saját mozgás vektor: AVG PM prob
    ### =(GYÖK(pmra A^2+pm dec A^2)+GYÖK(pm ra B^2+pm dec B^2))/2
    ### 3 sajátmozgás valószínűség: PM prob
    ### =HA(Vr A=0;0;ABS(1-(GYÖK((pm ra A-pm ra B)^2+(pm dec A-pm dec B)^2)/2 pont eredménye))*0.15)
    ### distance = math.sqrt((pmraa - pmrab)**2 + (pmdeca - pmdecb)**2)
    ### pmfac = abs(1 - (distance / 2)) * 0.15
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
        radius = math.sqrt(luminosity / ((float(teffStar) / teffSun) ** 4))
    else:
        radius = 'Cannot be calculated, T(eff) or Luminosity is missing'
    return radius

# Function to calculate Harshaw probapility of duplicity based on the parallax and proper motion factors
def calcHarshaw(parallaxFactor, pmFactor):
    ### =HA(VAGY(A plx<0;B plx <0);1 pont eredménye;SZUM(3 pont eredménye+1 pont eredménye))
    HarshawFactor = (parallaxFactor * 0.75) + (pmFactor * 0.15)
    return HarshawFactor

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
    radspeeddif = math.fabs(float(radvela) - float(radvelb))
    sumspeeddiff = math.sqrt(tanspeeddiff ** 2 + radspeeddif ** 2)
    return sumspeeddiff


# Function to calculate the Escape velocity of the system, separation should be calculated in parsec!

def calcEscapevelocity(mass_a, mass_b, separation, gravconst):
    if bool(mass_a) and bool(mass_b) and bool(separation) and bool(gravconst):
        # print('calcEscapevelocity: ' + str(mass_a) + ' ' + str(mass_b) + ' ' + str(separation) + ' ' + str(gravconst))
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
    if pmfact >= 80.0:
        pmCommon = 'CPM'
    elif 40.0 <= pmfact < 80.0:
        pmCommon = 'SPM'
    elif pmfact < 40.0:
        pmCommon = 'DPM'
    return pmCommon

# Function to search coordinates in the WDS catalog file
def search_in_wds(ra_source, dec_source, wdsTable, wds_catalog):
    coordinates = SkyCoord(ra=ra_source, dec=dec_source)
    idx, d2d, d3d = coordinates.match_to_catalog_sky(wds_catalog)
    star_coord = SkyCoord(ra=ra_source, dec=dec_source)
    print('Coordinates: ' + str(ra_source) + '(Ra), ' + str(dec_source) + '(Dec), Separation: ' + str(d2d))
    print(wdsTable[idx]['2000 Coord'])
    print(wdsTable[np.where(wds_catalog['2000 Coord'] == wdsTable[idx]['2000 Coord'])])

def create_unique_id(wds_id, discov):
    id_array = [i + j for i, j in zip(wds_id, discov)]
    return id_array

def delete_invalid_lines_wds(catalog):
    rows_to_delete_ra = np.where(catalog['Coord (RA) hms'] == '.')
    catalog.remove_rows(rows_to_delete_ra)
    rows_to_delete_dec = np.where(catalog['Coord (DEC) dms'] == '.')
    catalog.remove_rows(rows_to_delete_dec)
    return catalog

def calculate_wds_ra_hourangle(wds_ra_array):
    wds_ra_hms = []
    for star in wds_ra_array:
        if len(str(star)) == 9:
            wds_ra_hms.append(str(star[0:2]) + 'h' + str(star[2:4]) + 'm' + str(star[4:9]) + 's')
        else:
            wds_ra_hms.append('.')
    return wds_ra_hms

def calculate_wds_dec_hourangle(wds_dec_array):
    wds_dec_dms = []
    for star in wds_dec_array:
        if len(str(star)) == 9:
            wds_dec_dms.append(str(star[0:3]) + 'd' + str(star[3:5]) + 'm' + str(star[5:9]) + 's')
        else:
            wds_dec_dms.append('.')
    return wds_dec_dms

# Create HRD plot of the double stars based on Hipparcos
def hrdPlot(pairname, working_directory, mag_abs_a, mag_abs_b, bv_a, bv_b):
    print(pairname, mag_abs_a, mag_abs_b, bv_a, bv_b)
    if pairname and mag_abs_a and mag_abs_b and bv_a and bv_b:
        hipparcos_abs_mag = hipparcos_file['Abs_mag']
        hipparcos_bv_index = hipparcos_file['B-V']
        colors = (hipparcos_bv_index)
        plt.scatter(hipparcos_bv_index, hipparcos_abs_mag, c=colors, s=0.5, alpha=0.1, cmap='RdYlBu_r', vmax=1.9, vmin=-0.4) #, 
        plt.scatter(bv_a, mag_abs_a, s=14, color="blue", label='Main star') # s= 1 / mag_abs_a
        plt.scatter(bv_b, mag_abs_b, s=7, color="red", label='Companion star') # s= 1 / mag_abs_a
        plt.legend(loc="upper left")
        plt.axis((-0.4,1.9,15,-10))
        plt.title('Double Star ' + pairname + ' H-R Diagram')
        plt.xlabel('B-V index')
        plt.ylabel('Absolute magnitude')
        plt.gca().set_aspect(0.1)
        #savename = str(working_directory + '/' + pairname + '_hrd.jpg').replace(' ', '')
        savename = str(working_directory + '/' + pairname + '_hrd.jpg').replace(' ', '')
        plt.savefig(savename, bbox_inches='tight', dpi=300.0)
        plt.close()
    else:
        print('Data is missiong, HRD plot cannot be created!')


# Create Image plot of the double stars
def imagePlot(filename, working_directory, pairname, raa, deca, rab, decb):
    # data[0]
    #image_data = fits.open(working_directory + '/' + filename)
    image_data = fits.open(filename)
    print('IMAGE DATA:', working_directory, '/', filename)
    header = image_data[0].header
    wcs_helix = WCS(image_data[0].header, naxis=2)
    image = image_data[0].data
    image_height = header['NAXIS2']
    star_a = SkyCoord(raa * u.deg, deca * u.deg, frame='icrs')
    star_b = SkyCoord(rab * u.deg, decb * u.deg, frame='icrs')
    star_a_pix = utils.skycoord_to_pixel(star_a, wcs_helix)
    star_b_pix = utils.skycoord_to_pixel(star_b, wcs_helix)
    plt.scatter(star_a_pix[0] + 40, star_a_pix[1], marker="_", s=50, color="grey")
    plt.scatter(star_a_pix[0], star_a_pix[1] + 40, marker="|", s=50, color="grey")
    plt.scatter(star_b_pix[0] + 40, star_b_pix[1], marker="_", s=50, color="grey")
    plt.scatter(star_b_pix[0], star_b_pix[1] + 40, marker="|", s=50, color="grey")
    plt.title(pairname)
    plt.imshow(image, origin='lower',cmap='Greys', aspect='equal', vmax=image_limit, vmin=0) # , cmap='cividis'
    plt.savefig(str(working_directory + '/' + pairname + '_img.jpg').replace(' ', ''),dpi=300.0, bbox_inches='tight', pad_inches=0.2)
    plt.close()

def get_gaia_dr3_data(doublestars):
    pairACoord = SkyCoord(ra=doublestars['ra_deg_1'].mean(), dec=doublestars['dec_deg_1'].mean(), unit=(u.degree, u.degree), frame='icrs')
    pairBCoord = SkyCoord(ra=doublestars['ra_deg_2'].mean(), dec=doublestars['dec_deg_2'].mean(), unit=(u.degree, u.degree), frame='icrs')
    a = Gaia.cone_search_async(pairACoord, radius=u.Quantity(search_cone, u.deg))
    b = Gaia.cone_search_async(pairBCoord, radius=u.Quantity(search_cone, u.deg))
    a_query = a.get_results()
    b_query = b.get_results()
    return a_query, b_query

def calc_average_distance(star_a_par, star_a_par_err, star_b_par, star_b_par_err, sep):
    # Calcualte average distance in lightyears
    # Excel formula: = 1000 / (((1-star A parallax) * star A parallax error in %)+((1-star B parallax) * star B parallax error in %))/((1-star A parallax error in %)+(1-star B parallax error in %))
    err_a_per = star_a_par_err / star_a_par
    err_b_per = star_b_par_err / star_b_par
    wtd_px = (((1-err_a_per)*star_a_par)+((1-err_b_per)*star_b_par))/((1-err_a_per)+(1-err_b_per))
    wtd_dist = 1000 / wtd_px
    wtd_sep = wtd_dist * sep
    return wtd_sep, wtd_dist, wtd_sep
    


# Calculate maximum orbit speed
def calc_historic_orbit(massa, massb, sep_pc, avg_distance, dr3_rho, wds_first_rho, measured_rho, wds_first_theta, measured_theta, wds_first_obs, obs_time):
    # Check data types
    print('\n ##### Calc history orbit #####\n')
    print('massa: ', massa, ' (', type(massa), ')')
    print('massb: ', massb, ' (', type(massb), ')')
    print('sep_pc: ', sep_pc, ' (', type(sep_pc), ')')
    print('avg_distance: ', avg_distance, ' (', type(avg_distance), ')')
    print('dr3_rho: ', dr3_rho, ' (', type(dr3_rho), ')')
    print('wds_first_rho: ', wds_first_rho, ' (', type(wds_first_rho), ')')
    print('measured_rho: ', measured_rho, ' (', type(measured_rho), ')')
    print('wds_first_theta: ', wds_first_theta, ' (', type(wds_first_theta), ')')
    print('measured_theta: ', measured_theta, ' (', type(measured_theta), ')')
    print('wds_first_obs: ', wds_first_obs, ' (', type(wds_first_obs), ')')
    print('obs_time: ', obs_time, ' (', type(obs_time), ')')

    # Calculate historical delta position angle (theta in decimal degree)
    delta_theta = math.fabs(wds_first_theta - measured_theta)
    
    # Calculate historical delta separation (rho in arcsconds)
    delta_rho = math.fabs(wds_first_rho - measured_rho)
    
    # Calculate historical delta time (years)
    delta_time = float(obs_time) - float(wds_first_obs)
    
    # Calculate half axis
    half_axis = avg_distance * (1.26 * dr3_rho)
    print('half_axis: ', half_axis, ' (', type(half_axis), ')')
    
    # Calculate maximum orbital velocity
    # GYÖK(0.0043*($M$20+$N$20)*(2/($N$11*0.00000485)-1/(N25*0.00000485)))
    max_orbit_velolicy = math.sqrt(gravConst * ((massa + massb) * (2 / (sep_pc)) - 1 / (half_axis * 0.00000485)))
    #max_orbit_velolicy_alt = math.sqrt(0.0043 * (massa + massb) * (2 / (sep_pc) - 1 / (half_axis * 0.00000485)))
    
    # GYÖK((P33*(P35/$I$11))^2+(P34/$I$11)^2)
    relative_velocity = math.sqrt((measured_rho * (delta_theta / delta_time)) ** 2 + (delta_rho / delta_time) ** 2)
    #relative_velocity_alt = math.sqrt((measured_rho * (delta_theta / delta_time)) ** 2 + (delta_rho / delta_time) ** 2)
    
    # 0.0474*$N$10*P36
    observed_velocity = 0.0474 * avg_distance * relative_velocity
    historic_criterion = ''
    input_data_variables = {
        "massa" : massa,
        "massb" : massb,
        "sep_pc" : sep_pc, 
        "avg_distance" : avg_distance, 
        "dr3_rho" : dr3_rho, 
        "wds_first_rho" : wds_first_rho, 
        "measured_rho" : measured_rho, 
        "wds_first_theta" : wds_first_theta, 
        "measured_theta" : measured_theta, 
        "wds_first_obs" : wds_first_obs, 
        "obs_time" : obs_time
    }
    
    if observed_velocity < max_orbit_velolicy:
        historic_criterion = 'Physical'
    else:
        historic_criterion = 'Optical'   
    return historic_criterion, max_orbit_velolicy, observed_velocity, input_data_variables, delta_theta, delta_rho, delta_time

def calculate_photo_center(wcs, header):
    photo_left_upper = SkyCoord.from_pixel(0, 0, wcs, origin=0, mode='all')
    photo_right_lower = SkyCoord.from_pixel(header['NAXIS2'], header['NAXIS1'], wcs, origin=0, mode='all')
    center = SkyCoord(header['CRVAL1'] * u.degree, header['CRVAL2'] * u.degree)
    radius = photo_left_upper.separation(photo_right_lower) / 2
    print('Center of photo: ', center.to_string('hmsdms'), '/', center.to_string('decimal'),
      '\nRadius of photo: ', radius)
    return center, radius

def get_objects_from_catalog(catalog, photo_center, photo_radius):
    d2d = photo_center.separation(catalog)
    catalog_mask = d2d < photo_radius
    return catalog_mask

def plot_image_with_frame(image_array, frame_edges, output_filename): #
    print('frame edges: ', frame_edges)
    # Create a figure and axis to plot the image
    fig, ax = plt.subplots()
    
    # Display the image
    ax.imshow(image_array, cmap='gray', vmax=image_limit, vmin=0)
    
    # Get the frame edges from the dictionary
    left = frame_edges[0][0]
    right = frame_edges[0][1]
    top = frame_edges[1][0]
    bottom = frame_edges[1][1]
    
    # Plot the frame as a rectangle
    rect = plt.Rectangle((left, top), right-left, bottom-top, edgecolor='red', facecolor='none', linewidth=2)
    ax.add_patch(rect)

    # Plot small circles at specified coordinates
    '''for coord in circle_coords:
        circle = plt.Circle(coord, radius=5, color='blue', fill=True)
        ax.add_patch(circle)'''
    
    # Adjust the axis limits to match the image size
    ax.set_xlim(0, image_array.shape[1])
    ax.set_ylim(image_array.shape[0], 0)
    
    # Add axis labels
    ax.set_xlabel('X Coordinate (pixels)')
    ax.set_ylabel('Y Coordinate (pixels)')
    
    # Show x and y scales with ticks
    ax.set_xticks(np.arange(0, image_array.shape[1], step=image_array.shape[1] // 10))
    ax.set_yticks(np.arange(0, image_array.shape[0], step=image_array.shape[0] // 10))
        
    # Save the result as a JPG file
    plt.savefig(str(workingDirectory + '/' + output_filename), format='jpg', bbox_inches='tight', pad_inches=0)
    
    # Close the plot to free up memory
    plt.close()

def define_image_plane(wcs, header):
    # Define inner rectangle coordinates
    # Define the inner rectangle percent in image pixels measured from the edge of the image

    frame_sink = (frame_percentage / 100) * header['NAXIS2']

    photo_x_coords = [frame_sink, header['NAXIS1'] - frame_sink]
    photo_y_coords = [frame_sink, header['NAXIS2'] - frame_sink]


    photo_center = SkyCoord(header['CRVAL1'] * u.degree, header['CRVAL2'] * u.degree)
    photo_left_upper = SkyCoord.from_pixel(photo_y_coords[0], photo_x_coords[0], wcs, origin=0, mode='all')
    photo_left_lower = SkyCoord.from_pixel(photo_y_coords[1], photo_x_coords[0], wcs, origin=0, mode='all')
    photo_right_upper = SkyCoord.from_pixel(photo_y_coords[0], photo_x_coords[1], wcs, origin=0, mode='all')
    photo_right_lower = SkyCoord.from_pixel(photo_y_coords[1], photo_x_coords[1], wcs, origin=0, mode='all')
    # print('photo_coords in pixel: ', photo_x_coords, photo_y_coords)
    # print('photo_coords in degree: ', SkyCoord.to_pixel(photo_left_upper, wcs, origin=0, mode='all'), SkyCoord.to_pixel(photo_left_lower, wcs, origin=0, mode='all'), SkyCoord.to_pixel(photo_right_upper, wcs, origin=0, mode='all'))
    half_width = photo_left_upper.separation(photo_right_upper) / 2
    half_height = photo_left_upper.separation(photo_left_lower) / 2
    corners = SkyCoord(
        ra=[photo_center.ra - half_width, photo_center.ra + half_width,
            photo_center.ra + half_width, photo_center.ra - half_width],
        dec=[photo_center.dec - half_height, photo_center.dec - half_height,
            photo_center.dec + half_height, photo_center.dec + half_height]
        )
    frame = [photo_x_coords, photo_y_coords]

    return corners, frame


def catalog_search_in_image(wcs, header, center, radius, wds_catalog_list):
    #image_region = define_image_region(wcs, header)
    d2d = center.separation(wds_catalog_list)
    catalogmsk = d2d < radius
    idxcatalog = np.where(catalogmsk)[0]
    image_plane = define_image_plane(wcs, header)[0]

    in_fov = []
    for obj in wds_catalog_list[idxcatalog]:
        if (np.min(image_plane.ra) <= obj.ra <= np.max(image_plane.ra)) and (np.min(image_plane.dec) <= obj.dec <= np.max(image_plane.dec)):
            in_fov.append(wds_catalog_list)

    return wdsTable[idxcatalog]

def exit_if_no_doubles_found(doubles_table):
    if len(doubles_table) == 0:
        print("Double stars not found on images. Exiting the program.")
        sys.exit()

## Add general functions of ds measurements

# Set observation date and time
def set_observation_date(file_header):
    fitsFileDate = ''
    key_to_lookup_a = 'DATE-OBS'
    key_to_lookup_b = 'DATE'
    if key_to_lookup_a in file_header:
        fitsFileDate = file_header['DATE-OBS']
    elif key_to_lookup_b in file_header:
        fitsFileDate = file_header['DATE']
    else:
        fitsFileDate = np.nan
    
    return fitsFileDate

def get_sources_from_image(image_data, image_wcs, dao_sigma, dao_fwhm, dao_threshold, fits_file_name, fits_file_date):
    mean, median, std = sigma_clipped_stats(image_data, sigma = dao_sigma)  

    daofind = DAOStarFinder(fwhm=dao_fwhm, threshold=dao_threshold*std)  
    sources = daofind(image_data - median)
    ra2, dec2 = image_wcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 1)
    sources.add_column(ra2, name='ra_deg') 
    sources.add_column(dec2, name='dec_deg')
    sources.add_column(fits_file_date, name='image_date')
    sources.add_column(fits_file_name, name='file')

    return sources

def get_gaia_dr3_data_offline(doublestars, segment_lib):
    pairACoord = SkyCoord(ra=doublestars['ra_deg_1'].mean(), dec=doublestars['dec_deg_1'].mean(), unit=(u.degree, u.degree), frame='icrs')
    pairBCoord = SkyCoord(ra=doublestars['ra_deg_2'].mean(), dec=doublestars['dec_deg_2'].mean(), unit=(u.degree, u.degree), frame='icrs')
    ds_catalog = SkyCoord(ra=doublestars['ra_deg_1'], dec=doublestars['dec_deg_1'], unit=(u.degree, u.degree), frame='icrs')
    print(ds_catalog.ra.deg, ds_catalog.dec.deg)

    segments = []
    for star in ds_catalog:
        #ra, dec = wcs.all_pix2world([[star ['xcentroid'], star ['ycentroid']]], 0)[0]
        ra, dec = star.ra.deg, star.dec.deg
        segmentRaCalc = int((float(ra) // 5) + 1)
        segmentDecCalc = int((float(dec) // 5) + 1)
        segmentName = f"{segmentRaCalc}-{segmentDecCalc}.csv"
        if segmentName not in segments:
            segments.append(segmentName)
    print('### Segments are:',segments)

    # Read all segments into an array
    gaiaStars = np.empty((0, 20))

    # Add all segments to the numpy array
    for seg in segments:
        if len(segments) > 1:
            #segmentpart = Table.read(f"C:\Astro\catalogs\GaiaDR3\gaiadr3_15mag_catalog\{seg}", format='ascii')
            segmentpart = Table.read(segment_lib + str(seg), format='ascii')
            gaiaStars = vstack([gaiaStars, segmentpart])
        else:
            #segmentpart = Table.read(f"C:\Astro\catalogs\GaiaDR3\gaiadr3_15mag_catalog\{seg}", format='ascii')
            segmentpart = Table.read(segment_lib + str(seg), format='ascii')
            gaiaStars = segmentpart

    gaia_catalog = SkyCoord(ra=gaiaStars['ra'], dec=gaiaStars['dec'], unit='degree, degree', frame="icrs")
    idx_a, d2d_a, d3d_a = match_coordinates_sky(pairACoord, gaia_catalog) #, 3*u.arcseconds)
    idx_b, d2d_b, d3d_b = match_coordinates_sky(pairBCoord, gaia_catalog) #, , 3*u.arcseconds

    return gaiaStars[idx_a], gaiaStars[idx_b]

def create_wds_table(wdsdata):
    wdsTable = hstack([wdsdata, calculate_wds_ra_hourangle(wdsdata['Coord (RA)'])])
    wdsTable.rename_column('col0', 'Coord (RA) hms')
    wdsTable = hstack([wdsTable, calculate_wds_dec_hourangle(wdsdata['Coord (DEC)'])])
    wdsTable.rename_column('col0', 'Coord (DEC) dms')
    wdsTable = hstack([wdsTable, create_unique_id(wdsdata['2000 Coord'], wdsdata['Discov'])])
    wdsTable.rename_column('col0', 'Unique ID')
    wdsTable = delete_invalid_lines_wds(wdsTable)

    return wdsTable

def get_fits_data(fits_file):
    hdu = fits.open(fits_file)
    wcs = WCS(hdu[0].header, naxis=2)
    fits_header = hdu[0].header
    fits_data = hdu[0].data
    
	# Set observation date and time
    fits_file_date = ''
    key_to_lookup_a = 'DATE-OBS'
    key_to_lookup_b = 'DATE'
    if key_to_lookup_a in fits_header:
        fits_file_date = fits_header['DATE-OBS']
    elif key_to_lookup_b in fits_header:
        fits_file_date = fits_header['DATE']
    else:
        fits_file_date = np.nan

    return hdu, wcs, fits_header, fits_data, fits_file_date

def get_sources_from_image(sources_ds, wds_catalog, fits_data, fits_header, fits_file_name, fits_file_date, image_wcs, wds_table, dao_sigma, dao_fwhm, dao_threshold):
    ds_sources = Table()
    mean, median, std = sigma_clipped_stats(fits_data, sigma = dao_sigma)  
    daofind = DAOStarFinder(fwhm=dao_fwhm, threshold=dao_threshold*std)  
    sources = daofind(fits_data - median)
    ra2, dec2 = image_wcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 1)
    sources.add_column(ra2, name='ra_deg') 
    sources.add_column(dec2, name='dec_deg')
    sources.add_column(fits_file_date, name='image_date')
    sources.add_column(fits_file_name, name='file')
    photo_center, photo_radius = calculate_photo_center(image_wcs, fits_header)
    sources_catalog = SkyCoord(ra=sources['ra_deg']*u.degree, dec=sources['dec_deg']*u.degree, frame='fk5')
    idxw, idxs, wsd2d, wsd3d = search_around_sky(wds_catalog, sources_catalog, search_cone*u.deg)
    composit_catalog = hstack([wds_table[idxw]['2000 Coord', 'Discov', 'Comp', 'Date (first)', 'PA_f', 'PA_l', 'Sep_f', 'Sep_l', 'Mag_A', 'Mag_B'], sources[idxs]['id', 'mag', 'ra_deg', 'dec_deg']])
    companion_catalog = SkyCoord(ra=composit_catalog['ra_deg'] * u.degree, dec=composit_catalog['dec_deg'] * u.degree).directional_offset_by(composit_catalog['PA_l']*u.degree, composit_catalog['Sep_l']*u.arcsec)
    idxs2, d2ds2, d3ds2 = match_coordinates_sky(companion_catalog, sources_catalog)
    composit_catalog2 = hstack([composit_catalog, sources[idxs2]]) #['id', 'mag', 'ra_deg', 'dec_deg']

    sources_pa = SkyCoord(ra=composit_catalog2['ra_deg_1']*u.degree, dec=composit_catalog2['dec_deg_1']*u.degree).position_angle(SkyCoord(ra=composit_catalog2['ra_deg_2']*u.degree, dec=composit_catalog2['dec_deg_2']*u.degree)).to(u.deg)
    sources_sep = SkyCoord(ra=composit_catalog2['ra_deg_1']*u.degree, dec=composit_catalog2['dec_deg_1']*u.degree).separation(SkyCoord(ra=composit_catalog2['ra_deg_2']*u.degree, dec=composit_catalog2['dec_deg_2']*u.degree)).to(u.arcsec)
    sources_mag_diff = composit_catalog2['mag_2'] - composit_catalog2['mag_1']
    
    composit_catalog2.add_column(sources_pa, name='theta_measured')
    composit_catalog2.add_column(sources_sep, name='rho_measured')
    composit_catalog2.add_column(sources_mag_diff, name='mag_diff')
    ds_sources = vstack([sources_ds, composit_catalog2])

    return ds_sources

# Define double star class
class gaia_calculated_attributes:
    def __init__(self, pairObjectId, pairDistanceMinA, pairDistanceMinB, starRa1, starDec1, starRa2, starDec2, starCoord1, starCoord2, starParallax1, starParallaxError1, possSep1, rhoStar, starName1,starName2, starParallax2, starParallaxError2, starActualRa1, starActualDec1, starActualRa2, starActualDec2, starActualCoord1, starActualCoord2, rhoActual, starDistance1, starDistanceMax1, starDistanceMin1, starDistanceRange1, starDistance2, starDistanceMax2, starDistanceMin2, starDistanceRange2, distanceCommon, pairParallaxFactor, pairPmFactor, pairPmCommon, pairAbsMag1, pairAbsMag2, pairLum1, pairLum2, pairAltLum1, pairAltLum2, pairRad1, pairRad2, pairDR3Theta, pairDR3Rho, pairMass1, pairMass2, pairBVIndexA, pairBVIndexB, pairSepPar2, pairDistance, pairSepPar, pairEscapeVelocity, pairRelativeVelocity, pairHarshawFactor, pairHarshawPhysicality, pairBinarity, pairDesignationA, pairDesignationB, pairRaA, pairDecA, pairRaB, pairDecB, pairMagA, pairMagB, pairMeanTheta, pairMeanThetaErr, pairMeanRho, pairMeanRhoErr, pairMagnitudeA, pairMagnitudeB, pairGMagDiff, pairMagDiff, pairMagDiffErr, pairRadVelA, pairRadVelErrA, pairRadVelB, pairRadVelErrB, pairRadVelRatioA, pairRadVelRatioB, pairDesA, pairDesB, dateOfObservation, pairACurrentCoord, pairBCurrentCoord, pairAMeasuredCoord, pairBMeasuredCoord, pairACoordErr, pairBCoordErr, preciseCoord, pair_orbit, gaiaData):
        self.pairObjectId = pairObjectId
        self.pairDistanceMinA = pairDistanceMinA
        self.pairDistanceMinB = pairDistanceMinB
        self.starRa1 = starRa1
        self.starDec1 = starDec1
        self.starRa2 = starRa2
        self.starDec2 = starDec2
        self.starCoord1 = starCoord1
        self.starCoord2 = starCoord2
        self.starParallax1 = starParallax1
        self.starParallaxError1 = starParallaxError1
        self.possSep1 = possSep1
        self.rhoStar = rhoStar
        self.starName1 = starName1
        self.starName2 = starName2
        self.starParallax2 = starParallax2
        self.starParallaxError2 = starParallaxError2
        self.starActualRa1 = starActualRa1
        self.starActualDec1 = starActualDec1
        self.starActualRa2 = starActualRa2
        self.starActualDec2 = starActualDec2
        self.starActualCoord1 = starActualCoord1
        self.starActualCoord2 = starActualCoord2
        self.rhoActual = rhoActual
        self.starDistance1 = starDistance1
        self.starDistanceMax1 = starDistanceMax1
        self.starDistanceMin1 = starDistanceMin1
        self.starDistanceRange1 = starDistanceRange1
        self.starDistance2 = starDistance2
        self.starDistanceMax2 = starDistanceMax2
        self.starDistanceMin2 = starDistanceMin2
        self.starDistanceRange2 = starDistanceRange2
        self.distanceCommon = distanceCommon
        self.pairParallaxFactor = pairParallaxFactor
        self.pairPmFactor = pairPmFactor
        self.pairPmCommon = pairPmCommon
        self.pairAbsMag1 = pairAbsMag1
        self.pairAbsMag2 = pairAbsMag2
        self.pairLum1 = pairLum1
        self.pairLum2 = pairLum2
        self.pairAltLum1 = pairAltLum1
        self.pairAltLum2 = pairAltLum2
        self.pairRad1 = pairRad1
        self.pairRad2 = pairRad2
        self.pairDR3Theta = pairDR3Theta
        self.pairDR3Rho = pairDR3Rho
        self.pairMass1 = pairMass1
        self.pairMass2 = pairMass2
        self.pairBVIndexA = pairBVIndexA
        self.pairBVIndexB = pairBVIndexB
        self.pairSepPar2 = pairSepPar2
        self.pairDistance = pairDistance
        self.pairSepPar = pairSepPar
        self.pairEscapeVelocity = pairEscapeVelocity
        self.pairRelativeVelocity = pairRelativeVelocity
        self.pairHarshawFactor = pairHarshawFactor
        self.pairHarshawPhysicality = pairHarshawPhysicality
        self.pairBinarity = pairBinarity
        self.pairDesignationA = pairDesignationA
        self.pairDesignationB = pairDesignationB
        self.pairRaA = pairRaA
        self.pairDecA = pairDecA
        self.pairRaB = pairRaB
        self.pairDecB = pairDecB
        self.pairMagA = pairMagA
        self.pairMagB = pairMagB
        self.pairMeanTheta = pairMeanTheta
        self.pairMeanThetaErr = pairMeanThetaErr
        self.pairMeanRho = pairMeanRho
        self.pairMeanRhoErr = pairMeanRhoErr
        self.pairMagnitudeA = pairMagnitudeA
        self.pairMagnitudeB = pairMagnitudeB
        self.pairGMagDiff = pairGMagDiff
        self.pairMagDiff = pairMagDiff
        self.pairMagDiffErr = pairMagDiffErr
        self.pairRadVelA = pairRadVelA
        self.pairRadVelErrA = pairRadVelErrA
        self.pairRadVelB = pairRadVelB
        self.pairRadVelErrB = pairRadVelErrB
        self.pairRadVelRatioA = pairRadVelRatioA
        self.pairRadVelRatioB = pairRadVelRatioB
        self.pairDesA = pairDesA
        self.pairDesB = pairDesB
        self.dateOfObservation = dateOfObservation
        self.pairACurrentCoord = pairACurrentCoord
        self.pairBCurrentCoord = pairBCurrentCoord
        self.pairAMeasuredCoord = pairAMeasuredCoord
        self.pairBMeasuredCoord = pairBMeasuredCoord
        self.pairACoordErr = pairACoordErr
        self.pairBCoordErr = pairBCoordErr
        self.preciseCoord = preciseCoord
        self.pair_orbit = pair_orbit
        self.gaiaData = gaiaData

class wds_attributes:
    def __init__(self, pairObjectId, starActualRa1, starActualDec1, starActualRa2, starActualDec2, pairMeanTheta, pairMeanThetaErr, pairMeanRho, pairMeanRhoErr, pairMagnitudeA, pairMagnitudeB, pairMagDiff, pairMagDiffErr, dateOfObservation, pairAMeasuredCoord, pairBMeasuredCoord, preciseCoord):
        self.pairObjectId = pairObjectId
        self.starActualRa1 = starActualRa1
        self.starActualDec1 = starActualDec1
        self.starActualRa2 = starActualRa2
        self.starActualDec2 = starActualDec2
        self.pairMeanTheta = pairMeanTheta
        self.pairMeanThetaErr = pairMeanThetaErr
        self.pairMeanRho = pairMeanRho
        self.pairMeanRhoErr = pairMeanRhoErr
        self.pairMagnitudeA = pairMagnitudeA
        self.pairMagnitudeB = pairMagnitudeB
        self.pairMagDiff = pairMagDiff
        self.pairMagDiffErr = pairMagDiffErr
        self.dateOfObservation = dateOfObservation
        self.pairAMeasuredCoord = pairAMeasuredCoord
        self.pairBMeasuredCoord = pairBMeasuredCoord
        self.preciseCoord = preciseCoord

def gaia_calculations(gaia_star_a, gaia_star_b, double_star):
    pairObjectId = double_star[0]['2000 Coord'] + double_star[0]['Discov'] + str(double_star[0]['Comp'])
    pairDistanceMinA = calcDistanceMin(float(gaia_star_a['parallax']), float(gaia_star_a['parallax_error']))
    pairDistanceMinB = calcDistanceMin(float(gaia_star_b['parallax']), float(gaia_star_b['parallax_error']))
    
    # Calculate physical binarity

    # Creating empty arrays for Star related calculations
    #Set input data
    starRa1 = float(gaia_star_a['ra'])
    starDec1 = float(gaia_star_a['dec'])
    starRa2 = float(gaia_star_b['ra'])
    starDec2 = float(gaia_star_b['dec'])
    starCoord1 = SkyCoord(starRa1, starDec1, unit="deg")
    starCoord2 = SkyCoord(starRa2, starDec2, unit="deg")
    starParallax1 = float(gaia_star_a['parallax'])
    starParallaxError1 = float(gaia_star_a['parallax_error'])
                
    # Calculate the widest possible separation for StarA
    possSep1 = possible_distance / calcDistanceMax(starParallax1, starParallaxError1)
    # rhoStar = rhoCalc(starRa1, starDec1, starRa2, starDec2)
    rhoStar = starCoord1.separation(starCoord2).arcsecond
    if possSep1 > rhoStar:
        #starId1 = gaia_star_a['solution_id']
        starName1 = gaia_star_a['designation']
        #starId2 = gaia_star_b['solution_id']
        starName2 = gaia_star_b['designation']
        starParallax2 = float(gaia_star_b['parallax'])
        starParallaxError2 = float(gaia_star_b['parallax_error'])
        starActualRa1 = float(double_star['ra_deg_1'].mean())
        starActualDec1 = float(double_star['dec_deg_1'].mean())
        starActualRa2 = float(double_star['ra_deg_2'].mean())
        starActualDec2 = float(double_star['dec_deg_2'].mean())
        starActualCoord1 = SkyCoord(starActualRa1, starActualDec1, unit="deg")
        starActualCoord2 = SkyCoord(starActualRa2, starActualDec2, unit="deg")

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
    
    # Calculate attributes
    pairParallaxFactor, pairPmFactor, pairPmFactor, pairPmCommon, pairAbsMag1, pairAbsMag2, pairLum1, pairLum2, pairRad1, pairRad2, pairDR3Theta, pairDR3Rho, pairMass1, pairMass2, pairBVIndexA, pairBVIndexB, pairSepPar, pairEscapeVelocity, pairRelativeVelocity, pairHarshawFactor, pairHarshawPhysicality, pairBinarity = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 

    pairParallaxFactor = (calcParallaxFactor(gaia_star_a['parallax'], gaia_star_b['parallax'])) * 100
    pairPmFactor = (calcPmFactor(gaia_star_a['pmra'], gaia_star_a['pmdec'], gaia_star_b['pmra'], gaia_star_b['pmdec']))
    pairPmCommon = calcPmCategory(pairPmFactor * 100)
    pairAbsMag1 = calcAbsMag(float(gaia_star_a['phot_g_mean_mag']), float(gaia_star_a['parallax'])) # Calculate Absolute magnitude
    pairAbsMag2 = calcAbsMag(float(gaia_star_b['phot_g_mean_mag']), float(gaia_star_b['parallax'])) # Calculate Absolute magnitude
    pairLum1 = calcLuminosity(pairAbsMag1)
    pairLum2 = calcLuminosity(pairAbsMag2)
    pairAltLum1 = calcLuminosityAlternate(pairAbsMag1)
    pairAltLum2 = calcLuminosityAlternate(pairAbsMag2)
    pairRad1 = calcRadius(pairLum1, gaia_star_a['teff_gspphot'])
    pairRad2 = calcRadius(pairLum2, gaia_star_b['teff_gspphot'])
    # pairDR3Theta = thetaCalc(deltaRa(gaia_star_a['ra'], gaia_star_b['ra'], gaia_star_b['dec']), deltaDec(gaia_star_b['dec'], gaia_star_a['dec'])) + addThetaValue
    pairDR3Theta = starCoord1.position_angle(starCoord2).degree
    # pairDR3Rho = rhoCalc(gaia_star_a['ra'], gaia_star_a['dec'], gaia_star_b['ra'], gaia_star_b['dec'])
    pairDR3Rho = rhoStar
    pairMass1 = calcMass(pairLum1)
    pairMass2 = calcMass(pairLum2)
    pairBVIndexA = float(gaia_star_a['phot_bp_mean_mag']) - float(gaia_star_a['phot_rp_mean_mag'])
    pairBVIndexB = float(gaia_star_b['phot_bp_mean_mag']) - float(gaia_star_b['phot_rp_mean_mag'])
    pairSepPar2 = sepCalc(pairDistanceMinA, pairDistanceMinB, rhoStar) # Separation of the pairs in parsecs
    pairDistance = calc_average_distance(float(gaia_star_a['parallax']), float(gaia_star_a['parallax_error']), float(gaia_star_b['parallax']), float(gaia_star_b['parallax_error']), pairDR3Rho)
    pairSepPar = pairDistance[2] * auToParsec
    print('pairSepPar: ', pairSepPar)
    print('pairSepPar2: ', pairSepPar2)
    print('pairDistance: ', pairDistance[1])

    pairEscapeVelocity = calcEscapevelocity(pairMass1, pairMass2, pairSepPar, gravConst)
    pairRelativeVelocity = calcRelativeVelocity(float(gaia_star_a['pmra']), float(gaia_star_a['pmdec']), float(gaia_star_b['pmra']), float(gaia_star_b['pmdec']), gaia_star_a['radial_velocity'], gaia_star_b['radial_velocity'], pairDistanceMinA, pairDistanceMinB)
    pairHarshawFactor = calcHarshaw((pairParallaxFactor) / 100, (pairPmFactor))
    pairHarshawPhysicality = calcHarshawPhysicality(pairHarshawFactor * 100)
    pairBinarity = calcBinarity(pairRelativeVelocity, pairEscapeVelocity)
    
    # Calculate values for each pair based on the groups
    pairDesignationA = gaia_star_a['designation']
    pairDesignationB = gaia_star_b['designation']
    pairRaA = gaia_star_a['ra']
    pairDecA = gaia_star_a['dec']
    pairRaB = gaia_star_b['ra']
    pairDecB = gaia_star_b['dec']
    pairMagA = gaia_star_a['phot_g_mean_mag']
    pairMagB = gaia_star_b['phot_g_mean_mag']
    pairMeanTheta = float(double_star['theta_measured'].degree.mean())
    pairMeanThetaErr = float(double_star['theta_measured'].degree.std())
    pairMeanRho = float(double_star['rho_measured'].arcsec.mean())
    pairMeanRhoErr = float(double_star['rho_measured'].arcsec.std())
    pairMagnitudeA = double_star['mag_1']
    pairMagnitudeB = double_star['mag_2']
    pairGMagDiff = float(gaia_star_b['phot_g_mean_mag']) - float(gaia_star_a['phot_g_mean_mag'])
    pairMagDiff = float((double_star['mag_diff']).mean())
    pairMagDiffErr = (double_star['mag_diff']).std()
    pairRadVelA = convertStringToNan(gaia_star_a['radial_velocity'])
    pairRadVelErrA = convertStringToNan(gaia_star_a['radial_velocity_error'])
    pairRadVelB = convertStringToNan(gaia_star_b['radial_velocity'])
    pairRadVelErrB = convertStringToNan(gaia_star_b['radial_velocity_error'])
    pairRadVelRatioA = math.fabs(float(convertStringToNan(gaia_star_a['radial_velocity_error'])) / float(convertStringToNan(gaia_star_a['radial_velocity']))) * 100
    pairRadVelRatioB = math.fabs(float(convertStringToNan(gaia_star_b['radial_velocity_error'])) / float(convertStringToNan(gaia_star_b['radial_velocity']))) * 100
    pairDesA = str(gaia_star_a['designation'])
    pairDesB = str(gaia_star_b['designation'])
    
    # dateOfObservation = getUTC(fitsFileDate)
    print('dateOfObservation: ', Time(double_star['image_date'].data))
    print('dateOfObservationMean: ', Time(double_star['image_date'].data).mean())
    #dateOfObservation = getUTC(Time(double_star['image_date']).mean())
    dateOfObservation = getUTC(Time(double_star['image_date'].data).mean())
    
    pairACurrentCoord = calcCurrentDR3Coord(dateOfObservation, pairRaA, pairDecA, float(gaia_star_a['pmra']), float(gaia_star_a['pmdec']))
    pairBCurrentCoord = calcCurrentDR3Coord(dateOfObservation, pairRaB, pairDecB, float(gaia_star_b['pmra']), float(gaia_star_b['pmdec']))
    pairAMeasuredCoord = SkyCoord(ra=double_star['ra_deg_1'].groups.aggregate(np.mean) * u.deg, dec=double_star['dec_deg_1'].groups.aggregate(np.mean) * u.deg)
    pairBMeasuredCoord = SkyCoord(ra=double_star['ra_deg_2'].groups.aggregate(np.mean) * u.deg, dec=double_star['dec_deg_2'].groups.aggregate(np.mean) * u.deg)
    pairACoordErr = pairACurrentCoord.separation(pairAMeasuredCoord)
    pairBCoordErr = pairBCurrentCoord.separation(pairBMeasuredCoord)
    # Caculate the common distance from Earth
    
    if (pairMass1 is not None and pairMass2 and pairSepPar is not None and pairDistance[1] is not None and pairACurrentCoord.separation(pairBCurrentCoord).arcsecond is not None and double_star[0]['Sep_f'] is not None and pairMeanRho is not None and double_star[0]['PA_f'] is not None and pairMeanTheta is not None and double_star[0]['Date (first)'] is not None and dateOfObservation):
        # OLD function to calculate the historical orbit values based on the first measurements found in WDS
        # pair_orbit = calc_historic_orbit(pairMass1, pairMass2, pairSepPar, pairDistance[1], pairACurrentCoord.separation(pairBCurrentCoord).arcsecond, double_star[0]['Sep_f'], pairMeanRho, double_star[0]['PA_f'], pairMeanTheta, double_star[0]['Date (first)'], dateOfObservation)
        
        # NEW function to calculate the historical orbit values based on the calculated PA and SEP from Gaia DR3 on epoch 2016
        pair_orbit = calc_historic_orbit(pairMass1, pairMass2, pairSepPar, pairDistance[1], pairACurrentCoord.separation(pairBCurrentCoord).arcsecond, SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).arcsecond, pairMeanRho, SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).degree, pairMeanTheta, gaia_dr3_epoch, dateOfObservation)
    else:
        pair_orbit = ['Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.']
    
    preciseCoord = str(getPreciseCoord(pairRaA, pairDecA, Time(double_star['image_date'].data).mean()))
    gaiaData = str(double_star[0]['2000 Coord']) + ',' + str(double_star[0]['Discov']) + ',' + str(gaia_star_a['pmra']) + ',' + str(gaia_star_a['pmdec']) + ',' + str(gaia_star_b['pmra']) + ',' + str(gaia_star_b['pmdec']) + ',' + str(gaia_star_a['parallax']) + ',' + str(gaia_star_b['parallax']) + ',' + str(calcDistance(gaia_star_a['parallax'])) + ',' + str(calcDistance(gaia_star_b['parallax'])) + ',' + str(gaia_star_a['radial_velocity']) + ',' + str(gaia_star_b['radial_velocity']) + ',' + 'pairRad1' + ',' + 'pairRad2' + ',' + str(pairLum1) + ',' + str(pairLum2) + ',' + str(gaia_star_a['teff_gspphot']) + ',' + str(gaia_star_b['teff_gspphot']) + ',' + str(gaia_star_a['phot_g_mean_mag']) + ',' + str(gaia_star_b['phot_g_mean_mag']) + ',' + str(gaia_star_a['phot_bp_mean_mag']) + ',' + str(gaia_star_b['phot_bp_mean_mag']) + ',' + str(gaia_star_a['phot_rp_mean_mag']) + ',' + str(gaia_star_b['phot_rp_mean_mag']) + ',' + str(pairDR3Theta) + ',' + str(pairDR3Rho) + ',' + str(gaia_star_a['ra']) + ',' + str(gaia_star_a['dec']) + ',' + str(gaia_star_b['ra']) + ',' + str(gaia_star_b['dec']) + ',' + str(gaia_star_a['parallax_error']) + ',' + str(gaia_star_b['parallax_error'])

    

    double_star_calculation_results = gaia_calculated_attributes(pairObjectId, pairDistanceMinA, pairDistanceMinB, starRa1, starDec1, starRa2, starDec2, starCoord1, starCoord2, starParallax1, starParallaxError1, possSep1, rhoStar, starName1,starName2, starParallax2, starParallaxError2, starActualRa1, starActualDec1, starActualRa2, starActualDec2, starActualCoord1, starActualCoord2, rhoActual, starDistance1, starDistanceMax1, starDistanceMin1, starDistanceRange1, starDistance2, starDistanceMax2, starDistanceMin2, starDistanceRange2, distanceCommon, pairParallaxFactor, pairPmFactor, pairPmCommon, pairAbsMag1, pairAbsMag2, pairLum1, pairLum2, pairAltLum1, pairAltLum2, pairRad1, pairRad2, pairDR3Theta, pairDR3Rho, pairMass1, pairMass2, pairBVIndexA, pairBVIndexB, pairSepPar2, pairDistance, pairSepPar, pairEscapeVelocity, pairRelativeVelocity, pairHarshawFactor, pairHarshawPhysicality, pairBinarity, pairDesignationA, pairDesignationB, pairRaA, pairDecA, pairRaB, pairDecB, pairMagA, pairMagB, pairMeanTheta, pairMeanThetaErr, pairMeanRho, pairMeanRhoErr, pairMagnitudeA, pairMagnitudeB, pairGMagDiff, pairMagDiff, pairMagDiffErr, pairRadVelA, pairRadVelErrA, pairRadVelB, pairRadVelErrB, pairRadVelRatioA, pairRadVelRatioB, pairDesA, pairDesB, dateOfObservation, pairACurrentCoord, pairBCurrentCoord, pairAMeasuredCoord, pairBMeasuredCoord, pairACoordErr, pairBCoordErr, preciseCoord, pair_orbit, gaiaData)

    print('### Gaia Calculations object ###')
    pprint.pprint(vars(double_star_calculation_results))

    return double_star_calculation_results

def wds_measurement(double_star):
    pairObjectId = double_star[0]['2000 Coord'] + double_star[0]['Discov'] + str(double_star[0]['Comp'])
    starActualRa1 = float(double_star['ra_deg_1'].mean())
    starActualDec1 = float(double_star['dec_deg_1'].mean())
    starActualRa2 = float(double_star['ra_deg_2'].mean())
    starActualDec2 = float(double_star['dec_deg_2'].mean())
    pairMeanTheta = float(double_star['theta_measured'].degree.mean())
    pairMeanThetaErr = float(double_star['theta_measured'].degree.std())
    pairMeanRho = float(double_star['rho_measured'].arcsec.mean())
    pairMeanRhoErr = float(double_star['rho_measured'].arcsec.std())
    pairMagnitudeA = double_star['mag_1']
    pairMagnitudeB = double_star['mag_2']
    pairMagDiff = float((double_star['mag_diff']).mean())
    pairMagDiffErr = (double_star['mag_diff']).std()
    dateOfObservation = getUTC(Time(double_star['image_date'].data).mean())
    pairAMeasuredCoord = SkyCoord(ra=double_star['ra_deg_1'].groups.aggregate(np.mean) * u.deg, dec=double_star['dec_deg_1'].groups.aggregate(np.mean) * u.deg)
    pairBMeasuredCoord = SkyCoord(ra=double_star['ra_deg_2'].groups.aggregate(np.mean) * u.deg, dec=double_star['dec_deg_2'].groups.aggregate(np.mean) * u.deg)
    preciseCoord = str(getPreciseCoord(starActualRa1, starActualDec1, Time(double_star['image_date'].data).mean()))
    wds_measurements = wds_attributes(pairObjectId, starActualRa1, starActualDec1, starActualRa2, starActualDec2, pairMeanTheta, pairMeanThetaErr, pairMeanRho, pairMeanRhoErr, pairMagnitudeA, pairMagnitudeB, pairMagDiff, pairMagDiffErr, dateOfObservation, pairAMeasuredCoord, pairBMeasuredCoord, preciseCoord)

    return wds_measurements