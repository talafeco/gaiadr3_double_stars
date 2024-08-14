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
hipparcos_file = Table.read(f"/usr/share/dr3map/hipparcos/I_239_selection.csv", format='ascii')

# Insert the downloaded wds file path here
wds_file = "/usr/share/dr3map/wds/wdsweb_summ2.txt"


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
        radius = math.sqrt(luminosity / ((teffStar / teffSun) ** 4))
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
    radspeeddif = math.fabs(radvela - radvelb)
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
        savename = str(working_directory + '/' + pairname + '_hrd.jpg').replace(' ', '')
        plt.savefig(savename, bbox_inches='tight', dpi=300.0)
        plt.close()
    else:
        print('Data is missiong, HRD plot cannot be created!')


# Create Image plot of the double stars
def imagePlot(filename, working_directory, pairname, raa, deca, rab, decb):
    # data[0]
    image_data = fits.open(working_directory + '/' + filename)
    #image_data = fits.open(filename)
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