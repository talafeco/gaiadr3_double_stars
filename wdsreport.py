#! /usr/bin/python3
# WDS Report tool to measure double stars on astronomical images based on Gaia DR3 data
# Version: 1.0
# Usage: wdsreport <wds_file> <image_folder>

import csv
import os
import sys
import numpy as np
import datetime
import math
import sys
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky
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
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import warnings
from io import StringIO
from astropy.io import ascii
warnings.filterwarnings("ignore")
from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default

# Constant variables

# WDS table to be used to identify double stars on the image
# wdsTable = Table.read(sys.argv[1], delimiter=',', format='ascii')
wdsTable = Table.read(f"/usr/share/dr3map/dr3-wds/wdsweb_summ2.dat", format='ascii')
print(wdsTable.info)

# Set working directory to read Double Star images
workingDirectory = sys.argv[1]

#########################
### Declare functions ###
#########################

def convertStringToNan(str):
    if str == 'null' or str == '' or str == '.':
        str = np.nan
    return str

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
gravConst = 0.0043009 # Gravitational constant is convenient if measure distances in parsecs (pc), velocities in kilometres per second (km/s) and masses in solar units M
def calcEscapevelocity(mass_a, mass_b, separation, gravconst):
    escvel = math.sqrt((2 * gravconst * (mass_a + mass_b)) / separation)
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
def searchWds(pairadesig):
    wdsIdx = np.where(pairadesig == wdsFile['designation'])
    wdsRow = wdsFile[wdsIdx]
    if wdsRow:
        wdsPair = wdsRow[0]
    elif not wdsRow:
        wdsPair = 'Not found'
    return wdsPair

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
reportTable = QTable([reportw_identifier, reportdate, reporttheta, reportthetaerr, reportrho, reportrhoerr, reportmag_pri, reportmag_prierr, reportmag_sec, reportmag_secerr, reportfilter, reportfilterfwhm, reporttelescopeap, reportnights, reportrefcode, reporttech, reportcat], names=('wds_identifier', 'date_of_obs', 'mean_theta', 'mean_theta_err', 'mean_rho', 'mean_rho_err', 'mag_pri', 'mag_pri_err', 'mag_sec', 'mag_sec_err', 'filter', 'filter_fwhm', 'telescope_ap', 'nights_of_obs', 'reference_code', 'tech_code', 'catalog_code'), meta={'name': 'report table'})

print('\n### Creating filelist ###')


directoryContent = os.listdir(workingDirectory)
print('Working directory: ', workingDirectory)

files = [f for f in directoryContent if os.path.isfile(workingDirectory+'/'+f) and f.endswith('.new')]
print('Files:', files)

wdscatalog = SkyCoord(ra=wdsTable['ra_deg']*u.degree, dec=wdsTable['dec_deg']*u.degree)
print('WDS Catalog:\n', wdscatalog)

sources_ds = Table()

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
    ra2, dec2 = mywcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 1)
    sources.add_column(ra2, name='ra_deg') 
    sources.add_column(dec2, name='dec_deg')
    
    #print(sources.info)
    #print(sources)
    
    wds_catalog = SkyCoord(ra=wdsTable['ra_deg']*u.degree, dec=wdsTable['dec_deg']*u.degree)
    sources_catalog = SkyCoord(ra=sources['ra_deg']*u.degree, dec=sources['dec_deg']*u.degree)
    idxw, idxs, wsd2d, wsd3d = search_around_sky(wds_catalog, sources_catalog, 0.001*u.deg)
    composit_catalog = hstack([wdsTable[idxw]['wds_identifier', 'discovr', 'comp', 'theta', 'rho'], sources[idxs]['id', 'mag', 'ra_deg', 'dec_deg']])
    companion_catalog = SkyCoord(ra=composit_catalog['ra_deg'] * u.degree, dec=composit_catalog['dec_deg'] * u.degree).directional_offset_by(composit_catalog['theta'] * u.degree, composit_catalog['rho'] * u.arcsec)
    
    #print('Companion catalog\n', companion_catalog)
    #print('Composit catalog\n', composit_catalog)
    idxs2, d2ds2, d3ds2 = match_coordinates_sky(companion_catalog, sources_catalog)
    #print('Companion catalog idxs2\n', sources[idxs2])
    #print(idxs2)

    composit_catalog2 = hstack([composit_catalog, sources[idxs2]]) #['id', 'mag', 'ra_deg', 'dec_deg']
    # print(composit_catalog2.info)
    #print('Composit catalog 2.\n', composit_catalog2)

    sources_pa = SkyCoord(ra=composit_catalog2['ra_deg_1']*u.degree, dec=composit_catalog2['dec_deg_1']*u.degree).position_angle(SkyCoord(ra=composit_catalog2['ra_deg_2']*u.degree, dec=composit_catalog2['dec_deg_2']*u.degree)).to(u.deg)
    sources_sep = SkyCoord(ra=composit_catalog2['ra_deg_1']*u.degree, dec=composit_catalog2['dec_deg_1']*u.degree).separation(SkyCoord(ra=composit_catalog2['ra_deg_2']*u.degree, dec=composit_catalog2['dec_deg_2']*u.degree)).to(u.arcsec)
    sources_mag_diff = composit_catalog2['mag_2'] - composit_catalog2['mag_1']
    
    composit_catalog2.add_column(sources_pa, name='theta_measured')
    composit_catalog2.add_column(sources_sep, name='rho_measured')
    composit_catalog2.add_column(sources_mag_diff, name='mag_diff')


    #print('Matching sources list')
    #print(sources[idxs])
    #print(wdsTable[idxw])
    #print('Composit catalog 2. extended\n', composit_catalog2)

    sources_ds = vstack([sources_ds, composit_catalog2])

print('### Sources DS ###')
print(sources_ds)

upd_sources_ds = sources_ds[sources_ds['rho_measured'] != 0]
#upd_sources_ds.add_column(str(sources_ds['wds_identifier']) + '_' + str(sources_ds['discovr']) + '_' + str(sources_ds['comp']), name='object_id')
upd_sources_ds_by_object = upd_sources_ds.group_by(['wds_identifier', 'discovr', 'comp'])

# Create object ID array
# upd_sources_ds_object_id = str(sources_ds['wds_identifier']) + '_' + str(sources_ds['discovr']) + '_' + str(sources_ds['comp'])
#upd_sources_ds_object_id = np.empty(0,dtype=str)

#for line in sources_ds:
#    upd_sources_ds_object_id_instance = str(line['wds_identifier']) + '_' + str(line['discovr']) + '_' + str(line['comp'])
#    np.char.replace(upd_sources_ds_object_id_instance,' ','_')
#    print(upd_sources_ds_object_id_instance)
#    np.append(upd_sources_ds_object_id, upd_sources_ds_object_id_instance)

print(upd_sources_ds.info)
print('### Updated sources DS table grouped by WDS Identifier, Discoverer and Components ###')
print(upd_sources_ds_by_object)

#Ide kell egy source_id generáló algoritmus!!!

#print('\n### Double stars ###')
#print(dsTable)
#print('\n### Double stars table info ###')
#print(dsTable.info)       
        
### Search double stars on the image sequence
#dsTable_by_object = dsTable.group_by('object_id')
#print('\n### Report Table by object ###')
#print(dsTable_by_object)


objectMean = upd_sources_ds_by_object.groups.aggregate(np.mean)
#print(objectMean)

count = 1
for ds in upd_sources_ds_by_object.groups:
    print('\n### Group index:', count, '###')
    print(ds)
    count = count + 1
    pairObjectId = ds[0]['wds_identifier'] + ds[0]['discovr'] + str(ds[0]['comp'])
    
    # Search component in the Gaia DR3 database
    pairACoord = SkyCoord(ra=ds[0]['ra_deg_1'], dec=ds[0]['dec_deg_1'], unit=(u.degree, u.degree), frame='icrs')
    pairBCoord = SkyCoord(ra=ds[0]['ra_deg_2'], dec=ds[0]['dec_deg_2'], unit=(u.degree, u.degree), frame='icrs')
    a = Gaia.cone_search_async(pairACoord, radius=u.Quantity(0.001, u.deg))
    b = Gaia.cone_search_async(pairBCoord, radius=u.Quantity(0.001, u.deg))
    gaiaAStar = a.get_results()
    gaiaBStar = b.get_results()
    print(gaiaAStar[0]['DESIGNATION'])
    print(gaiaBStar[0]['DESIGNATION'])
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
                
    # Calculate the widest possible separation for StarA
    possSep1 = 10000 / calcDistanceMax(starParallax1, starParallaxError1)
    rhoStar = rhoCalc(starRa1, starDec1, starRa2, starDec2)
    if possSep1 > rhoStar:
        starId1 = gaiaAStar[0]['solution_id']
        starName1 = gaiaAStar[0]['DESIGNATION']
        starId2 = gaiaBStar[0]['solution_id']
        starName2 = gaiaBStar[0]['DESIGNATION']
        starParallax2 = float(gaiaBStar[0]['parallax'])
        starParallaxError2 = float(gaiaBStar[0]['parallax_error'])
        #starPmRa1 = float(gaiaAStar['pmra'])
        #starPmDec1 = float(gaiaAStar['pmdec'])
        #starPmRa2 = float(gaiaBStar['pmra'])
        #starPmDec2 = float(gaiaBStar['pmdec'])
        #starGMag1 = float(gaiaAStar['phot_g_mean_mag'])
        #starGMag2 = float(gaiaBStar['phot_g_mean_mag'])
        #starBpMag1 = float(gaiaAStar['phot_bp_mean_mag'])
        #starBpMag2 = float(gaiaBStar['phot_bp_mean_mag'])
        #starRpMag1 = float(gaiaAStar['phot_rp_mean_mag'])
        #starRpMag2 = float(gaiaBStar['phot_rp_mean_mag'])
        #starRadVel1 = float(gaiaAStar['radial_velocity'])
        #starRadVelErr1 = float(gaiaAStar['radial_velocity_error'])
        #starRadVel2 = float(gaiaBStar['radial_velocity'])
        #starRadVelErr2 = float(gaiaBStar['radial_velocity_error'])
        #starTemp1 = float(gaiaAStar['teff_gspphot'])
        #starTemp2 = float(gaiaBStar['teff_gspphot'])
        #starImageIdA = float(gaiaBStar[15])
        #starImageIdB = float(gaiaBStar[15])
        starActualRa1 = float(ds['ra_deg_1'].mean())
        starActualDec1 = float(ds[0]['dec_deg_1'].mean())
        starActualRa2 = float(ds[0]['ra_deg_2'].mean())
        starActualDec2 = float(ds[0]['dec_deg_2'].mean())
        #starObjectId = (str(starId1) + '_' + str(starId2))

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

        # Check if stars shares a common distance range

        distanceCommon = ()
        if starDistanceMin1 < starDistanceMin2 < starDistanceMax1 or starDistanceMin2 < starDistanceMin1 < starDistanceMax2:
            distanceCommon = 'overlapping'
        else:
            distanceCommon = 'no'
        
        # Check if the pair is a Common Proper Motion pairs (CPM), Similar Proper Motion (SPM) or Different Proper Motion (DPM)


        #Print data, if stars are close and share a common distance range
        #if distanceCommon == 'overlapping':
        #    reportTable.add_row([star[0], starId1, starName1, starRa1, starDec1, starParallax1, starParallaxError1, starPmRa1, starPmDec1, starGMag1, starBpMag1, starRpMag1, starRadVel1, starRadVelErr1, starTemp1, starActualRa1, starActualDec1, starActualMag1, starId2, starName2, starRa2, starDec2, starParallax2, starParallaxError2, starPmRa2, starPmDec2, starGMag2, starBpMag2, starRpMag2, starRadVel2, starRadVelErr2, starTemp2, starActualRa2, starActualDec2, starActualMag2, thetaStar, thetaActual, rhoStar, rhoActual, starObjectId])
    
    
    pairParallaxFactor = (calcParallaxFactor(gaiaAStar[0]['parallax'], gaiaBStar[0]['parallax'])) * 100
    pairPmFactor = (calcPmFactor(gaiaAStar[0]['pmra'], gaiaAStar[0]['pmdec'], gaiaBStar[0]['pmra'], gaiaBStar[0]['pmdec'])) * 100
    pairPmCommon = calcPmCategory(pairPmFactor)
    pairAbsMag1 = calcAbsMag(gaiaAStar[0]['phot_g_mean_mag'], gaiaAStar[0]['parallax']) # Calculate Absolute magnitude
    pairAbsMag2 = calcAbsMag(gaiaBStar[0]['phot_g_mean_mag'], gaiaBStar[0]['parallax']) # Calculate Absolute magnitude
    pairLum1 = calcLuminosity(pairAbsMag1)
    pairLum2 = calcLuminosity(pairAbsMag2)
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
    pairMeanTheta = ds['theta_measured'].degree.mean()
    pairMeanThetaErr = ds['theta_measured'].degree.std()
    pairMeanRho = ds['rho_measured'].arcsec.mean()
    pairMeanRhoErr = ds['rho_measured'].arcsec.std()
    pairMagnitudeA = ds[0]['mag_1']
    pairMagnitudeB = ds[0]['mag_2']
    pairMagDiff = (ds['mag_diff']).mean()
    pairMagDiffErr = (ds['mag_diff']).std()
    reportName = (workingDirectory + '/' + pairObjectId + '.txt')
    reportFile = open(reportName, "a")
    
    # Print temp data
    print('### COMPONENTS ###')
    print('\nWDS Identifier:', ds[0]['wds_identifier'], ds[0]['discovr'], ds[0]['comp'])
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
    print('Proper motion factor:', pairPmFactor, '%')
    print('Proper motion category:', pairPmCommon)
    print('Absolute magnitude A:', pairAbsMag1)
    print('Absolute magnitude B:', pairAbsMag2)
    print('Luminosity A:', pairLum1)
    print('Luminosity B:', pairLum2)
    print('Mass A:', pairMass1)
    print('Mass B:', pairMass2)
    print('BV index A:', pairBVIndexA, 'B:', pairBVIndexB)
    print('Radial velocity of the stars', 'A:', gaiaAStar[0]['radial_velocity'], 'km/s (Err:', gaiaAStar[0]['radial_velocity_error'], 'km/s)', 'B:', gaiaBStar[0]['radial_velocity'], 'km/s (Err:', gaiaBStar[0]['radial_velocity_error'], 'km/s)')
    print('Radial velocity ratio A:', math.fabs(gaiaAStar[0]['radial_velocity_error'] / gaiaAStar[0]['radial_velocity']) * 100, '%')
    print('Radial velocity ratio B:', math.fabs(gaiaBStar[0]['radial_velocity_error'] / gaiaBStar[0]['radial_velocity']) * 100, '%')
    print('Separation:', pairSepPar, 'parsec,', pairSepPar * 206265, 'AU')
    print('Pair Escape velocity:', pairEscapeVelocity, 'km/s')
    print('Pair Relative velocity:', pairRelativeVelocity, 'km/s')
    print('Pair Harshaw factor:', pairHarshawFactor)
    print('Pair Harshaw physicality:', pairHarshawPhysicality)
    print('Pair binarity:', pairBinarity)
    
    # Write results to file
    reportTable.add_row([ds[0]['wds_identifier'] + ds[0]['discovr'] + str(ds[0]['comp']), 'Date of observation', pairMeanTheta, pairMeanThetaErr, pairMeanRho, pairMeanRhoErr, np.nan, np.nan, pairMagDiff, pairMagDiffErr, 'Filter wawelenght', 'filter FWHM', '0.2', '1', 'TAL_2022', 'C', '7'])
    reportFile.write('\n\nWDS Identifier: ' + ds[0]['wds_identifier'])
    reportFile.write('\nDiscoverer and components: ' + str(ds[0]['discovr']) + ' ' + str(ds[0]['comp']))
    reportFile.write('\nMagnitude(s) (Pri / Sec): ' + str(ds[0]['mag_1']) + ' / ' +  str(ds[0]['mag_2']))
    reportFile.write('\nPA, Sep: ' + str(ds[0]['theta']) + ' / ' +  str(ds[0]['rho']))
    reportFile.write('\nPosition angle:')
    reportFile.write('\nTheta measurements' + str(ds['theta_measured'].degree))
    reportFile.write('\nMean: ' + str(pairMeanTheta))
    reportFile.write('\nError: ' + str(pairMeanThetaErr))
    reportFile.write('\nSeparation:')
    reportFile.write('\nRho measurements\n' + str(ds['rho_measured'].arcsec))
    reportFile.write('\nMean: ' + str(pairMeanRho))
    reportFile.write('\nError: ' + str(pairMeanRhoErr))
    reportFile.write('\nMagnitude measurements\n'  + str(ds['mag_diff']))
    reportFile.write('\nMean: ' + str(pairMagDiff))
    reportFile.write('\nError: ' + str(pairMagDiffErr))
    reportFile.write('\n\nParallax factor: ' + str(pairParallaxFactor) + ' %')
    reportFile.write('\nProper motion factor: ' + str(pairPmFactor) + ' %')
    reportFile.write('\nProper motion category: '+ str(pairPmCommon))
    reportFile.write('\nAbsolute magnitude A: ' + str(pairAbsMag1))
    reportFile.write('\nAbsolute magnitude B: ' + str(pairAbsMag2))
    reportFile.write('\nLuminosity A: ' + str(pairLum1))
    reportFile.write('\nLuminosity B: ' + str(pairLum2))
    reportFile.write('\nMass A: ' + str(pairMass1))
    reportFile.write('\nMass B: ' + str(pairMass2))
    reportFile.write('\nBV index A: ' + str(pairBVIndexA) + 'B: ' + str(pairBVIndexB))
    reportFile.write('\nRadial velocity of the stars A: ' + str(gaiaAStar[0]['radial_velocity']) + ' km/s (Err: ' + str(gaiaAStar[0]['radial_velocity_error']) + ' km/s) B: ' + str(gaiaBStar[0]['radial_velocity']) + ' km/s (Err: ' + str(gaiaBStar[0]['radial_velocity_error']) + ' km/s)')
    reportFile.write('\nRadial velocity ratio A: ' + str((math.fabs(gaiaAStar[0]['radial_velocity_error'] / gaiaAStar[0]['radial_velocity'])) * 100) + ' %')
    reportFile.write('\nRadial velocity ratio B: ' + str((math.fabs(gaiaBStar[0]['radial_velocity_error'] / gaiaBStar[0]['radial_velocity'])) * 100) + ' %')
    reportFile.write('\nSeparation: ' + str(pairSepPar) + ' parsec, ' + str((pairSepPar * 206265)) + ' AU')
    reportFile.write('\nPair Escape velocity: ' + str(pairEscapeVelocity) + ' km/s')
    reportFile.write('\nPair Relative velocity: ' + str(pairRelativeVelocity) + ' km/s')
    reportFile.write('\nPair Harshaw factor: ' + str(pairHarshawFactor))
    reportFile.write('\nPair Harshaw physicality: ' + str(pairHarshawPhysicality))
    reportFile.write('\nPair binarity: ' + str(pairBinarity))
    reportFile.write('\n\n### WDS form:\n')
    wdsform = str(ds[0]['wds_identifier']) + ',' + 'Date of observation' + ',' +  str(pairMeanTheta) + ',' +  str(pairMeanThetaErr) + ',' +  str(pairMeanRho) + ',' +  str(pairMeanRhoErr) + ',' +  'nan' + ',' +  'nan' + ',' +  str(pairMagDiff) + ',' +  str(pairMagDiffErr) + ',' + 'Filter wawelenght' + ',' + 'filter FWHM' + ',' + '0.2' + ',' + '1' + ',' + 'TAL_2022' + ',' +  'C' + ',' + '7'
    #print(str(wdsform))
    reportFile.write(str(wdsform))
    
#print(reportTable)
upd_sources_ds_by_object.write('double_stars.txt', format='ascii', overwrite=True, delimiter=',')
reportTable.write('double_stars_wds_format.txt', format='ascii', overwrite=True, delimiter=',')