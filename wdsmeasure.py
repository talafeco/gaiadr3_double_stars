#!/usr/bin/python3
# WDS Report tool to measure double stars on astronomical images based on Gaia DR3 data
# Version: 1.0
# Usage: wdsreport <image_folder>

# Tasks
# Calculate the double star's precise coordinates for current epoch
# Get all double stars from wds, which should be found on the images
# Search double stars based on magnitude difference, separation and position angle, if not found based on the coordinates


# 2024-07-27: updated only data[0] is counted due to rgb images
# check color fits file, then slice and measure independently: https://astronomy.stackexchange.com/questions/32222/converting-an-rgb-image-to-fits-astropy
# Check, which double stars should be on the image?
# it there should be a double star, but not found, search for similar separation + similar mag difference doubles
# design the measurement and calculation functions independently, process the images and write measurements, if there is no response from gaia in 10s or the star not found
# Upodate gaia search, to get a minimal distance in arcsecs

# 40-50 ivmásodperccel beljebb mérjen
# Listázza ki fileba a képen szereplő kettősöket

import os
import sys
import numpy as np
import sys
import datetime
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy.coordinates import match_coordinates_sky, search_around_sky
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack, hstack
from photutils.detection import DAOStarFinder
from astropy.wcs import WCS
import astropy.units as u
import warnings
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta
warnings.filterwarnings("ignore")
from regions import RegularPolygonPixelRegion, RectangleSkyRegion, RectanglePixelRegion
 
# Configuration for the SeeStar camera

'''dao_sigma = 3.0
dao_fwhm = 14.0
dao_threshold = 5.0'''

# Configuration for the CANON camera

'''dao_sigma = 3.0
dao_fwhm = 8.0
dao_threshold = 9.0'''


# Configuration for the ATIK camera

dao_sigma = 3.0
dao_fwhm = 8.0
dao_threshold = 12.0


# Configurations for calculations
search_cone = 0.001 # Decimal degree
dummyObservationDate = "2022-01-01T12:00:00"

# Gravitational constant is convenient if measure distances in parsecs (pc), velocities in kilometres per second (km/s) and masses in solar units M
image_limit = 2000

# Insert the downloaded wds file path here
wds_file = "/usr/share/dr3map/wds/wdsweb_summ2.txt"

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

# Function to calculate the Standard error in RA/DEC measurements
def calcStandardError(arr):
    stderr = np.std(arr)
    return stderr

    
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

# Create Image plot of the double stars
def imagePlot(filename, pairname, raa, deca, rab, decb):
    image_data = fits.open(filename)
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
    plt.savefig(str(workingDirectory + '/' + pairname + '_img.jpg').replace(' ', ''),dpi=300.0, bbox_inches='tight', pad_inches=0.2)
    plt.close()

def calculate_photo_center(wcs, header):
    photo_left_upper = SkyCoord.from_pixel(0, 0, wcs, origin=0, mode='all')
    photo_right_lower = SkyCoord.from_pixel(header['NAXIS2'], header['NAXIS1'], wcs, origin=0, mode='all')
    center = SkyCoord(header['CRVAL1'] * u.degree, header['CRVAL2'] * u.degree)
    radius = photo_left_upper.separation(photo_right_lower) / 2
    print('Center of photo: ', center.to_string('hmsdms'), '/', center.to_string('decimal'),
      '\nRadius of photo: ', radius)
    return center, radius


def define_image_plane(wcs, header):
    photo_center = SkyCoord(header['CRVAL1'] * u.degree, header['CRVAL2'] * u.degree)
    photo_left_upper = SkyCoord.from_pixel(0, 0, wcs, origin=0, mode='all')
    photo_left_lower = SkyCoord.from_pixel(header['NAXIS2'], 0, wcs, origin=0, mode='all')
    photo_right_upper = SkyCoord.from_pixel(0, header['NAXIS1'], wcs, origin=0, mode='all')
    photo_right_lower = SkyCoord.from_pixel(header['NAXIS2'], header['NAXIS1'], wcs, origin=0, mode='all')
    half_width = photo_left_upper.separation(photo_right_upper) / 2
    half_height = photo_left_upper.separation(photo_left_lower) / 2
    corners = SkyCoord(
        ra=[photo_center.ra - half_width, photo_center.ra + half_width,
            photo_center.ra + half_width, photo_center.ra - half_width],
        dec=[photo_center.dec - half_height, photo_center.dec - half_height,
            photo_center.dec + half_height, photo_center.dec + half_height]
        )
    
    return corners


def catalog_search_in_image(wcs, header, center, radius, wds_catalog_list):
    #image_region = define_image_region(wcs, header)
    d2d = center.separation(wds_catalog_list)
    catalogmsk = d2d < radius
    idxcatalog = np.where(catalogmsk)[0]
    image_plane = define_image_plane(wcs, header)

    in_fov = []
    for obj in wds_catalog_list[idxcatalog]:
        if (np.min(image_plane.ra) <= obj.ra <= np.max(image_plane.ra)) and (np.min(image_plane.dec) <= obj.dec <= np.max(image_plane.dec)):
            in_fov.append(wds_catalog_list)
    
    return wdsTable[idxcatalog]

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

wds_catalog = SkyCoord(ra=wdsTable['Coord (RA) hms'], dec=Angle(wdsTable['Coord (DEC) dms']), unit='hour, degree', frame="icrs") #

print('WDS Catalog created successfully - ' , datetime.datetime.now())

sources_ds = Table()


file_counter = 0

### Run source detection, collect star data to Qtable
print('\n### Running source detection ###\n' , datetime.datetime.now())
for fitsFile in files:
    # 1. Read the list of sources extracted from an image (fits) file
    file_counter = file_counter + 1
    print('\n\n### Processing file', file_counter, 'out of', len(files),': ', fitsFile, '###')
    fitsFileName = workingDirectory + '/' + fitsFile
    hdu = fits.open(fitsFileName)
    mywcs = WCS(hdu[0].header, naxis=2)
    file_header = hdu[0].header
    
	# Set observation date and time
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
    # data[0]
    data = hdu[0].data
    mean, median, std = sigma_clipped_stats(data, sigma = dao_sigma)  

    daofind = DAOStarFinder(fwhm=dao_fwhm, threshold=dao_threshold*std)  
    sources = daofind(data - median)
    ra2, dec2 = mywcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 1)
    sources.add_column(ra2, name='ra_deg') 
    sources.add_column(dec2, name='dec_deg')
    sources.add_column(fitsFileDate, name='image_date')
    sources.add_column(fitsFileName, name='file')

    
    photo_center, photo_radius = calculate_photo_center(mywcs, file_header)
    # doubles_on_photo = get_objects_from_catalog(wds_catalog, photo_center, photo_radius)
    # print(wdsTable[doubles_on_photo])

    print(catalog_search_in_image(mywcs, file_header, photo_center, photo_radius, wds_catalog))

    sources_catalog = SkyCoord(ra=sources['ra_deg']*u.degree, dec=sources['dec_deg']*u.degree, frame='fk5')
    idxw, idxs, wsd2d, wsd3d = search_around_sky(wds_catalog, sources_catalog, search_cone*u.deg)
    composit_catalog = hstack([wdsTable[idxw]['2000 Coord', 'Discov', 'Comp', 'Date (first)', 'PA_f', 'PA_l', 'Sep_f', 'Sep_l', 'Mag_A', 'Mag_B'], sources[idxs]['id', 'mag', 'ra_deg', 'dec_deg']])
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

upd_sources_ds = sources_ds[sources_ds['rho_measured'] != 0]
upd_sources_ds_by_object = upd_sources_ds.group_by(['2000 Coord', 'Discov', 'Comp'])

print('### Updated sources DS table grouped by WDS Identifier, Discoverer and Components ###')

upd_sources_ds_by_object.write(workingDirectory + '/double_stars.csv', format='ascii', overwrite=True, delimiter=',')

objectMean = upd_sources_ds_by_object.groups.aggregate(np.mean)

count = 1
for ds in upd_sources_ds_by_object.groups:
    print('\n#---------------------------------------------------------------------------------------------------------------------#')
    print('\n### Group index:', count, '###')
    # print(ds)
    count = count + 1
    pairObjectId = ds[0]['2000 Coord'] + ds[0]['Discov'] + str(ds[0]['Comp'])

    starActualRa1 = float(ds['ra_deg_1'].mean())
    starActualDec1 = float(ds['dec_deg_1'].mean())
    starActualRa2 = float(ds['ra_deg_2'].mean())
    starActualDec2 = float(ds['dec_deg_2'].mean())
         
    # Calculate attributes
    pairParallaxFactor, pairPmFactor, pairPmFactor, pairPmCommon, pairAbsMag1, pairAbsMag2, pairLum1, pairLum2, pairRad1, pairRad2, pairDR3Theta, pairDR3Rho, pairMass1, pairMass2, pairBVIndexA, pairBVIndexB, pairSepPar, pairEscapeVelocity, pairRelativeVelocity, pairHarshawFactor, pairHarshawPhysicality, pairBinarity = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 

    pairMeanTheta = float(ds['theta_measured'].degree.mean())
    pairMeanThetaErr = float(ds['theta_measured'].degree.std())
    pairMeanRho = float(ds['rho_measured'].arcsec.mean())
    pairMeanRhoErr = float(ds['rho_measured'].arcsec.std())
    pairMagnitudeA = ds['mag_1']
    pairMagnitudeB = ds['mag_2']
    pairMagDiff = float((ds['mag_diff']).mean())
    pairMagDiffErr = (ds['mag_diff']).std()
    dateOfObservation = getUTC(Time(ds['image_date'].data).mean())
    pairAMeasuredCoord = SkyCoord(ra=ds['ra_deg_1'].groups.aggregate(np.mean) * u.deg, dec=ds['dec_deg_1'].groups.aggregate(np.mean) * u.deg)
    pairBMeasuredCoord = SkyCoord(ra=ds['ra_deg_2'].groups.aggregate(np.mean) * u.deg, dec=ds['dec_deg_2'].groups.aggregate(np.mean) * u.deg)
    
    preciseCoord = str(getPreciseCoord(starActualRa1, starActualDec1, fitsFileDate))
    reportName = (workingDirectory + '/' + pairObjectId + '.txt').replace(' ', '')
    reportFile = open(reportName, "a")

    firstFitsImageFileName = ds['file'][0]
    imagePlot(firstFitsImageFileName, pairObjectId, starActualRa1, starActualDec1, starActualRa2, starActualDec2)

    # Print temp data
    print('\n### COMPONENTS ###')
    print('WDS Identifier:', ds[0]['2000 Coord'], ds[0]['Discov'], ds[0]['Comp'])
    print('Magnitude Pri: ' + str(ds[0]['Mag_A']))
    print('Magnitude Sec: ' + str(ds[0]['Mag_B']))
    print('PA last: ', str(ds[0]['PA_l']))
    print('Sep last: ',  str(ds[0]['Sep_l']))
    print('Date of observation (human readable): ', Time(ds['image_date'].data).mean())
    print('Date of observation (Julian date): ' + dateOfObservation)
    print('Precise coordinates (J2000): ' + preciseCoord)
    print('\nTheta measurements\n') # , ds['dspaactual']
    print('Mean:', pairMeanTheta)
    print('Error:', pairMeanThetaErr)
    print('\nRho measurements\n') # , ds['dssepactual']
    print('Mean:', pairMeanRho)
    print('Error:', pairMeanRhoErr)
    print('\nMagnitude difference measurements\n') # , ds['dsmagdiff']
    print('Mean:', pairMagDiff)
    print('Error:', pairMagDiffErr)
    
    # Write results to file
    reportTable.add_row([ds[0]['2000 Coord'] + ds[0]['Discov'] + str(ds[0]['Comp']), dateOfObservation, pairMeanTheta, pairMeanThetaErr, pairMeanRho, pairMeanRhoErr, np.nan, np.nan, pairMagDiff, pairMagDiffErr, 'Filter wawelenght', 'filter FWHM', '0.2', '1', 'TLB_2024', 'C', '7', preciseCoord])
    reportFile.write('### WDS Data ###')
    reportFile.write('\nWDS Identifier: ' + ds[0]['2000 Coord'])
    reportFile.write('\nDiscoverer and components: ' + str(ds[0]['Discov']) + ' ' + str(ds[0]['Comp']))
    reportFile.write('\nMagnitude Pri: ' + str(ds[0]['Mag_A']))
    reportFile.write('\nMagnitude Sec: ' + str(ds[0]['Mag_B']))
    reportFile.write('\nPA last: ' + str(ds[0]['PA_l']))
    reportFile.write('\nSep last: ' +  str(ds[0]['Sep_l']))
    reportFile.write('\n\n### Measurements ###')
    reportFile.write('\nDate of observation (human readable): ' + str(Time(ds['image_date'].data).mean()))
    reportFile.write('\nDate of observation: ' + dateOfObservation)
    reportFile.write('\nPrecise coordinates (J2000): ' + preciseCoord)
    reportFile.write('\n\nPosition angle:')
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

    
    reportFile.write('\n\n### WDS form:\n')
    wdsform = str(ds[0]['2000 Coord']) + ',' + dateOfObservation + ',' +  str(roundNumber(pairMeanTheta)) + ',' +  str(roundNumber(pairMeanThetaErr)) + ',' +  str(roundNumber(pairMeanRho)) + ',' +  str(roundNumber(pairMeanRhoErr)) + ',' +  'nan' + ',' +  'nan' + ',' +  str(roundNumber(pairMagDiff)) + ',' +  str(roundNumber(pairMagDiffErr)) + ',' + 'Filter wawelenght' + ',' + 'filter FWHM' + ',' + '0.2' + ',' + '1' + ',' + 'TLB_2023' + ',' +  'C' + ',' + '7'+ ',' + str(getPreciseCoord(starActualRa1, starActualDec1, fitsFileDate))
    reportFile.write(str(wdsform))

reportTable.write(workingDirectory + '/double_stars_wds_format.txt', format='ascii', overwrite=True, delimiter=',')