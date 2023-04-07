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

# Constant variables

# WDS table to be used to identify double stars on the image
# wdsTable = Table.read(sys.argv[1], delimiter=',', format='ascii')
wdsTable = Table.read(f"/usr/share/dr3map/dr3-wds/wdsweb_summ2.dat", format='ascii')
print(wdsTable.info)

# Set working directory to read Double Star images
workingDirectory = sys.argv[1]

def convertStringToNan(str):
    if str == '' or str == '.':
        str = np.nan
    return str

def rhoCalc(raa, deca, rab, decb):
    rhocalc = math.sqrt(((raa-rab) * math.cos(math.radians(deca))) ** 2 + (deca - decb) ** 2) * 3600
    return rhocalc


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

    composit_catalog2.add_column(sources_pa, name='theta_measured')
    composit_catalog2.add_column(sources_sep, name='rho_measured')


    #print('Matching sources list')
    #print(sources[idxs])
    #print(wdsTable[idxw])
    #print('Composit catalog 2. extended\n', composit_catalog2)

    sources_ds = vstack([sources_ds, composit_catalog2])

    for star in sources:
        mainstar = SkyCoord(ra=star['ra_deg']*u.degree, dec=star['dec_deg']*u.degree)
        idx, d2d, d3d = match_coordinates_sky(mainstar, wdscatalog)
        catalogstar = SkyCoord(ra=wdsTable[idx]['ra_deg']*u.degree, dec=wdsTable[idx]['dec_deg']*u.degree)
        sep = catalogstar.separation(mainstar)
        if sep < Angle(0.002 * u.deg):
            companion = mainstar.directional_offset_by(wdsTable[idx][3] * u.degree, (wdsTable[idx][4] / 3600) * u.deg)
            separation = mainstar.separation(companion)
            for star2 in sources:
                ra3, dec3 = mywcs.all_pix2world([[star2['xcentroid'], star2['ycentroid']]], 0)[0]   
                compstar = SkyCoord(ra=ra3*u.degree, dec=dec3*u.degree)  
                sep2 = compstar.separation(companion)
                if sep2 <= Angle(0.002 * u.deg) and star != star2:
                    #print('\nWDS Identifier:', wdsTable[idx][0], wdsTable[idx][1], wdsTable[idx][2])
                    #print('Main star:', wdsTable[idx][0], wdsTable[idx][1], wdsTable[idx][3], wdsTable[idx][4])
                    #print('Main star data:\n', star)
                    #print('Separation:', Angle(sep))
                    #print('Comapnion star data:\n', star2)
                    paactual = mainstar.position_angle(compstar).to(u.deg)
                    sepactual = mainstar.separation(compstar) * 3600
                    #print('\nPA:', paactual.degree, 'Sep:', (sepactual.degree))
                    objectid = str(wdsTable[idx][0]) + '_' + str(wdsTable[idx][1]) + '_' + str(wdsTable[idx][2])
                    wdsmag_diff =  float(wdsTable[idx]['mag_sec']) - float(wdsTable[idx]['mag_pri'])
                    magdiff = star2['mag'] - star['mag']
                    dsTable.add_row([wdsTable[idx][0], wdsTable[idx][1], wdsTable[idx][2], wdsTable[idx][3], wdsTable[idx][4], wdsTable[idx][5], wdsTable[idx][6], wdsmag_diff, wdsTable[idx][7], wdsTable[idx][8], wdsTable[idx][9], wdsTable[idx][10], wdsTable[idx][11], str(wdsTable[idx][12]), str(wdsTable[idx][13]), str(wdsTable[idx][14]), str(wdsTable[idx][15]), objectid, paactual.degree, sepactual.degree, magdiff])

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

#print('### upd_sources_ds_object_id ###')
#print(upd_sources_ds_object_id)
# Add object ID array to upd_sources_ds table
# upd_sources_ds.add_column(upd_sources_ds_object_id, name='object_id')

print('### Updated sources DS ###')
print(upd_sources_ds)
print('### Updated sources DS table info ###')
print(upd_sources_ds.info)
print('### Updated sources DS table grouped by WDS Identifier, Discoverer and Components ###')
print(upd_sources_ds_by_object)

#Ide kell egy source_id generáló algoritmus!!!

print('\n### Double stars ###')
print(dsTable)
print('\n### Double stars table info ###')
print(dsTable.info)       
        
### Search double stars on the image sequence
dsTable_by_object = dsTable.group_by('object_id')
#print('\n### Report Table by object ###')
#print(dsTable_by_object)


objectMean = upd_sources_ds_by_object.groups.aggregate(np.mean)
#print(objectMean)

count = 1
for ds in upd_sources_ds_by_object.groups:
    #print('\n### Group index:', count, '###')
    #print(ds)
    count = count + 1
    pairObjectId = ds[0]['wds_identifier'] + ds[0]['discovr'] + ds[0]['comp']
    print(ds[0]['wds_identifier'] + ds[0]['discovr'] + ds[0]['comp'])
    print('Theta: ', ds['theta_measured'].degree)
    print('Rho: ', ds['rho_measured'].arcsec)
    print('Mag A: ', ds[0]['mag_1'])
    print('Mag B: ', ds[0]['mag_2'])
    print('Mag Diff.: ', ds[0]['mag_2'] - ds[0]['mag_1'])

    pairMeanTheta = ds['theta_measured'].degree.mean()
    pairMeanThetaErr = ds['theta_measured'].degree.std()
    pairMeanRho = ds['rho_measured'].arcsec.mean()
    pairMeanRhoErr = ds['rho_measured'].arcsec.std()
    pairMagnitudeA = ds[0]['mag_1']
    pairMagnitudeB = ds[0]['mag_2']
    pairMagDiff = (ds[0]['mag_2'] - ds[0]['mag_1']).mean()
    pairMagDiffErr = (ds[0]['mag_2'] - ds[0]['mag_1']).mean()
    pairMagDiffWDS = float(ds[0]['theta']) - ds[0]['rho']
    reportName = (workingDirectory + '/' + pairObjectId + '.txt')
    reportFile = open(reportName, "a")
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
    reportTable.add_row([ds[0]['wds_identifier'], 'Date of observation', pairMeanTheta, pairMeanThetaErr, pairMeanRho, pairMeanRhoErr, np.nan, np.nan, pairMagDiff, pairMagDiffErr, 'Filter wawelenght', 'filter FWHM', '0.2', '1', 'TAL_2022', 'C', '7'])
    reportFile.write('\n\nWDS Identifier: ' + ds[0]['wds_identifier'])
    reportFile.write('\nDiscoverer and components: ' + str(ds[0]['discovr']) + ' ' + str(ds[0]['comp']))
    reportFile.write('\nMagnitude(s) (Pri / Sec): ' + str(ds[0]['mag_1']) + ' / ' +  str(ds[0]['mag_2']))
    reportFile.write('\nPA, Sep: ' + str(ds[0]['theta']) + ' / ' +  str(ds[0]['rho']))
    reportFile.write('\nPosition angle: \n')
    reportFile.write('\nMean: ' + str(pairMeanTheta))
    reportFile.write('\nError: ' + str(pairMeanThetaErr))
    reportFile.write('\nSeparation: \n')
    reportFile.write('\nMean: ' + str(pairMeanRho))
    reportFile.write('\nError: ' + str(pairMeanRhoErr))
    reportFile.write('\nMagnitude measurements\n')
    reportFile.write('\nMean: ' + str(pairMagDiff))
    reportFile.write('\nError: ' + str(pairMagDiffErr))
    reportFile.write('\n\n### WDS form:\n')
    wdsform = str(ds[0]['wds_identifier']) + ',' + 'Date of observation' + ',' +  str(pairMeanTheta) + ',' +  str(pairMeanThetaErr) + ',' +  str(pairMeanRho) + ',' +  str(pairMeanRhoErr) + ',' +  'nan' + ',' +  'nan' + ',' +  str(pairMagDiff) + ',' +  str(pairMagDiffErr) + ',' + 'Filter wawelenght' + ',' + 'filter FWHM' + ',' + '0.2' + ',' + '1' + ',' + 'TAL_2022' + ',' +  'C' + ',' + '7'
    #print(str(wdsform))
    reportFile.write(str(wdsform))
    
#print(reportTable)
dsTable_by_object.write('double_stars.txt', format='ascii', overwrite=True, delimiter=',')
reportTable.write('double_stars_wds_format.txt', format='ascii', overwrite=True, delimiter=',')