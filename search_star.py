# importing the module
import sys
import csv
import numpy as np

gaiaData = open(sys.argv[1], 'r')
file1 = csv.DictReader(gaiaData)
fieldnames = ['designation', ['ra'], ['dec'], ['parallax'], ['parallax_error'], ['pmra'], ['pmdec'], ['phot_g_mean_mag']]

designationList = open(sys.argv[2], 'r')
file2 = csv.DictReader(designationList)
fieldnames = [['designation']]

starList = list(file1)
gaiaList = list(file2)

#print(starList)
#print(gaiaList)

for gstar in gaiaList:
    StarA = (gstar['designation'])
    for lstar in starList:
        StarB = (lstar['designation'])
        if StarA == StarB:
            print(StarA, ',', lstar['ra'], ',', lstar['dec'], ',', lstar['parallax'], ',', lstar['parallax_error'], ',', lstar['pmra'], ',', lstar['pmdec'], ',', lstar['phot_g_mean_mag'], lstar['phot_bp_mean_mag'], ',', lstar['phot_rp_mean_mag'], ',', lstar['radial_velocity'], ',', lstar['radial_velocity_error'], ',', lstar['non_single_star'], ',', lstar['teff_gspphot'], ',', lstar['teff_gspphot_lower'], ',', lstar['teff_gspphot_upper'], ',', lstar['logg_gspphot'], ',', lstar['logg_gspphot_lower'], ',', lstar['logg_gspphot_upper'], ',', lstar['distance_gspphot'], ',', lstar['distance_gspphot_lower'], ',', lstar['distance_gspphot_upper'])

#starIndex = ()
#for gstar in gaiaList:
#    starIndex = starList.index(gstar['designation'])
#    print(starList([starIndex]))

#while gaiaList['designation']:
#	filterObject = filter(lambda a: gaiaList in a, starList)
#print(filter_object)

#for gstar in gaiaList:
#     StarA = (gstar['designation'])
#     filterObject = filter(lambda a: StarA in a, starList)
#     print(filterObject)#
