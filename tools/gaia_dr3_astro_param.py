#Search for double star data in Gaia DR3 Astrophysical parameters
#Importing modules
import csv
import sys
import numpy as np
import datetime

print('Program start: ', datetime.datetime.now())

# Import csv file of Gaia DR3 to a numpy array
gaiaFileName = open(sys.argv[1], 'r')
gaiaStars = np.genfromtxt(gaiaFileName, delimiter=",")
print('Gaia Stars: ', gaiaStars)
#print(gaiaStars.dtype)
print('Gaia file import done: ', datetime.datetime.now())

# Read file contains double stars to a numpy array
doubleStarFileName = open(sys.argv[2], 'r')
doubleStars = np.genfromtxt(doubleStarFileName, delimiter=",")
print('Double Stars: ', doubleStars)
#print(doubleStars.dtype)
print('DS file import done: ', datetime.datetime.now())

# Check if any of the double stars are in the Gaia DR3 file based on the Source ID, collect findings into an array
doubleStarsAstroParam = np.empty((0, 226), float)

for ds in doubleStars:
    doubleStar = np.where(gaiaStars[:, 1] == ds)
    #print(gaiaStars[doubleStar])
    doubleStarsAstroParam = np.append(doubleStarsAstroParam, gaiaStars[doubleStar], axis=0)


print('The complete array:')
#doubleStarsAstroParam[:, 0:2] = doubleStarsAstroParam[:, 0:2].astype(int)
print(doubleStarsAstroParam)

# Write data to file
resultsHeader = ('solution_id,source_id,classprob_dsc_combmod_quasar,classprob_dsc_combmod_galaxy,classprob_dsc_combmod_star,classprob_dsc_combmod_whitedwarf,classprob_dsc_combmod_binarystar,classprob_dsc_specmod_quasar,classprob_dsc_specmod_galaxy,classprob_dsc_specmod_star,classprob_dsc_specmod_whitedwarf,classprob_dsc_specmod_binarystar,classprob_dsc_allosmod_quasar,classprob_dsc_allosmod_galaxy,classprob_dsc_allosmod_star,teff_gspphot,teff_gspphot_lower,teff_gspphot_upper,logg_gspphot,logg_gspphot_lower,logg_gspphot_upper,mh_gspphot,mh_gspphot_lower,mh_gspphot_upper,distance_gspphot,distance_gspphot_lower,distance_gspphot_upper,azero_gspphot,azero_gspphot_lower,azero_gspphot_upper,ag_gspphot,ag_gspphot_lower,ag_gspphot_upper,abp_gspphot,abp_gspphot_lower,abp_gspphot_upper,arp_gspphot,arp_gspphot_lower,arp_gspphot_upper,ebpminrp_gspphot,ebpminrp_gspphot_lower,ebpminrp_gspphot_upper,mg_gspphot,mg_gspphot_lower,mg_gspphot_upper,radius_gspphot,radius_gspphot_lower,radius_gspphot_upper,logposterior_gspphot,mcmcaccept_gspphot,libname_gspphot,teff_gspspec,teff_gspspec_lower,teff_gspspec_upper,logg_gspspec,logg_gspspec_lower,logg_gspspec_upper,mh_gspspec,mh_gspspec_lower,mh_gspspec_upper,alphafe_gspspec,alphafe_gspspec_lower,alphafe_gspspec_upper,fem_gspspec,fem_gspspec_lower,fem_gspspec_upper,fem_gspspec_nlines,fem_gspspec_linescatter,sife_gspspec,sife_gspspec_lower,sife_gspspec_upper,sife_gspspec_nlines,sife_gspspec_linescatter,cafe_gspspec,cafe_gspspec_lower,cafe_gspspec_upper,cafe_gspspec_nlines,cafe_gspspec_linescatter,tife_gspspec,tife_gspspec_lower,tife_gspspec_upper,tife_gspspec_nlines,tife_gspspec_linescatter,mgfe_gspspec,mgfe_gspspec_lower,mgfe_gspspec_upper,mgfe_gspspec_nlines,mgfe_gspspec_linescatter,ndfe_gspspec,ndfe_gspspec_lower,ndfe_gspspec_upper,ndfe_gspspec_nlines,ndfe_gspspec_linescatter,feiim_gspspec,feiim_gspspec_lower,feiim_gspspec_upper,feiim_gspspec_nlines,feiim_gspspec_linescatter,sfe_gspspec,sfe_gspspec_lower,sfe_gspspec_upper,sfe_gspspec_nlines,sfe_gspspec_linescatter,zrfe_gspspec,zrfe_gspspec_lower,zrfe_gspspec_upper,zrfe_gspspec_nlines,zrfe_gspspec_linescatter,nfe_gspspec,nfe_gspspec_lower,nfe_gspspec_upper,nfe_gspspec_nlines,nfe_gspspec_linescatter,crfe_gspspec,crfe_gspspec_lower,crfe_gspspec_upper,crfe_gspspec_nlines,crfe_gspspec_linescatter,cefe_gspspec,cefe_gspspec_lower,cefe_gspspec_upper,cefe_gspspec_nlines,cefe_gspspec_linescatter,nife_gspspec,nife_gspspec_lower,nife_gspspec_upper,nife_gspspec_nlines,nife_gspspec_linescatter,cn0ew_gspspec,cn0ew_gspspec_uncertainty,cn0_gspspec_centralline,cn0_gspspec_width,dib_gspspec_lambda,dib_gspspec_lambda_uncertainty,dibew_gspspec,dibew_gspspec_uncertainty,dibewnoise_gspspec_uncertainty,dibp0_gspspec,dibp2_gspspec,dibp2_gspspec_uncertainty,dibqf_gspspec,flags_gspspec,logchisq_gspspec,ew_espels_halpha,ew_espels_halpha_uncertainty,ew_espels_halpha_flag,ew_espels_halpha_model,classlabel_espels,classlabel_espels_flag,classprob_espels_wcstar,classprob_espels_wnstar,classprob_espels_bestar,classprob_espels_ttauristar,classprob_espels_herbigstar,classprob_espels_dmestar,classprob_espels_pne,azero_esphs,azero_esphs_uncertainty,ag_esphs,ag_esphs_uncertainty,ebpminrp_esphs,ebpminrp_esphs_uncertainty,teff_esphs,teff_esphs_uncertainty,logg_esphs,logg_esphs_uncertainty,vsini_esphs,vsini_esphs_uncertainty,flags_esphs,spectraltype_esphs,activityindex_espcs,activityindex_espcs_uncertainty,activityindex_espcs_input,teff_espucd,teff_espucd_uncertainty,flags_espucd,radius_flame,radius_flame_lower,radius_flame_upper,lum_flame,lum_flame_lower,lum_flame_upper,mass_flame,mass_flame_lower,mass_flame_upper,age_flame,age_flame_lower,age_flame_upper,flags_flame,evolstage_flame,gravredshift_flame,gravredshift_flame_lower,gravredshift_flame_upper,bc_flame,mh_msc,mh_msc_upper,mh_msc_lower,azero_msc,azero_msc_upper,azero_msc_lower,distance_msc,distance_msc_upper,distance_msc_lower,teff_msc1,teff_msc1_upper,teff_msc1_lower,teff_msc2,teff_msc2_upper,teff_msc2_lower,logg_msc1,logg_msc1_upper,logg_msc1_lower,logg_msc2,logg_msc2_upper,logg_msc2_lower,ag_msc,ag_msc_upper,ag_msc_lower,logposterior_msc,mcmcaccept_msc,mcmcdrift_msc,flags_msc,neuron_oa_id,neuron_oa_dist,neuron_oa_dist_percentile_rank,flags_oa')
np.savetxt("result.csv", doubleStarsAstroParam, delimiter = ",", header=resultsHeader, fmt='%s')
print('DS data process done: ', datetime.datetime.now())