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

# 1. Read the list of sources extracted from an image (fits) file
fitsFileName = sys.argv[1]
hdu = fits.open(fitsFileName)
mywcs = WCS(hdu[0].header)

# Estimate the background and background noise
data = hdu[0].data
mean, median, std = sigma_clipped_stats(data, sigma=5.0)  

daofind = DAOStarFinder(fwhm=10.0, threshold=18.0*std)  
sources = daofind(data - median)
#print(sources)

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

#print(segments)
#print(sources)

# Define Qtable for results
dr3Designation = np.array([], dtype=str)
dr3Ra = np.array([], dtype=np.float64)
dr3Dec = np.array([], dtype=np.float64)
dr3Parallax = np.array([], dtype=np.float64)
dr3ParallaxError = np.array([], dtype=np.float64)
dr3PmRa = np.array([], dtype=np.float64)
dr3PmDec = np.array([], dtype=np.float64)
imageId = np.array([], dtype=np.int32)
sourceId = np.array([], dtype=np.int64)
gMag = np.array([], dtype=np.float64)
sourceRa = np.array([], dtype=np.float64)
sourceDec = np.array([], dtype=np.float64)
sourceMag = np.array([], dtype=np.float64)
sourceTable = QTable([sourceId, dr3Designation, dr3Ra, dr3Dec, dr3Parallax, dr3ParallaxError, dr3PmRa, dr3PmDec, imageId, gMag, sourceRa, sourceDec, sourceMag], names=('source_id', 'designation', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'image_id', 'phot_g_mean_mag', 'source_ra', 'source_dec', 'source_mag'), meta={'name': 'first table'})

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
        #print("\n", star)
        #print('Separation based on catalog function: ', Angle(d2d))
        #print('Separation based on SkyCoord calculation: ', sep.arcsecond)
        #print('Image data: ', ra2, dec2)
        #print('Gaia data : ', gaiaStars[idx + 1][5], gaiaStars[idx + 1][7], str(gaiaStars[idx + 1][2]))
        sourceTable.add_row([gaiaStars[idx + 1][2], 'Gaia DR3 ' + str(int(gaiaStars[idx + 1][2])), gaiaStars[idx + 1][5], gaiaStars[idx + 1][7], gaiaStars[idx + 1][9], gaiaStars[idx + 1][10], gaiaStars[idx + 1][12], gaiaStars[idx + 1][14], star['id'], gaiaStars[idx + 1][69], ra2, dec2, star['mag']])

# Write found sources into file
tableFileName = (str(sys.argv[1][:-4] + '.csv'))
sourceTable.write(tableFileName, format='ascii', overwrite=True, delimiter=',')
print(sourceTable)






























    # dir_path = os.getcwd()
"""     print(star['designation'], dir_path, segmentName)
    isFile = os.path.isfile(f"{dir_path}/{segmentName}")
    print(isFile)
    if isFile == True:
        with open(f"{segmentName}", 'a') as f:
            mywriter = csv.DictWriter(f, delimiter=',', fieldnames = ['solution_id', 'designation', 'source_id', 'random_index', 'ref_epoch', 'ra', 'ra_error', 'dec', 'dec_error', 'parallax', 'parallax_error', 'parallax_over_error', 'pm', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 'ra_dec_corr', 'ra_parallax_corr', 'ra_pmra_corr', 'ra_pmdec_corr', 'dec_parallax_corr', 'dec_pmra_corr', 'dec_pmdec_corr', 'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr', 'astrometric_n_obs_al', 'astrometric_n_obs_ac', 'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_params_solved', 'astrometric_primary_flag', 'nu_eff_used_in_astrometry', 'pseudocolour', 'pseudocolour_error', 'ra_pseudocolour_corr', 'dec_pseudocolour_corr', 'parallax_pseudocolour_corr', 'pmra_pseudocolour_corr', 'pmdec_pseudocolour_corr', 'astrometric_matched_transits', 'visibility_periods_used', 'astrometric_sigma5d_max', 'matched_transits', 'new_matched_transits', 'matched_transits_removed', 'ipd_gof_harmonic_amplitude', 'ipd_gof_harmonic_phase', 'ipd_frac_multi_peak', 'ipd_frac_odd_win', 'ruwe', 'scan_direction_strength_k1', 'scan_direction_strength_k2', 'scan_direction_strength_k3', 'scan_direction_strength_k4', 'scan_direction_mean_k1', 'scan_direction_mean_k2', 'scan_direction_mean_k3', 'scan_direction_mean_k4', 'duplicated_source', 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag', 'phot_bp_n_obs', 'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag', 'phot_rp_n_obs', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error', 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag', 'phot_bp_rp_excess_factor', 'phot_bp_n_contaminated_transits', 'phot_bp_n_blended_transits', 'phot_rp_n_contaminated_transits', 'phot_rp_n_blended_transits', 'phot_proc_mode', 'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error', 'rv_method_used', 'rv_nb_transits', 'rv_nb_deblended_transits', 'rv_visibility_periods_used', 'rv_expected_sig_to_noise', 'rv_renormalised_gof', 'rv_chisq_pvalue', 'rv_time_duration', 'rv_amplitude_robust', 'rv_template_teff', 'rv_template_logg', 'rv_template_fe_h', 'rv_atm_param_origin', 'vbroad', 'vbroad_error', 'vbroad_nb_transits', 'grvs_mag', 'grvs_mag_error', 'grvs_mag_nb_transits', 'rvs_spec_sig_to_noise', 'phot_variable_flag', 'l', 'b', 'ecl_lon', 'ecl_lat', 'in_qso_candidates', 'in_galaxy_candidates', 'non_single_star', 'has_xp_continuous', 'has_xp_sampled', 'has_rvs', 'has_epoch_photometry', 'has_epoch_rv', 'has_mcmc_gspphot', 'has_mcmc_msc', 'in_andromeda_survey', 'classprob_dsc_combmod_quasar', 'classprob_dsc_combmod_galaxy', 'classprob_dsc_combmod_star', 'teff_gspphot', 'teff_gspphot_lower', 'teff_gspphot_upper', 'logg_gspphot', 'logg_gspphot_lower', 'logg_gspphot_upper', 'mh_gspphot', 'mh_gspphot_lower', 'mh_gspphot_upper', 'distance_gspphot', 'distance_gspphot_lower', 'distance_gspphot_upper', 'azero_gspphot', 'azero_gspphot_lower', 'azero_gspphot_upper', 'ag_gspphot', 'ag_gspphot_lower', 'ag_gspphot_upper', 'ebpminrp_gspphot', 'ebpminrp_gspphot_lower', 'ebpminrp_gspphot_upper', 'libname_gspphot'])
            mywriter.writerow(star) # Add selected rows to file
            print('file exists')
            # Close the file object
            f.close()
    else:
        with open(f"{segmentName}", 'w') as f:
            mywriter = csv.DictWriter(f, delimiter=',', fieldnames = ['solution_id', 'designation', 'source_id', 'random_index', 'ref_epoch', 'ra', 'ra_error', 'dec', 'dec_error', 'parallax', 'parallax_error', 'parallax_over_error', 'pm', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 'ra_dec_corr', 'ra_parallax_corr', 'ra_pmra_corr', 'ra_pmdec_corr', 'dec_parallax_corr', 'dec_pmra_corr', 'dec_pmdec_corr', 'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr', 'astrometric_n_obs_al', 'astrometric_n_obs_ac', 'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_params_solved', 'astrometric_primary_flag', 'nu_eff_used_in_astrometry', 'pseudocolour', 'pseudocolour_error', 'ra_pseudocolour_corr', 'dec_pseudocolour_corr', 'parallax_pseudocolour_corr', 'pmra_pseudocolour_corr', 'pmdec_pseudocolour_corr', 'astrometric_matched_transits', 'visibility_periods_used', 'astrometric_sigma5d_max', 'matched_transits', 'new_matched_transits', 'matched_transits_removed', 'ipd_gof_harmonic_amplitude', 'ipd_gof_harmonic_phase', 'ipd_frac_multi_peak', 'ipd_frac_odd_win', 'ruwe', 'scan_direction_strength_k1', 'scan_direction_strength_k2', 'scan_direction_strength_k3', 'scan_direction_strength_k4', 'scan_direction_mean_k1', 'scan_direction_mean_k2', 'scan_direction_mean_k3', 'scan_direction_mean_k4', 'duplicated_source', 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag', 'phot_bp_n_obs', 'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag', 'phot_rp_n_obs', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error', 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag', 'phot_bp_rp_excess_factor', 'phot_bp_n_contaminated_transits', 'phot_bp_n_blended_transits', 'phot_rp_n_contaminated_transits', 'phot_rp_n_blended_transits', 'phot_proc_mode', 'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error', 'rv_method_used', 'rv_nb_transits', 'rv_nb_deblended_transits', 'rv_visibility_periods_used', 'rv_expected_sig_to_noise', 'rv_renormalised_gof', 'rv_chisq_pvalue', 'rv_time_duration', 'rv_amplitude_robust', 'rv_template_teff', 'rv_template_logg', 'rv_template_fe_h', 'rv_atm_param_origin', 'vbroad', 'vbroad_error', 'vbroad_nb_transits', 'grvs_mag', 'grvs_mag_error', 'grvs_mag_nb_transits', 'rvs_spec_sig_to_noise', 'phot_variable_flag', 'l', 'b', 'ecl_lon', 'ecl_lat', 'in_qso_candidates', 'in_galaxy_candidates', 'non_single_star', 'has_xp_continuous', 'has_xp_sampled', 'has_rvs', 'has_epoch_photometry', 'has_epoch_rv', 'has_mcmc_gspphot', 'has_mcmc_msc', 'in_andromeda_survey', 'classprob_dsc_combmod_quasar', 'classprob_dsc_combmod_galaxy', 'classprob_dsc_combmod_star', 'teff_gspphot', 'teff_gspphot_lower', 'teff_gspphot_upper', 'logg_gspphot', 'logg_gspphot_lower', 'logg_gspphot_upper', 'mh_gspphot', 'mh_gspphot_lower', 'mh_gspphot_upper', 'distance_gspphot', 'distance_gspphot_lower', 'distance_gspphot_upper', 'azero_gspphot', 'azero_gspphot_lower', 'azero_gspphot_upper', 'ag_gspphot', 'ag_gspphot_lower', 'ag_gspphot_upper', 'ebpminrp_gspphot', 'ebpminrp_gspphot_lower', 'ebpminrp_gspphot_upper', 'libname_gspphot'])
            mywriter.writeheader()    #Add header
            mywriter.writerow(star) # Add selected rows to file
            print('no file exists, creating')
            # Close the file object
            f.close() """

# 3. Search sources in the catalog








# newFileName = str(filename[-28:])

""" for line in file:
    starList.append(line) """

# 4. Write matching data into file.

# printList = []




            # Send the final star list to an array
            # printList.append(star)
        #print(line['parallax'], line['parallax_error'])

""" with open('result.csv', 'w') as f:
    mywriter = csv.DictWriter(f, delimiter=',', fieldnames = ['solution_id', 'designation', 'source_id', 'random_index', 'ref_epoch', 'ra', 'ra_error', 'dec', 'dec_error', 'parallax', 'parallax_error', 'parallax_over_error', 'pm', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 'ra_dec_corr', 'ra_parallax_corr', 'ra_pmra_corr', 'ra_pmdec_corr', 'dec_parallax_corr', 'dec_pmra_corr', 'dec_pmdec_corr', 'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr', 'astrometric_n_obs_al', 'astrometric_n_obs_ac', 'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_params_solved', 'astrometric_primary_flag', 'nu_eff_used_in_astrometry', 'pseudocolour', 'pseudocolour_error', 'ra_pseudocolour_corr', 'dec_pseudocolour_corr', 'parallax_pseudocolour_corr', 'pmra_pseudocolour_corr', 'pmdec_pseudocolour_corr', 'astrometric_matched_transits', 'visibility_periods_used', 'astrometric_sigma5d_max', 'matched_transits', 'new_matched_transits', 'matched_transits_removed', 'ipd_gof_harmonic_amplitude', 'ipd_gof_harmonic_phase', 'ipd_frac_multi_peak', 'ipd_frac_odd_win', 'ruwe', 'scan_direction_strength_k1', 'scan_direction_strength_k2', 'scan_direction_strength_k3', 'scan_direction_strength_k4', 'scan_direction_mean_k1', 'scan_direction_mean_k2', 'scan_direction_mean_k3', 'scan_direction_mean_k4', 'duplicated_source', 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag', 'phot_bp_n_obs', 'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag', 'phot_rp_n_obs', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error', 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag', 'phot_bp_rp_excess_factor', 'phot_bp_n_contaminated_transits', 'phot_bp_n_blended_transits', 'phot_rp_n_contaminated_transits', 'phot_rp_n_blended_transits', 'phot_proc_mode', 'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error', 'rv_method_used', 'rv_nb_transits', 'rv_nb_deblended_transits', 'rv_visibility_periods_used', 'rv_expected_sig_to_noise', 'rv_renormalised_gof', 'rv_chisq_pvalue', 'rv_time_duration', 'rv_amplitude_robust', 'rv_template_teff', 'rv_template_logg', 'rv_template_fe_h', 'rv_atm_param_origin', 'vbroad', 'vbroad_error', 'vbroad_nb_transits', 'grvs_mag', 'grvs_mag_error', 'grvs_mag_nb_transits', 'rvs_spec_sig_to_noise', 'phot_variable_flag', 'l', 'b', 'ecl_lon', 'ecl_lat', 'in_qso_candidates', 'in_galaxy_candidates', 'non_single_star', 'has_xp_continuous', 'has_xp_sampled', 'has_rvs', 'has_epoch_photometry', 'has_epoch_rv', 'has_mcmc_gspphot', 'has_mcmc_msc', 'in_andromeda_survey', 'classprob_dsc_combmod_quasar', 'classprob_dsc_combmod_galaxy', 'classprob_dsc_combmod_star', 'teff_gspphot', 'teff_gspphot_lower', 'teff_gspphot_upper', 'logg_gspphot', 'logg_gspphot_lower', 'logg_gspphot_upper', 'mh_gspphot', 'mh_gspphot_lower', 'mh_gspphot_upper', 'distance_gspphot', 'distance_gspphot_lower', 'distance_gspphot_upper', 'azero_gspphot', 'azero_gspphot_lower', 'azero_gspphot_upper', 'ag_gspphot', 'ag_gspphot_lower', 'ag_gspphot_upper', 'ebpminrp_gspphot', 'ebpminrp_gspphot_lower', 'ebpminrp_gspphot_upper', 'libname_gspphot'])
    mywriter.writeheader()    #Add header
    mywriter.writerows(printList) # Add selected rows to file

    # Close the file object
    f.close() """

# print(starList)

# NOTES:
# Search for coordinates in a catalog: https://keflavich-astropy.readthedocs.io/en/latest/coordinates/matchsep.html

## Example