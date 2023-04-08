# Gaia DR3 Double Stars
Search visual double star candidates in Gaia DR3. Mag &lt; 15, par > 0.5.

## Fits double star measurement
(fits_double_measurement.py)
1. Read image file, apply dark+flat+bias
2. Run astrometry on image files
3. Search sources on the files
4. Find pairs of interest, record ra, dec positions
5. Calculate sep and pa based on positions for the selected sources
6. Write data into ascii text file
7. Create plotted image file about the pair

# Changelog 2023-04-07
(wdsreport.py)
1. All component are measured on images
2. Mean and error values are calculated for each component measurements: theta, rho, magnitude difference

## TBD
1. Clean up code
2. Extend data table with additional useful columns
3. Extend reports: list of individual measurements
    - Added, to be formated
4. Add Gaia DR3 identification of the components
    - Done
5. Calculate Physical properties of the double star: mass, absolute magnitude, luminosity, gravitational bound
    - calculations added, report generation tbd
6. Create image of each components, comments on the image is optional
7. Create HRD diagram of the components
8. Update report file path to be created in the image directory
9. Fix magnitude error calculation algorithm
10. Add magnitude difference to the updated table
11. Extend star mass calculation with alternate equatione based on measured temperature