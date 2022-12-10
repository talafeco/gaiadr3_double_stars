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
8. tbd