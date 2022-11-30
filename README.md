# Gaia DR3 Double Stars
Search visual double star candidates in Gaia DR3. Mag &lt; 15, par > 0.5.

## Fits double star measurement
(fits_double_measurement.py)
1. read image file, apply dark+flat+bias
2. run astrometry on image files
3. search sources on the files
4. find pairs of interest, record ra, dec positions
5. calculate sep and pa based on positions for the selected sources
6. write data into ascii text file
7. create plotted image file about the pair