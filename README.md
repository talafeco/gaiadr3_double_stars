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

# Changelog 2023-04-08
(wdsreport.py)
1. All component are measured on images
    - Done
2. Mean and error values are calculated for each component measurements: theta, rho, magnitude difference
    - Done

## TBD
1. Clean up code
2. Extend data table with additional useful columns
3. Extend reports: list of individual measurements - Done
4. Add Gaia DR3 identification of the components - Done
5. Calculate Physical properties of the double star: mass, absolute magnitude, luminosity, gravitational bound - Done
6. Extend star mass calculation with alternate equation based on measured temperature
7. Create HRD diagram of the components - Done
8. Update report file path to be created in the image directory - Done  
9. Fix magnitude error calculation algorithm - Done
10. Add magnitude difference to the updated table - Done
11. Calculate Gaia star position for current epoch from 2016 + proper motion, calculate measurement difference
12. Add rp calculation
13. Query Tycho database for stars? (wdsreport, Astroquery)
14. Check pair magnitude difference (wds + gaia + measured), calculate difference
15. Calculate background, calculate star magnitudes with photometry
16. Remove spaces from filenames
17. Calclate sensor tilt
18. Calcuate fwhm of the image
19. Increase plot resolution and create colors for the stars, optionally color black the background


# New functions
## Create image of each components
1. Stack fits images, process dark, flat and bias images
2. Find center pixel of the double star
3. Crop image to a reasonable size of the double star
4. Draw lines towards components, write the component designation over the line
5. Mark north and east on the image (lower left corner)
6. Write the Double star name to the image (upper left corner)
7. Write observation date and time (upper right corner)
8. Write measurements of the components: designation, theta, rho (lower right corner)
9. Calculate the double star's precise coordinates for current epoch
10. Get all double stars from wds, which should be found on the images
11. Search double stars based on magnitude difference, separation and position angle, if not found based on the coordinates
12. Calculate double star orbit plane (inclination) for components, compare component planes.

## Function externsions
1. Gaia data collection in one line
    - Done
2. Calculate hrd placement based on temperature & luminosity
    - Done
3. Create HRD plot
    - Done
4. Calculate orbit period
5. Check, if absolute magnitude measurement is possible on the images
6. Calculate Gaia source positions to the current epoch and compare them to the measurement results.


# Zsolti's input
## Issues
1. Részben t toolhoz kikurkássza a Simbadból a hiányzókat. Pl. a Gaia nem adott semmit a csillagról, de a Hipparchos már egyszer megmérte a pm-t és a plx-t. Ha ez sem, akkor tényleg nincsen semmink róla. Tehát maradna az ismertelen eredmény.
3. Nem tudom megoldható-e, hogy egy többes rendszernél egy csoportba legyenek a komponensek? PL STF133 A; STF133 AC; ABH111 AD stb. A txt fájba összefűzheti ezeket?
4. Részben OK - Megoldható lenne-e, hogy ha több pár van egy képen, mert ez sűrűn előfordul, a CMD-ben elválassza őket mondjuk egy „-----„ sávval?
5. Megoldható-e, hogy a txt fájba ne ömlesztve jelenjenek meg a DR3 adatok? Nem sok kell belőle, de lehetne rendezni ezeket valamiféle kvázi táblázatba?
6. OK - Szerintem nem kell végtelen tizedesjelű szám. Tudom így pontos a dolog, de pl. a végső eredményt elég lenne 3 tizedesre kerekítve kiíratni vele.
7. Bevinni az rp kalkulációt
8. OK - fit és fits fileokra is fusson le
9. OK - radial velocity error nem jelenik meg a listában rendesen

# Tools
1. Convert date/time object to julian year
2. Run calculation functions manually based on coordinates / Gaia designations