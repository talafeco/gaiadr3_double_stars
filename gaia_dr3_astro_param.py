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
print(gaiaStars.dtype)
print('Gaia file import done: ', datetime.datetime.now())

# Read file contains double stars to a numpy array
doubleStarFileName = open(sys.argv[2], 'r')
doubleStars = np.genfromtxt(doubleStarFileName, delimiter=",")
print('Double Stars: ', doubleStars)
print(doubleStars.dtype)
print('DS file import done: ', datetime.datetime.now())

# Check if any of the double stars are in the Gaia DR3 file based on the Source ID, collect findings into an array
doubleStarsAstroParam = []
for ds in doubleStars:
    #print(ds)
    doubleStar = np.where(gaiaStars[:, 1] == ds)
    print(gaiaStars[doubleStar])
    #doubleStarsAstroParam = np.append(doubleStars, ds[doubleStar])
print(doubleStarsAstroParam)

print('DS data process done: ', datetime.datetime.now())
# Write findings into a csv file