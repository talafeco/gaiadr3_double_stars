#Search for double star data in Gaia DR3 Astrophysical parameters
#Importing modules
import csv
import sys
import numpy as np

# Import csv file of Gaia DR3 to a numpy array
gaiaFileName = open(sys.argv[1], 'r')
gaiaStars = np.genfromtxt(gaiaFileName, delimiter=",")
print('Gaia Stars: ', gaiaStars)
print(gaiaStars.dtype)

# Read file contains double stars to a numpy array
doubleStarFileName = open(sys.argv[2], 'r')
doubleStars = np.genfromtxt(doubleStarFileName, delimiter=",")
print('Double Stars: ', doubleStars)
print(doubleStars.dtype)

# Check if any of the double stars are in the Gaia DR3 file based on the Source ID, collect findings into an array

for ds in doubleStars:
    #print(ds)
    doubleStar = np.where(gaiaStars == ds)
    print(doubleStar)
# Write findings into a csv file