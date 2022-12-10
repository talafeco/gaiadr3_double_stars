# Search for double star data in Gaia DR3 Astrophysical parameters
# Importing modules
import csv
import sys
import numpy as np

# Import csv file of Gaia DR3 to a numpy array
gaiaFileName = open(sys.argv[1], 'r')
gaiaStars = np.genfromtxt(gaiaFileName, delimiter=",")

# Read file contains double stars to a numpy array
doubleStarFileName = open(sys.argv[1], 'r')
doubleStars = np.genfromtxt(doubleStarFileName, delimiter=",")

# Check if any of the double stars are in the Gaia DR3 file based on the Source ID, collect findings into an array
for ds in doubleStars:
    doubleStar = np.where(gaiaStars = ds)
    print(doubleStar)

# Write findings into a csv file