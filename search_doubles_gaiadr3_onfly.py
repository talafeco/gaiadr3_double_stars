# importing the module
import sys
import csv
import math

# Function to calculate Delta RA
def deltaRa(raa, rab, decb):
    deltara = (rab - raa) * math.cos(math.radians(decb))
    return deltara

# Function to calculate Delta DEC
def deltaDec(decb, deca):
    deltadec = decb - deca
    return deltadec

# Function to calculate Theta (position angle)
def thetaCalc(deltara,deltadec):
    thetacalc = math.degrees(math.atan(deltara/deltadec))
    return thetacalc

# Function to calculate Rho (separation)
def rhoCalc(raa, deca, rab, decb):
    rhocalc = math.sqrt(((raa-rab) * math.cos(math.radians(deca))) ** 2 + (deca - decb) ** 2) * 3600
    return rhocalc

# Function to calculate the distance of the star based on the parallax
def calcDistance(par):
    dist = 1 / (par/1000)
    return dist

# Function to calculate the minimum distance of the star based on the parallax error
def calcDistanceMin(par, err):
    distmin = 1 / ((par + err) / 1000)
    return distmin

# Function to calculate the maximum distance of the star based on the parallax error
def calcDistanceMax(par, err):
    distmax = 1 / ((par - err) / 1000)
    return distmax

# Function to calculate the parallax factor of the pair
# EXCEL formula =1-ABS((raA-raB)/(0,5*(raA+raB)))
def calcParallaxFactor(para, parb):
    parfac = 1 - math.fabs((para-parb)/(0.5*(para+parb)))
    return parfac

# Function to calculate the proper motion factor of the stars and define if it is a CPM (common proper motion) pair
# EXCEL formula =ABS(1-((SQRT(pmraA-pmraB)^2+(pmdecA-pmdecB)^2)/(SQRT(pmraA^2+pmdecA^2)+(pmraB^2+pmdecB^2)))))
def calcPmFactor(pmraa, pmdeca, pmrab, pmdecb):
    pmfac = math.fabs(1-((math.sqrt(((pmraa-pmrab) ** 2) + ((pmdeca-pmdecb) ** 2))/(math.sqrt((pmraa ** 2) + (pmdeca ** 2))+((pmrab ** 2) + pmdecb ** 2)))))
    return pmfac

# Richard W. Harshaw (2016)
# CCD Measurements of 141 Proper Motion Stars:
# The Autumn 2015 Observing Program at the
# Brilliant Sky Observatory, Part 3
# The proper motion of a star can be depicted as a
# vector. When the resultant of the two vectors is divided
# by the largest vector, the result will either be zero (or
# very near it) if the proper motions are identical,
# somewhere between 20% and 60% of the resultant of the
# vectors, or over 60% of the resultant. Pairs in the first
# category are classed as Common Proper Motion pairs,
# or CPM. Pairs in the second category are classed as
# Similar Proper Motion pairs (SPM), and those in the
# third category are classed as Different Proper Motion
# pairs (DPM).

filename = open(sys.argv[1], 'r')
file = csv.DictReader(filename)
fieldnames = ['designation', 'ra', 'dec', 'parallax', 'parallax_error']

# here a logic to be inserted based on modulo (a%b tp get the remains of the division) to decide, which array should be populated with the star
# modulo: (deg(ra or dec)%5)//1+1

# here an array definition needs to be inserted for the matrix
starList = []
for line in file:
    starList.append(line)

# Creating empty arrays for Star related calculations

StarA = []
StarB = []

# Print header for the table
# print('Star A name | Star B name | Theta  | Rho | Star A magnitude | Star B magnitude | Star A distance (pc) | Star A max distance | Star A min distance | Star A distance range | Star B distance (pc) | Star B max distance | Star B min distance | Star B distance range | distanceCommon | starParallaxFactor | starPmFactor | Proper Motion indicator')

# here a for loop to be inserted to go through the arrays somehow :)
for star in starList: ## modify according to arrays instead of starlist
    StarA = (star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'], star['phot_g_mean_mag'], star['source_ra'], star['source_dec'])
    for star in starList:
        StarB = (star['designation'], star['ra'], star['dec'], star['parallax'], star['parallax_error'], star['pmra'], star['pmdec'], star['phot_g_mean_mag'], star['source_ra'], star['source_dec'])
        if StarA != StarB and float(StarA[7]) < float(StarB[7]) and float(StarA[3]) != 0 and float(StarB[3]) != 0:
            #Set input data
            #print('StarA', StarA)
            #print('StarB', StarB)
            starRa1 = float(StarA[1])
            starDec1 = float(StarA[2])
            starRa2 = float(StarB[1])
            starDec2 = float(StarB[2])
            starParallax1 = float(StarA[3])
            starParallaxError1 = float(StarA[4])
                        
            # Calculate the widest possible separation for StarA
            possSep1 = 10000 / calcDistanceMax(starParallax1, starParallaxError1)
            rhoStar = rhoCalc(starRa1, starDec1, starRa2, starDec2)
            #print(possSep1, rhoStar)
            if possSep1 > rhoStar:
                starName1 = StarA[0]
                starName2 = StarB[0]
                starParallax2 = float(StarB[3])
                starParallaxError2 = float(StarB[4])
                starPmRa1 = float(StarA[5])
                starPmDec1 = float(StarA[6])
                starPmRa2 = float(StarB[5])
                starPmDec2 = float(StarB[6])
                starGMag1 = float(StarA[7])
                starGMag2 = float(StarB[7])
                starActualRa1 = float(StarA[8])
                starActualDec1 = float(StarA[9])
                starActualRa2 = float(StarB[8])
                starActualDec2 = float(StarB[9])

                # Value to modify Theta according to the appropriate quadrant
                addThetaValue = ()
                if deltaRa(starRa1, starRa2, starDec2) > 0 and deltaDec(starDec2, starDec1) > 0:
                    addThetaValue = 0
                elif deltaRa(starRa1, starRa2, starDec2) > 0 and deltaDec(starDec2, starDec1) < 0:
                    addThetaValue = 180
                elif deltaRa(starRa1, starRa2, starDec2) < 0 and deltaDec(starDec2, starDec1) < 0:
                    addThetaValue = 180
                elif deltaRa(starRa1, starRa2, starDec2) < 0 and deltaDec(starDec2, starDec1) > 0:
                    addThetaValue = 360
                
                # Calculate actual data based on functions
                thetaStar = thetaCalc(deltaRa(starRa1, starRa2, starDec2), deltaDec(starDec2, starDec1)) + addThetaValue
                thetaActual = thetaCalc(deltaRa(starActualRa1, starActualRa2, starActualDec2), deltaDec(starActualDec2, starActualDec1)) + addThetaValue
                rhoActual = rhoCalc(starActualRa1, starActualDec1, starActualRa2, starActualDec2)
                # rhoStar = rhoCalc(starRa1, starDec1, starRa2, starDec2) # Already calculated
                starDistance1 = calcDistance(starParallax1)
                starDistanceMax1 = calcDistanceMax(starParallax1, starParallaxError1)
                starDistanceMin1 = calcDistanceMin(starParallax1, starParallaxError1)
                starDistanceRange1 = starDistanceMax1 - starDistanceMin1
                starDistance2 = calcDistance(starParallax2)
                starDistanceMax2 = calcDistanceMax(starParallax2, starParallaxError2)
                starDistanceMin2 = calcDistanceMin(starParallax2, starParallaxError2)
                starDistanceRange2 = starDistanceMax2 - starDistanceMin2
                starParallaxFactor = calcParallaxFactor(starParallax1, starParallax2)
                starPmFactor = calcPmFactor(starPmRa1, starPmDec1, starPmRa2, starPmDec2)
                
                # Check if stars shares a common distance range
                distanceCommon = ()
                if starDistanceMin1 < starDistanceMin2 < starDistanceMax1 or starDistanceMin2 < starDistanceMin1 < starDistanceMax2:
                    distanceCommon = 'overlapping'
                else:
                    distanceCommon = 'no'
                
                # Check if the pair is a Common Proper Motion pairs (CPM), Similar Proper Motion (SPM) or Different Proper Motion (DPM)
                pmCommon = ()
                if starPmFactor >= 0.8:
                    pmCommon = 'CPM'
                elif 0.4 <= starPmFactor < 0.8:
                    pmCommon = 'SPM'
                elif starPmFactor < 0.4:
                    pmCommon = 'DPM'
                
                #Print data, if stars are close and share a common distance range
                if distanceCommon == 'overlapping':
                    print(str(sys.argv[1]),'|',starName1,'|',starName2,'|',thetaStar,'|',rhoStar,'|',starGMag1,'|',starGMag2,'|',starDistance1,'|',starDistanceMax1,'|',starDistanceMin1,'|',starDistanceRange1,'|',starDistance2,'|',starDistanceMax2,'|',starDistanceMin2,'|',starDistanceRange2,'|',distanceCommon,'|',starParallaxFactor,'|',starPmFactor,'|',pmCommon,'|',thetaActual,'|',rhoActual)


# Iterating over each row and append values to empty list
""" for col in file:
    starDesig.append(col['designation'])
    starRaArray.append(col['ra'])
    starDecArray.append(col['dec'])
    starParalArray.append(col['parallax'])
    starParalErrArray.append(col['parallax_error']) """

# Printing arrays
""" print('Star designations:', len(starDesig))
print('Stars RA:', len(starRaArray))
print('Stars Dec:', len(starDecArray))
print('Stars parallax:', len(starParalArray))
print('Stars parallax error:', len(starParalErrArray)) """

""" class StarAttributes:
    def __init__(mystar, designation, ra, dec, parallax, parallax_error):
        mystar.designation = designation
        mystar.ra = ra
        mystar.dec = dec
        mystar.parallax = parallax
        mystar.parallax_error = parallax_error
    
    def calcDistance(mystar.parallax):
        dist = 1 / (par/1000)
        return dist

    def calcDistanceMin(mystar.parallax, mystar.parallax_error):
        distmin = 1 / ((mystar.parallax + mystar.parallax_error) / 1000)
        return distmin

    def calcDistanceMax(mystar.parallax, mystar.parallax_error):
        distmax = 1 / ((mystar.parallax - mystar.parallax_error) / 1000)
        return distmax """

