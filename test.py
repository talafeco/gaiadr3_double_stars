import numpy as np
import math

magnitude = 7.523838
parallax = 2.655862264105706

def calcDistance(par):
    dist = 1 / (math.fabs(par/1000))
    return dist

def calcAbsMag(gmag, par):
    dist = calcDistance(par)
    absmag = gmag - 5 * math.log(dist, 10) + 5
    return absmag

def roundNumber(num):
    if type(num) == float or type(num) == int:
        finalNumber = round(num, 3)
    else:
        finalNumber = np.nan
    return finalNumber


print(str(roundNumber(calcAbsMag(magnitude, parallax))))