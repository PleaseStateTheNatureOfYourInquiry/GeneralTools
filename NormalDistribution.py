
'''
Three functions to calculate the theoretical cumulative normal distribution for a given average (mu) and standard deviation (sigma)
'''

import numpy as np


def getNormalDistributionValue (xi,mu,sigma):

    normalValue = np.exp ( -0.5 * (xi - mu) * (xi - mu) / sigma / sigma )

    return normalValue


def getNormalDistribution (mu,sigma):

    cnst = 1 / (sigma * np.sqrt (2 * np.pi))

    # x a list of the variable centred at its average and ranging from  -4 * sigma  to  +4 * sigma
    x = [ (mu + i * sigma / 10)  for i in range (-40, 40) ]

    normalDistribution = [ cnst * getNormalDistributionValue (x[i], mu, sigma)  for i in range(len(x))]

    return x, normalDistribution


def getCumulativeNormalDistribution (mu,sigma):

    x, normalDistribution = getNormalDistribution (mu,sigma)

    cumulativeNormalDistribution = []

    cumulativeNormalDistribution.append (normalDistribution[0])

    for i in range (1, len (x) - 1):

        cumulativeNormalDistribution.append ( cumulativeNormalDistribution [-1] +
                                              ( x [i] - x [i-1] ) * ( normalDistribution [i] + normalDistribution [i-1] ) / 2 )

    cumulativeNormalDistribution.append (cumulativeNormalDistribution[-1])

    return x, cumulativeNormalDistribution


