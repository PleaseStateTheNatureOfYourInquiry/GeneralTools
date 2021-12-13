from readtable import ReadTable

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import sys
import re

currentVersionHandyTools = '20211213'


class HandyTools:
    """
    HandyTools version 20211213

    Functions:
    filesInFolderTree
    getDateAndTimeForHeader
    getUncertaintyElectrogramADU
    QQPlot
    linearLeastSquare
    
    """


    def filesInFolderTree (startPath, extension = '.annotation', checkStartPathOnly = False):
        """
        Determine a list of strings that contain all the files paths/names ending on  extension  down the directory tree starting at  startPath .
         startPath  is the point in the directory tree where to start. The function then looks down the tree and seeks out all the files.
        """

        dirPath = Path (startPath)

        listOfFileNames = []
        if checkStartPathOnly:

            listOfAllFileNames = sorted ( os.listdir (startPath) )         
            if extension == '':
 
                # the file .DS_Store is a MAC OS administration file, not of interest to keep
                listOfFileNames = [ fileName  for fileName in listOfAllFileNames  if fileName != '.DS_Store']
           
            
            else:
                    
                listOfFileNames = [ fileName  for fileName in listOfAllFileNames if fileName.endswith (extension) ]

        
        else:

            if extension == '':

                # the file .DS_Store is a MAC OS administration file, not of interest to keep
                listOfFileNames = [ os.path.join (root, name)
                                    for root, dirs, files in os.walk (dirPath)
                                    for name in sorted (files)
                                    if name != '.DS_Store'
                ]
        
        
            else:

                listOfFileNames = [ os.path.join (root, name)
                                    for root, dirs, files in os.walk (dirPath)
                                    for name in sorted (files)
                                    if name.endswith (extension)
                ]


        return listOfFileNames



    def getDateAndTimeForHeader ():
        """
        Get a date and time string for printing in a file header to indicate the time of creation of the file
        """
        
        timeOfExport = datetime.datetime.now ()
        
        dateTimeString =  'on ' + str (timeOfExport.year) + '-' + str( timeOfExport.month) + '-' + str (timeOfExport.day) + ' at ' \
            + str (timeOfExport.hour) + 'h' + str (timeOfExport.minute) + 'm' + str (timeOfExport.second) +'s'
        
        return dateTimeString



    def getUncertaintyElectrogramADU (electrogramADU):
        """
        Determine the uncertainty level of the voltage measurements of an electrogram in ADU!!
        """

        numberOfSamples = len (electrogramADU)
        np.seterr (all='ignore')

        segments = []
        segments.append (electrogramADU [1] - electrogramADU [0])
        for iDataPoint in range (1, numberOfSamples-1):

            deltaVoltage = electrogramADU [iDataPoint+1] - electrogramADU [iDataPoint]

            if (segments [-1] > 0 and deltaVoltage >= 0) or (segments [-1] < 0 and deltaVoltage <= 0):

                segments [-1] += deltaVoltage

            else:

                segments.append (deltaVoltage)

        segments = np.array (segments)
        iSegmentsNoise = np.where (np.abs (segments) < 350)[0]

        uncertaintyLevelADU = np.std (segments [iSegmentsNoise])

        return uncertaintyLevelADU, segments



    def QQPlot (xInput, xlabelToPrint = 'input values', ylabelToPrint = 'z (sigma)',
                QQTitleToPrint = 'HandyTools version ' + currentVersionHandyTools + ': QQ-plot of input data',
                HistTitleToPrint = 'HandyTools version ' + currentVersionHandyTools + ': Histogram of input data', 
                plotTextAverageMedian = True):
        """
        Calculate and plot the QQ-plot (or Quantile-Quantile plot) and the histogram of a given list of input values
        """

        xInput = np.array (xInput)
        xInputMean = np.mean (xInput)
        xInputMedian = np.median (xInput)
        xInputStd = np.std (xInput)

        x, cumulativeNormalDistribution = getCumulativeNormalDistribution (xInputMean, xInputStd)

        xInputOrdered = xInput.copy ()
        xInputOrdered.sort ()

        xInputCumulativeNormalDistribution = [ (i-1) / len (xInputOrdered) for i in range (1,len (xInputOrdered)+1) ]

        xValues = []

        # Search for the indices in the true normally distributed cumulativeNormalDistribution between which
        #  the estimated cumulative distrbution value of the input (xInputCumulativeNormalDistribution) values lie
        #  and calculate the true normally distributed cumulative value that would correspond to the this iDataPoint
        #  by interpolation.
        for i, xInputCumulativeNormalDistributionValue in enumerate (xInputCumulativeNormalDistribution):

            cumulativeNormalDistributionAdded = cumulativeNormalDistribution.copy ()
            cumulativeNormalDistributionAdded.append (xInputCumulativeNormalDistributionValue)
            cumulativeNormalDistributionAdded.sort ()
            valueIndex = np.where ( np.array (cumulativeNormalDistributionAdded) == xInputCumulativeNormalDistributionValue )[0][-1]

            if valueIndex == 0:

                xValue = x [0]

            if valueIndex == len (x):

                xValue = x [-1]

            else:

                xValue = x [valueIndex-1] + (xInputCumulativeNormalDistributionValue - cumulativeNormalDistribution [valueIndex-1]) * \
                         (x [valueIndex] - x [valueIndex-1]) / \
                         (cumulativeNormalDistribution [valueIndex] - cumulativeNormalDistribution [valueIndex-1])



            xValues.append (xValue)


        zValues = (np.array (xValues) - xInputMean) / xInputStd
        phiValues = np.array (zValues) * xInputStd + xInputMean
        phiValuesStandard = (phiValues - xInputMean) / xInputStd

        zValuesxInputOrdered = (np.array (xInputOrdered) - xInputMean) / xInputStd

        # Plot the QQ-plot in figure 1 and the histogram in figure 2
        plt.figure (1)
        plt.clf ()

        labelMeanToPrint = 'mean, median = {:g}, {:g}'.format (xInputMean, xInputMedian)
        labelSTDToPrint = 'sd = {:g}'.format (xInputStd)

        plt.scatter (zValuesxInputOrdered, zValues, s=1, color = 'green')
        plt.xlim (-3,3)
        # Plot the vertical 1-sigma reference lines
        plt.plot ([-1,-1], [-3,3], color = 'green', linewidth = 0.5)
        plt.plot ([1,1], [-3,3], color = 'green', linewidth = 0.5)

        plt.plot (phiValuesStandard, zValues, color = 'orange')
        plt.ylim (-3,3)
        # Plot the horizontal 1-sigma reference lines
        plt.plot ([-3,3], [-1,-1], color = 'orange', linewidth = 0.5)
        plt.plot ([-3,3], [1,1], color = 'orange', linewidth = 0.5)

        plt.xlabel (xlabelToPrint + ' - normalised to N(mu,sigma)', fontsize = 12)
        plt.ylabel (ylabelToPrint, fontsize = 12)
        plt.title (QQTitleToPrint)

        if plotTextAverageMedian:

            plt.text (-2.8,2.6, labelMeanToPrint, color = 'green')
            plt.text (-2.8,2.3, labelSTDToPrint, color = 'green')


        plt.figure (2)
        plt.clf ()

        xInputMin = np.min (xInput)
        xInputMax = np.max (xInput)
        xInputRange = xInputMax - xInputMin
        xInputBinSize = xInputRange / 100

        bins = [ (xInputMin + i*xInputBinSize) for i in range (-1, 102) ]
        xInputHistogram = np.histogram (xInput, bins = bins)
        xInputHistogramMax = np.max (xInputHistogram [0])

        plt.hist (xInput, bins = bins, color = 'green')
        plt.xlim (xInputMin - 10*xInputBinSize, xInputMax + 10*xInputBinSize)
        plt.ylim (0, xInputHistogramMax*1.17)
        plt.xlabel (xlabelToPrint, fontsize = 12)
        plt.ylabel ('frequency', fontsize = 12)

        plt.title (HistTitleToPrint)

        if plotTextAverageMedian:
        
            plt.text (xInputMin - 5*xInputBinSize, 1.10*xInputHistogramMax , labelMeanToPrint, color = 'green')
            plt.text (xInputMin - 5*xInputBinSize, 1.04*xInputHistogramMax, labelSTDToPrint, color = 'green')

        plt.show ()



    def linearLeastSquare (xInput, yInput):
        '''
        Determine the values of the variables a and b for the linear least square solution y  =  a * x  +  b
        as well as their uncertainties.
        '''
    
        xInput = np.asarray (xInput)
        yInput = np.asarray (yInput)
        
        # do not take into account any NaN values in yInput
        iValid = np.where ( np.isfinite (yInput) )
        
        x = xInput [iValid]
        y = yInput [iValid]

        numberOfValues = len (x)
        
        sumX = np.sum (x)
        sumX2 = np.sum (x * x)
        
        sumY = np.sum (y)
        sumXY = np.sum (x * y)
        
        determinant = numberOfValues * sumX2 - sumX * sumX
        
        a = ( numberOfValues * sumXY - sumX * sumY ) / determinant
        b = ( sumX2 * sumY - sumX * sumXY ) / determinant

        
        residu = y - ( a * x + b )
        sumResidu2 = np.sum (residu * residu)
        
        averageErrorY = sumResidu2 / (numberOfValues - 2)
        averageErrorA = averageErrorY * numberOfValues / determinant
        averageErrorB = averageErrorY * sumX2 / determinant
        
        uncertaintyA = np.sqrt (averageErrorA)
        uncertaintyB = np.sqrt (averageErrorB)
        
        return a, b, uncertaintyA, uncertaintyB 



#------------------------------------------------
# Three functions needed for the QQplot function


def getNormalDistributionValue (xi,mu,sigma):
    """
    The value of the gaussian normal distribution defined by N(mu, sigma) at xi
    """

    normalValue = np.exp ( -0.5 * (xi - mu) * (xi - mu) / sigma / sigma )

    return normalValue


def getNormalDistribution (mu,sigma):
    """
    Calculate the gaussian normal distribution curve N(mu,sigma) between -5*sigma and +5*sigma at stepsize 0.05*sigma
    """

    cnst = 1 / (sigma * np.sqrt (2 * np.pi))

    x = [ (mu + i*sigma / 20) for i in range (-100, 101) ]

    normalDistribution = [cnst * getNormalDistributionValue (x [i], mu, sigma)  for i in range (len (x))]

    return x, normalDistribution


def getCumulativeNormalDistribution (mu,sigma):
    """
    Calculate the gaussian normal cumulative distribution curve N(mu,sigma) between -5*sigma and +5*sigma at stepsize 0.05*sigma
    """

    x, normalDistribution = getNormalDistribution (mu, sigma)

    cumulativeNormalDistribution = []

    cumulativeNormalDistribution.append (normalDistribution [0])

    for i in range (1,len (x)-1):

        cumulativeNormalDistribution.append ( cumulativeNormalDistribution [-1] + \
                                              (x [i]-x [i-1]) * ( normalDistribution [i] + normalDistribution [i-1]) / 2 \
                                            )

    cumulativeNormalDistribution.append (cumulativeNormalDistribution [-1])

    return x, cumulativeNormalDistribution


#------------------------------------------------
