# DataTools: a Python pseudo-class
# Author = Maarten Roos

DataToolsVersion = '20240117'

# Standard imports.
import os
import sys
import shutil
separatorCharacter = '\\' if sys.platform == 'win32' else '/'


from pathlib import Path

import datetime
import time

import matplotlib.pyplot as plt
import numpy as np

from scipy import signal


# Custom imports

import DataWranglingToolsPYtoCPP
import FilterToolsPYtoCPP




class DataTools:
    """
    DataTools is a pseudo-class (no instantiation, no 'self'), bundling a couple of functions to work with all sorts of data.
    """


    # Used in the  getUncertaintyLevelInElectrode  method and by  AnnotationTool  in te NoiseViewer method.
    # Determine the list of amplitude segments for the electrogram.
    def getListOfAmplitudeSegmentsFromDataValues (dataValues, PYtoCPP = True):
        '''    
        dataValues: 
    
            converted to a short (np.short / int16)

    
        returns tuple:
    
            [0]  numberOfSegments
            [1]  segmentStartIndices [0:numberOfSegments]
            [2]  segmentAmplitudes [0:numberOfSegments]
            [3]  segmentSlopes [0:numberOfSegments]
            [4]  segmentDurations [0:numberOfSegments]
            [5]  numberOfSegmentsNegative
            [6]  segmentStartIndicesNegative [0:numberOfSegmentsNegative]
            [7]  iSteepestNegativeSlopeSegment
            [8]  iSegmentStartIndicesSteepestNegativeSlope
            [9]  numberOfSegmentsPositive
            [10] segmentStartIndicesPositive [0:numberOfSegmentsPositive]
            [11] iSteepestPositiveSlopeSegment
            [12] iSegmentStartIndicesSteepestPositiveSlope
            
        **Description:**        
        Determine the list of amplitude segments from a list of data values.
        Note that values will be converted to int16 (short)!
        
        Also note that the deltaValue HAS TO BE an int32 at least, and not int16 hence the int-function in the Python version (same)
        '''

        # Run this function with the C++ core per default.
        if PYtoCPP:

            return DataWranglingToolsPYtoCPP.getListOfAmplitudeSegmentsPYtoCPP (dataValues)


        # Run the Python version (which needs updating to match the C++ code!).
        else:
                
            np.seterr (all = 'ignore')

            segmentAmplitudes = []
            segmentAmplitudes.append (int (dataValues [1]) - int (dataValues [0]))

            segmentStartindices = []
            segmentStartindices.append (0)

            for iSample in range (1, len (dataValues) - 1):
            
                deltaValue = int (dataValues [iSample + 1]) - int (dataValues [iSample])

                if (segmentAmplitudes [-1] > 0 and deltaValue >= 0) or (segmentAmplitudes [-1] < 0 and deltaValue <= 0):

                    segmentAmplitudes [-1] += deltaValue

                else:

                    segmentAmplitudes.append (deltaValue)
                    segmentStartindices.append (iSample)


            segmentAmplitudes = np.asarray (segmentAmplitudes)        
            segmentStartindices = np.asarray (segmentStartindices)
                
            return segmentAmplitudes, segmentStartindices


    #
    def passAverageFilter (listOfNumbers, widthOfWindow):
        '''
        '''
        

        return FilterToolsPYtoCPP.passAverageFilterPYtoCPP (listOfNumbers, widthOfWindow)
        
       

    # Pass a given input signal through a notch filter.
    def getNotchFilteredSignal (inputSignal, samplingFrequency = 1000, notchFrequency = 50, qualityFactor = 20):
        '''
        Pass a given input signal through a notch filter.
        '''
         
        # Design a notch filter using the signal.iirnotch method (Infinite Impulse Response)
        bNotch, aNotch = signal.iirnotch (notchFrequency, qualityFactor, samplingFrequency)
 
        # Compute magnitude response of the designed filter
        filterFrequency, amplitudedB = signal.freqz (bNotch, aNotch, fs = 2 * np.pi)
        
        outputSignal = signal.filtfilt (bNotch, aNotch, inputSignal)

        return outputSignal, filterFrequency, amplitudedB



    # Determine the values of the variables a and b for the linear least square solution y  =  a * x  +  b.
    def linearLeastSquare (xInput, yInput):
        '''
        **Description:**
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


    # Calculate and plot the QQ-plot (or Quantile-Quantile plot) and the histogram of a given list of input values.
    def QQPlot (xInput, xlabelToPrint = 'input values', ylabelToPrint = 'z (sigma)',
                QQTitleToPrint = 'DataTools version ' + DataToolsVersion + ': QQ-plot of input data',
                HistTitleToPrint = 'DataTools version ' + DataToolsVersion + ': Histogram of input data', 
                plotTextAverageMedian = True):
        '''
        **Description:**
        Calculate and plot the QQ-plot (or Quantile-Quantile plot) and the histogram of a given list of input values.
        '''

        xInput = np.array (xInput)
        xInputMean = np.mean (xInput)
        xInputMedian = np.median (xInput)
        xInputStd = np.std (xInput)

        x, cumulativeNormalDistribution = DataTools.getCumulativeNormalDistribution (xInputMean, xInputStd)

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



    # Calculate the gaussian normal cumulative distribution curve N(mu,sigma) between -5*sigma and +5*sigma at stepsize 0.05*sigma.
    def getCumulativeNormalDistribution (mu,sigma):
        '''
        **Description:**
        Calculate the gaussian normal cumulative distribution curve N(mu,sigma) between -5*sigma and +5*sigma at stepsize 0.05*sigma.

        (Used by the QQPlot function)
        '''

        x, normalDistribution = DataTools.getNormalDistribution (mu, sigma)

        cumulativeNormalDistribution = []

        cumulativeNormalDistribution.append (normalDistribution [0])

        for i in range (1,len (x)-1):

            cumulativeNormalDistribution.append ( cumulativeNormalDistribution [-1] + \
                                                  (x [i]-x [i-1]) * ( normalDistribution [i] + normalDistribution [i-1]) / 2 \
                                                )

        cumulativeNormalDistribution.append (cumulativeNormalDistribution [-1])

        return x, cumulativeNormalDistribution



    # Calculate the gaussian normal distribution curve N(mu,sigma) between -5*sigma and +5*sigma at stepsize 0.05*sigma.
    def getNormalDistribution (mu,sigma):
        '''
        **Description:**
        Calculate the gaussian normal distribution curve N(mu,sigma) between -5*sigma and +5*sigma at stepsize 0.05*sigma.
        
        (Used by the getCumulativeNormalDistribution function).
        '''

        cnst = 1 / (sigma * np.sqrt (2 * np.pi))

        x = [ (mu + i*sigma / 20) for i in range (-100, 101) ]

        normalDistribution = [cnst * DataTools.getNormalDistributionValue (x [i], mu, sigma)  for i in range (len (x))]

        return x, normalDistribution



    # The value of the gaussian normal distribution defined by N(mu, sigma) at xi.
    def getNormalDistributionValue (xi,mu,sigma):
        '''
        **Description:**
        The value of the gaussian normal distribution defined by N(mu, sigma) at xi.
        
        (Used by the getNormalDistribution function).
        '''

        normalValue = np.exp ( -0.5 * (xi - mu) * (xi - mu) / sigma / sigma )

        return normalValue

