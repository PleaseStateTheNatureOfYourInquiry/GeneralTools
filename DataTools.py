# DataTools: a Python pseudo-class
# Author = Maarten Roos

DataToolsVersion = '20240212'

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


# Custom imports of the Python to C++ libraries.
import DataWranglingToolsPYtoCPP
import FilterToolsPYtoCPP



class DataTools:
    """
    DataTools is a pseudo-class (no instantiation, no 'self'), bundling a couple of functions to work with all sorts of data.
    """


    # Used in the  getUncertaintyLevelInElectrode  method and by  AnnotationTool  in te NoiseViewer method.
    # Determine the list of amplitude segments for the electrogram.
    @staticmethod
    def getSegmentSpecsFromDataValues (dataValues, PYtoCPP = True):
        '''    
        :param dataValues: list (one dimension) of data values. The values in the ``dataValues`` list are converted to int16 type (for now).
        :type dataValues: list

        :param PYtoCPP: at the time of writing, the user can still chose to use Python to analyse the ``dataValues`` by setting ``PYtoCPP`` to ``False``. The result is less comprehensive and the runtime significantly longer.
        :type PYtoCPP: boolean; default = True

        :return: all the characteristics of the segments comprised by the list of ``dataValues`` (see **Description** below) when ``PYtoCPP = True``, or list of amplitudes and start indices when ``PYtoCPP = False``.
        :rtype: tuple when ``PYtoCPP = True`` or two lists, first list [same type as ``dataValues``] and second list [int]


        **Description:**
        With this function a list of (wiggling) data values is analysed to determine all the specifications in terms of its segments. A segment in a wiggling data values line is a section of the line from a local minimum to the next local maximum or a local maximum to the next local minimum. The specifications are returned in terms of a tuple of NumPy arrays and numbers:
    
            | [0]  numberOfSegments
            | [1]  segmentStartIndices - the start indices of each segment (the local minima and maxima).
            | [2]  segmentAmplitudes - the corresponding amplitudes between local extremes in units of ADU.
            | [3]  segmentSlopes - the corresponding slopes in units of ADU / time sample.
            | [4]  segmentDurations - the number of time samples between local extremes.
            | [5]  numberOfSegmentsNegative - subsection of the list above with only the segments from local maxima to the next local minima.
            | [6]  segmentStartIndicesNegative [0:numberOfSegmentsNegative]
            | [7]  iSteepestNegativeSlopeSegment
            | [8]  iSegmentStartIndicesSteepestNegativeSlope
            | [9]  numberOfSegmentsPositive - subsection of the list above with only the segments from local minima to the next local maxima.
            | [10] segmentStartIndicesPositive [0:numberOfSegmentsPositive]
            | [11] iSteepestPositiveSlopeSegment
            | [12] iSegmentStartIndicesSteepestPositiveSlope
 
        When calling this function with :code:`PYtoCPP = False`, only list of segment amplitudes and segment start indices is being returned.           
        '''

        # Run this function with the C++ core per default.
        if PYtoCPP:

            return DataWranglingToolsPYtoCPP.getSegmentSpecsFromDataValuesPYtoCPP (dataValues)


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
    @staticmethod
    def passAverageFilter (dataValues, widthOfWindow):
        '''
        :param dataValues: list (one dimension) of data values that represent the signal to be filtered.
        :type dataValues: list 

        :param widthOfWindow: half the width of the window: number of values to the left and to the right of the central value. Total number of point in the window is :code:`2 * widthOfWindow + 1`.
        :type widthOfWindow: int

        :return: signal ``dataValues`` filtered with a running average. 
        :rtype: list 


        **Description:**
        Use this function to perform a running average filtering of a list (one dimension) of data values. The total number of point in the averaging window is :code:`2 * widthOfWindow + 1`. At the beginning and the end of the list, filtering is done with the available data values. For example, the result for the first data value is the average of the first data value and the :code:`width` data values to the right. 
        
        The noise level of the filtered set can be easily calculated from the noise level of the original data values, by dividing by the square root of the total number of data values in the averaging window. 
        '''
        

        return FilterToolsPYtoCPP.passAverageFilterPYtoCPP (dataValues, widthOfWindow)
        
       

    # Pass a given input signal through a notch filter.
    @staticmethod
    def getNotchFilteredSignal (inputSignal, samplingFrequency = 1000, notchFrequency = 50, qualityFactor = 20):
        '''
        :param inputSignal: list (one dimension) of data values that represent the signal to be filtered.
        :type inputSignal: list 

        :param samplingFrequency: 
        :type samplingFrequency: float

        :param notchFrequency: 
        :type notchFrequency: float

        :param qualityFactor: 
        :type qualityFactor: float

        :return: outputSignal, filterFrequency, amplitudedB
        :rtype: list, float, float
        
        
        **Description:**
        Pass a given input signal through a notch filter, using the signal.iirnotch, signal.freqz and the signal.filtfilt (back and forth filtering) methods.
        '''
         
        # Design a notch filter using the signal.iirnotch method (Infinite Impulse Response)
        bNotch, aNotch = signal.iirnotch (notchFrequency, qualityFactor, samplingFrequency)
 
        # Compute magnitude response of the designed filter
        filterFrequency, amplitudedB = signal.freqz (bNotch, aNotch, fs = 2 * np.pi)
        
        outputSignal = signal.filtfilt (bNotch, aNotch, inputSignal)

        return outputSignal, filterFrequency, amplitudedB



    # Determine the values of the variables a and b for the linear least square solution y  =  a * x  +  b.
    @staticmethod
    def linearLeastSquare (xInput, yInput):
        '''
        :param xInput: the x-values of the data to fit.
        :type xInput: list [float]

        :param yInput: the y-values of the data to fit.
        :type yInput: list [float]

        :return: a, b, uncertaintyA, uncertaintyB
        :rtype: float, float, float, float
        
        
        **Description:**
        Determine the values of the variables a and b for the linear least square solution y (x) =  ax + b, along with their uncertainties, for a set of data points described by the x- and y-values.
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
    @staticmethod
    def QQPlot (xInput, xlabelToPrint = 'input values', ylabelToPrint = 'z (sigma)',
                QQTitleToPrint = 'DataTools version ' + DataToolsVersion + ': QQ-plot of input data',
                HistTitleToPrint = 'DataTools version ' + DataToolsVersion + ': Histogram of input data', 
                plotTextAverageMedian = True):
        '''
        :param xInput: the list of data points.
        :type xInput: list [float]

        :param xlabelToPrint: the x-label to print on the QQ-plot.
        :type xlabelToPrint: str

        :param ylabelToPrint: the x-label to print on the QQ-plot.
        :type ylabelToPrint: str

        :param QQTitleToPrint: the title to print on the QQ-plot.
        :type QQTitleToPrint: str

        :param HistTitleToPrint: the title to print on the histogram.
        :type HistTitleToPrint: str

        :param plotTextAverageMedian: if True then plot the values of the average, median and standar deviation in the plot.
        :type plotTextAverageMedian: bool

                
        **Description:**
        Calculate and plot the QQ-plot (or Quantile-Quantile plot), as well as the histogram, of a given list of data values.
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
    @staticmethod
    def getCumulativeNormalDistribution (mu,sigma):
        '''
        :param mu: the average values of the data values.
        :type mu: list [float]

        :param sigma: the standard deviation of the data values.
        :type sigma: float.
        
        :return: x and 
        :rtype: list [float], list [float]
        
        
        **Description:**
        Calculate the gaussian normal cumulative distribution curve N(mu,sigma) between -5 x sigma and +5 x sigma at stepsize 0.05 x sigma.
        The results of the :py:meth:`~.getNormalDistribution` function call is used. 
        This function is called from the :py:meth:`~.QQPlot` function.
        '''

        x, normalDistribution = DataTools.getNormalDistribution (mu, sigma)

        cumulativeNormalDistribution = []

        cumulativeNormalDistribution.append ( normalDistribution [0] )

        for i in range ( 1, len (x) - 1 ):

            cumulativeNormalDistribution.append ( cumulativeNormalDistribution [-1] + 
                                                  ( x [i]- x [i-1] ) * ( normalDistribution [i] + normalDistribution [i-1] ) / 2 )

        cumulativeNormalDistribution.append ( cumulativeNormalDistribution [-1] )

        return x, cumulativeNormalDistribution



    # Calculate the gaussian normal distribution curve N(mu,sigma) between -5*sigma and +5*sigma at stepsize 0.05*sigma.
    @staticmethod
    def getNormalDistribution (mu,sigma):
        '''
        :param mu: the average values of the data values.
        :type mu: float

        :param sigma: the standard deviation of the data values.
        :type sigma: float
 
        :return: the equally separated x values along the gaussian distribution and the corresponding values of a pure normal distribution.
        :rtype: list [float], list [float]
               
        **Description:**
        Calculate the gaussian normal distribution curve N (mu,sigma) between -5 x sigma and +5 x sigma at stepsize 0.05 x sigma.
        This function calls the :py:meth:`~.getNormalDistributionValue` function and it is called by the :py:meth:`~.getCumulativeNormalDistribution` function.
        '''

        constant = 1 / (sigma * np.sqrt (2 * np.pi))

        x = [ (mu + i*sigma / 20) for i in range (-100, 101) ]

        normalDistribution = [ constant * DataTools.getNormalDistributionValue (x [i], mu, sigma)  for i in range ( len (x) ) ]

        return x, normalDistribution



    # The value of the gaussian normal distribution defined by N(mu, sigma) at xi.
    @staticmethod
    def getNormalDistributionValue (xi,mu,sigma):
        '''
        :param xi: the x-value to calculate the normal distribution value for.
        :type xi: float

        :param mu: the average values of the data values.
        :type mu: float

        :param sigma: the standard deviation of the data values.
        :type sigma: float
 
        :return: the normal distribution value for the point x-i in the normal distribution with average mu and standard deviation sigma.
        :rtype: float
        
        
        **Description:**
        The value of the gaussian normal distribution defined by N (mu, sigma) at xi.
        This function is called by the :py:meth:`~.getNormalDistribution` function.
        '''

        normalValue = np.exp ( -0.5 * (xi - mu) * (xi - mu) / sigma / sigma )

        return normalValue



    #
    @staticmethod
    def getAverageVarAndSDPYtoCPP (dataValues = [], removeNaN = False):
        '''
        :param dataValues: list of data values.
        :type dataValues: list [float] or NumPy array (one dimension)

        :param removeNaN: if True then remove any NaN values from the list of data values.
        :type removeNaN: bool
        
        :return: average, standard deviation and variance of the list of data values.
        :rtype: float, float, float.
        
        **Description:**
        Calculate the average, standard deviation and variance of the list of data values using C++ code.   
        If :code:`removeNaN = True`, then call the :py:meth:`~.getNanFreeNumpyArray` function to remove any NaN values from the data list.     
        '''

        if removeNaN:
              
            dataValues = DataTools.getNanFreeNumpyArray (dataValues)

    
        if len (dataValues):
                    
            averagevalues, standardDeviation, variance = DataWranglingToolsPYtoCPP.getAverageVarAndSDPYtoCPP (dataValues)
            return averagevalues, standardDeviation, variance

        else:
        
            return None, None, None


    #
    @staticmethod
    def getMedianAndQuantilesPYtoCPP (dataValues = [], lowerQuantilePercentage = 25, upperQuantilePercentage = 75, removeNaN = False):
        '''
        :param dataValues: list of data values.
        :type dataValues: list [float] or NumPy array (one dimension)

        :param lowerQuantilePercentage: the lower quantile.
        :type lowerQuantilePercentage: float

        :param upperQuantilePercentage: the upper quantile.
        :type upperQuantilePercentage: float

        :param removeNaN: if True then remove any NaN values from the list of data values.
        :type removeNaN: bool
        
        :return: average, standard deviation and variance of the list of data values.
        :rtype: float, float, float.
        
        **Description:**
        Calculate the median, lower and upper quantiles of the list of data values using C++ code.   
        If :code:`removeNaN = True`, then call the :py:meth:`~.getNanFreeNumpyArray` function to remove any NaN values from the data list.     
        '''

        if removeNaN:
              
            dataValues = DataTools.getNanFreeNumpyArray (dataValues)
 
    
        if len (dataValues):
    
            medianValue, lowerQuantileValue, upperQuantileValue = \
             DataWranglingToolsPYtoCPP.getMedianAndQuantilesPYtoCPP (dataValues, lowerQuantilePercentage / 100, upperQuantilePercentage / 100)   
            return medianValue, lowerQuantileValue, upperQuantileValue

        else:
        
            return None, None, None


    #
    @staticmethod
    def getNanFreeNumpyArray (dataValues):
        '''
        :param dataValues: complete list of data values.
        :type dataValues: list [float] or NumPy array (one dimension)
        
        :return: list of data values with any NaN values removed.
        :rtype: NumPy array 
        
        **Description:**
        Remove all the NaN values from a list of data values and return the cleansed list as a NumPy array (one dimension).
        '''
        
        if type (dataValues) == list:
                
            iNaNValues = np.isnan ( np.asarray (dataValues) )           
            return np.asarray (dataValues) [~iNaNValues]
        
        else:
        
            iNaNValues = np.isnan (dataValues)            
            return dataValues [~iNaNValues]
            


