# HandyTools: a Python pseudo-class
# Author = Maarten Roos

currentVersionHandyTools = '20230529'

import os
import sys
import shutil
separatorCharacter = '\\' if sys.platform == 'win32' else '/'

try:

    from archiveToolConfiguration import *
    EMCSystem = True

except:

    EMCSystem = False


import path
from pathlib import Path

import datetime
import time

import matplotlib.pyplot as plt
import numpy as np

from scipy import signal


class HandyTools:
    """
    HandyTools is a pseudo-class (no instantiation, no 'self'), bundling a couple of functions that can be "handy" at times.
    """


    # A handy function to get the list of absolute paths of all files of a certain extension (default .png) down a directory tree
    def getFilesInDirectoryTree (startPath, extension = '', stringsToExclude = [], checkStartPathOnly = False):
        '''
        **Description:**
        Determine a list of strings that contain all the files paths/names ending on  extension  down the directory tree starting at  startPath .
         startPath  is the point in the directory tree where to start. The function then looks down the tree and seeks out all the files.
        '''

        dirPath = Path (startPath)

        listOfFileNames = []
        if checkStartPathOnly:

            listOfAllFileNames = sorted ( os.listdir (startPath) )         
            if extension == '':

                # The files '.DS_Store' is a MAC OS administration file, that is not of interest to keep.
                listOfFileNames = [ os.path.abspath (startPath + separatorCharacter + fileName)  for fileName in listOfAllFileNames  if fileName != '.DS_Store']
       
        
            else:
                
                listOfFileNames = [ os.path.abspath (startPath + separatorCharacter + fileName)  for fileName in listOfAllFileNames 
                                    if fileName.endswith (extension)
                ]

    
        else:

            if extension == '':

                # The files '.DS_Store' is a MAC OS administration file, that is not of interest to keep.
                listOfFileNames = [ os.path.abspath ( os.path.join (root, name) )
                                    for root, dirs, files in os.walk (dirPath)
                                    for name in sorted (files)
                                    if name != '.DS_Store'
                ]
    
    
            else:

                listOfFileNames = [ os.path.abspath ( os.path.join (root, name) )
                                    for root, dirs, files in os.walk (dirPath)
                                    for name in sorted (files)
                                    if name.endswith (extension)
                ]


    
        # The user can specify certains string that, when part of a file name, this file deleted from the  listOfFileNames
        if len (stringsToExclude):
        
            iDeleteFromList = []
            for iFileName, fileName in enumerate (listOfFileNames):
            
                for stringToExclude in stringsToExclude:
                
                    if stringToExclude in fileName:
                    
                        iDeleteFromList.append (iFileName)  
            

            # If there are files to be deleted from the  listOfFileNames , then delete them.
            if len (iDeleteFromList):
            
                # Sort the indices to be removed in reverse order, larger at the start of the  iDeleteFromList  list. 
                iDeleteFromList = sorted (iDeleteFromList, reverse = True)                
                for iDelete in iDeleteFromList:
                
                    listOfFileNames.pop (iDelete)

        
        # If on the EMC system, then make sure that the absolute paths are all starting with the standardRootPath = '\\\\department-m.erasmusmc.nl\\card\\Data'
        if EMCSystem:
        
            for iFileName in range ( len (listOfFileNames) ):
            
                if 'V:' in listOfFileNames [iFileName]:
                
                    listOfFileNames [iFileName] = listOfFileNames [iFileName].replace ('V:', '\\\\department-m.erasmusmc.nl\\card\\Data')


                if '\\\\storage.erasmusmc.nl\\v\\vcl09\\CARD\\Data' in listOfFileNames [iFileName]:
                
                    listOfFileNames [iFileName] = listOfFileNames [iFileName].replace ('\\\\storage.erasmusmc.nl\\v\\vcl09\\CARD\\Data', '\\\\department-m.erasmusmc.nl\\card\\Data')
                    


        return listOfFileNames



    # Get the absolute path for a file.
    def getFileAndAbsolutePath (fileName):
        '''
        :param fileName: file name (and path) of the file.
        
        :return: absolute path and file name, absolute directory, just the file name.
        :rtype: str, str, str

        **Description:**    
        Get the absolute path, directory and root (just the file name) for a file.
        If '' are returned, then the file does not exist and an error message is printed to the Python console.
        '''
        
        absolutePath = ''
        fileRootName = ''
        absoluteDirectory = ''

        if os.path.isfile (fileName):
        
            absolutePath = os.path.abspath (fileName)
            fileRootName = absolutePath.split ('/') [-1]
            absoluteDirectory = os.path.dirname (absolutePath)
            
        else:
        
            print ()
            print ( ' WARNING - File  {}  does not exist'.format (fileName) )
            print ()


        return absolutePath, absoluteDirectory, fileRootName


  
    # Create the full path for a given file with full path or just full path.   
    def createPathToFile (fileNameAndFullPath = '', fullPath = ''):
        '''
        :param fileNameAndFullPath: file name with full path.
        :type fileNameAndFullPath: stringsToExclude
        
        :param fullPath: full Path
        :type fullPath: str
        
        :return: full path
        :rtype: str
        
        **Description:**
        Create the full path for a given file with full path (``fileNameAndFullPath`` string) or just full path (``fullPath`` string).
        Return the full path that has been created as a string.
        If not successful (full path already exists or cannot be created due to writing permission restrictions), then return an empty string.
        '''


        if fileNameAndFullPath:
        
            fullPathToCreate = fileNameAndFullPath.split ( fileNameAndFullPath.split (separatorCharacter)[-1] ) [0]


        if fullPath:
        
            fullPathToCreate = fullPath

        
        try:
        
            os.makedirs (fullPathToCreate)
            return fullPathToCreate
            
        except:
        
            return ''



    # Read and return the content of a text file.
    def getTextFileContent (textFileNameAndPath, stripLineFromBlanks = False):
        '''
        :param textFileNameAndPath:
        :type textFileNameAndPath: str
        
        :param stripLineFromBlanks:
        :type stripLineFromBlanks: bool
        
        :return: content of file with \n chopped off.
        :rtype: list (str)

        D**Description:**
        Open, read and return the content of a text file.
        Make sure to delete any \n characters at the end of lines that may exist.
        If the stripLineFromBlanks bpolean is True, then also strip any blank spaces at the beginning and end of a string.
        Returns empty list if the file does not exist or there is an error in the reading.
        '''

        if os.path.isfile (textFileNameAndPath):

            try:
                              
                fileOpen = open (textFileNameAndPath, 'r')
                fileContent = fileOpen.readlines ()
                fileOpen.close ()

                fileContentClean = []
                for fileLine in fileContent:
                                
                    fileContentClean.append ( fileLine [:-1] if fileLine [-1] == '\n'  else  fileLine )

                    if stripLineFromBlanks:
                    
                        fileContentClean [-1].strip ()

        
                return fileContentClean
                
            except:
                
                print ('')
                print ('---WARNING---')
                print (' From HandyTools.getTextFileContent: ')
                print ('  file {} cannot be opened and / or read correctly.'.format (textFileNameAndPath))                
                return []

            
        else:
 
            print ('')
            print ('---WARNING---')
            print (' From HandyTools.getTextFileContent: ')
            print ('  Warning: file {} does not exist!'.format (textFileNameAndPath))                
            return []
        
            
                   
    # Save content (list, dictionary, ...) to a numpy file with a custom extension.
    def saveContentToNumpyWithCustomExtension (contentToSave, fileName, extensionWithoutDot, overWrite = False):
        '''
        :param contentToSave:
        
        :param fileName:
        :type fileName: str
 
        :param extensionWithoutDot:
        :type extensionWithoutDot: str
        
        :param overWrite: 
        :type overWrite: bool
        
        :return: file saved successful, file already exists.
        :rtype: bool, bool
        
        **Description:**
        Save content (list, dictionary, ...) to a numpy file with a custom extension.
        '''
    
        fileNameWithExtension = fileName + '.' + extensionWithoutDot
        fileSaved = False

        # Only attempt to save if the file does not yet exist.
        if not os.path.isfile (fileNameWithExtension) or overWrite:
        
            fileAlreadyExists = False            

            try:
        
                np.save (fileNameWithExtension, contentToSave)
                
                # os.rename does not work on Windows when the destination file already exists.
                shutil.move (fileNameWithExtension + '.npy', fileNameWithExtension)
        
                fileSaved = True

            except:
        
                print ( 'Warning: file {} already exists, it was not overwritten!'.format (fileNameWithExtension) )
 
 
        else:
            
            fileAlreadyExists = True

            
        return fileSaved, fileAlreadyExists
    


    # Get the date and time now.
    def getDateAndTime ():
        '''
        **Description:**
        Get the date and time now.
        '''
        
        # This try - except loop is needed because apparently between python versions the date module has changed some of its structure
        try:

            dateAndTime = datetime.now ()

        except:

            dateAndTime = datetime.datetime.now ()

                
        return dateAndTime
 


    # Get the date and time now and return in a string format.
    def getDateAndTimeString (includeYMD = True, dateFormat = 'YMD', includeHMS = True):
        '''
        **Description:**
        Get the date and time now and return in a string format.
        '''
      
        dateAndTime = HandyTools.getDateAndTime ()  
                
        YMD = '{}-{}-{}'.format ( str (dateAndTime.year), str (dateAndTime.month).zfill(2), str (dateAndTime.day).zfill(2) )
        DMY = '{}-{}-{}'.format ( str (dateAndTime.day).zfill(2), str (dateAndTime.month).zfill(2), str (dateAndTime.year) )
        MDY = '{}-{}-{}'.format ( str (dateAndTime.month).zfill(2), str (dateAndTime.day).zfill(2), str (dateAndTime.year) )

        HMS = '{}:{}:{}'.format ( str (dateAndTime.hour).zfill(2), str (dateAndTime.minute).zfill(2), str (dateAndTime.second).zfill(2) )

        dateAndTimeString = ''
        if includeYMD:
        
            if dateFormat == 'YMD': dateAndTimeString += YMD
            if dateFormat == 'DMY': dateAndTimeString += DMY
            if dateFormat == 'MDY': dateAndTimeString += MDY

        if includeHMS and includeYMD: dateAndTimeString += ' at '
        
        if includeHMS: dateAndTimeString += HMS
        
        return dateAndTimeString
   


    # Calculate  Hours Minutes Seconds  from a total number of seconds
    def getHMSFromTotalNumberOfSeconds (numberOfSecondsTotal = 0, numberOfDigitsAccuracy = 4):
        '''
        **Description:**
        Calculate  Hours Minutes Seconds
        '''
    
        numberOfHours = int (numberOfSecondsTotal / 3600)
        numberOfMinutes = int ( (numberOfSecondsTotal - numberOfHours * 3600) / 60 )
        
        accuracyFactor = 10**numberOfDigitsAccuracy
        numberOfSeconds = int ( ( numberOfSecondsTotal - numberOfHours * 3600 - numberOfMinutes * 60 ) * accuracyFactor ) / accuracyFactor

        return [numberOfHours, numberOfMinutes, numberOfSeconds], \
               [ '{:02d}'.format (numberOfHours), '{:02d}'.format (numberOfMinutes), '{:7.4f}'.format (numberOfSeconds) ]                 


    
    # Calculate the total number of seconds from an Hours Minutes Seconds list.
    def getTotalNumberOfSecondsfromHMS (HMS = [0,0,0]):
        '''
        **Description:**
        Calculate the total number of seconds from an Hours Minutes Seconds list.
        HMS is a list of [int, int, float] or a list of three strings.
        '''
        
        try:
        
           numberOfHours = HMS [0]
           numberOfMinutes = HMS [1]
           numberOfSeconds = HMS [2]
           
        except:
        
           numberOfHours = int ( HMS [0] )
           numberOfMinutes = int ( HMS [1] )
           numberOfSeconds = float ( HMS [2] )
            
            
        return numberOfHours * 3600 + numberOfMinutes * 60 + numberOfSeconds



    # Determine and print execution of a piece of code.
    def getRunTime (functionName, startTime = None, indent = 0, printResult = True):
        '''
        Determine and print execution of a piece of code.
        '''
        
        # Get the start time of the function's execution. This is the default.
        if not startTime:
        
            print ( ' ' * indent * 3 + ' --> Start time of {} retrieved'.format (functionName) )
            return time.time ()
        

        # If  printResult  is True, then print the end time, given the  startTime.        
        elif printResult:
        
            print ( ' ' * indent * 3 + ' --> Run time for {} = {:9.4f}s'.format (functionName, time.time () - startTime) )
        
        # If  printResult  is False, then return the end time, given the  startTime.
        else:
        
            return time.time () - startTime
        


    # Used in the  getUncertaintyLevelInElectrode  method and by  AnnotationTool  in te NoiseViewer method.
    # Determine the list of amplitude segments for the electrogram.
    def getListOfAmplitudeSegmentsFromDataValuesList (dataValuesList):
        '''
        Determine the list of amplitude segments from a list of wobbling data values.
        '''

        np.seterr (all = 'ignore')

        segmentAmplitudes = []
        segmentAmplitudes.append (dataValuesList [1] - dataValuesList [0])

        segmentStartindices = []
        segmentStartindices.append (0)

        segmentStartIndex = 0
        for iSample in range (1, len (dataValuesList) - 1):

            deltaValue = dataValuesList [iSample + 1] - dataValuesList [iSample]

            if (segmentAmplitudes [-1] > 0 and deltaValue >= 0) or (segmentAmplitudes [-1] < 0 and deltaValue <= 0):

                segmentAmplitudes [-1] += deltaValue

            else:

                segmentAmplitudes.append (deltaValue)

                segmentStartindices.append (segmentStartIndex)
                segmentStartIndex = iSample


        segmentAmplitudes = np.asarray (segmentAmplitudes)        
        segmentStartindices = np.asarray (segmentStartindices)
                
        return segmentAmplitudes, segmentStartindices
   


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
                QQTitleToPrint = 'HandyTools version ' + currentVersionHandyTools + ': QQ-plot of input data',
                HistTitleToPrint = 'HandyTools version ' + currentVersionHandyTools + ': Histogram of input data', 
                plotTextAverageMedian = True):
        '''
        **Description:**
        Calculate and plot the QQ-plot (or Quantile-Quantile plot) and the histogram of a given list of input values.
        '''

        xInput = np.array (xInput)
        xInputMean = np.mean (xInput)
        xInputMedian = np.median (xInput)
        xInputStd = np.std (xInput)

        x, cumulativeNormalDistribution = HandyTools.getCumulativeNormalDistribution (xInputMean, xInputStd)

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

        x, normalDistribution = HandyTools.getNormalDistribution (mu, sigma)

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

        normalDistribution = [cnst * HandyTools.getNormalDistributionValue (x [i], mu, sigma)  for i in range (len (x))]

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

