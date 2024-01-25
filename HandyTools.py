# HandyTools: a Python pseudo-class
# Author = Maarten Roos

currentVersionHandyTools = '20240117'

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
try:

    from archiveToolConfiguration import *
    EMCSystem = True

except:

    EMCSystem = False



from DataTools import DataTools




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
                                    if fileName.endswith (extension) ]

    
        else:

            if extension == '':

                # The files '.DS_Store' is a MAC OS administration file, that is not of interest to keep.
                listOfFileNames = [ os.path.abspath ( os.path.join (root, name) )
                                    for root, dirs, files in os.walk (dirPath)
                                    for name in sorted (files)
                                    if name != '.DS_Store' ]
    
    
            else:

                listOfFileNames = [ os.path.abspath ( os.path.join (root, name) )
                                    for root, dirs, files in os.walk (dirPath)
                                    for name in sorted (files)
                                    if name.endswith (extension) ]


    
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
        If the stripLineFromBlanks boolean is True, then also strip any blank spaces at the beginning and end of a string.
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
        
            if printResult:
            
                print ( ' ' * indent * 3 + ' --> Start time of {} retrieved'.format (functionName) )
        
            return time.time ()
        

        # If  printResult  is True, then print the end time, given the  startTime.        
        elif printResult:
        
            print ( ' ' * indent * 3 + ' --> Run time for {} = {:9.4f}s'.format (functionName, time.time () - startTime) )
        
        # If  printResult  is False, then return the end time, given the  startTime.
        else:
        
            return time.time () - startTime
        
