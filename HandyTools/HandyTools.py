# HandyTools: a Python pseudo-class
# Author = Maarten Roos

currentVersionHandyTools = '20240501'

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



class HandyTools:
    """
    HandyTools is a pseudo-class (no instantiation, no 'self'), bundling a couple of functions that can be "handy" at times.
    """


    # A handy function to get the list of absolute paths of all files of a certain extension (default .png) down a directory tree
    @staticmethod
    def getFilesInDirectoryTree (startPath, extension = '', stringsToExclude = [], checkStartPathOnly = False):
        '''
        :param startPath: directory where to start the search.
        :type startPath: str

        :param extension: extension of files to return. If ``extension = ''``, then all files are returned. If ``extension = ''`` and ``checkStartPathOnly = True``, then all files and directories in ``startPath`` are returned.
        :type extension: str

        :param stringsToExclude: list of strings that when in the file name are excluded from the file list.
        :type stringsToExclude: list [str]

        :param checkStartPathOnly: if ``True``, then only search the startPath folder.
        :type checkStartPathOnly: bool; default = ``False`` 

        :return: list of file names and absolute paths of all the files ending on ``extension`` down the directory tree starting at ``startPath``.
        :rtype: list [str]  


        **Description:**
        This function is used to perform a recursively search for all files with a user defined ``extension``, starting at the user defined ``startPath`` in the directory tree. 
        Files that contain any of the strings in the  ``stringsToExclude`` list, are excluded from the returned file list. If the boolean `checkStartPathOnly`` is set to ``True``, 
        then the search is restricted to the ``startPath`` folder only.
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


        return listOfFileNames



    # Get the absolute path for a file.
    @staticmethod
    def getFileAndAbsolutePath (fileName):
        '''
        :param fileName: file name (and path) of the file.
        
        :return: absolute path and file name, absolute directory, file name without path.
        :rtype: str, str, str

        **Description:**    
        Get the full absolute path, absolute directory and root (= file name without path) of a file.
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
    @staticmethod
    def createPathToFile (fileNameAndFullPath = '', fullPath = ''):
        '''
        :param fileNameAndFullPath: file name with full path.
        :type fileNameAndFullPath: str
        
        :param fullPath: full path.
        :type fullPath: str
        
        :return: full path created by the system.
        :rtype: str
        
        **Description:**
        Create the full path (all the folders in the directory tree) for a given file with full path ``fileNameAndFullPath``, or just full path ``fullPath``.
        Return the full path that has been created as a string.
        If not successful, because full path already exists or cannot be created due to writing permission restrictions, then return an empty string.
        The ``os.makedirs`` method is used.
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
    @staticmethod
    def getTextFileContent (textFileNameAndPath, stripLineFromBlanks = False):
        '''
        :param textFileNameAndPath: file name (and path) of the text file to load.
        :type textFileNameAndPath: str
        
        :param stripLineFromBlanks: set to ``True`` to strip empty spaces from the beginning and end of each line.
        :type stripLineFromBlanks: bool; default = False
        
        :return: content of file with *\\\\n* chopped off.
        :rtype: list (str)

        **Description:**
        Open, read and return the content of a text file.
        Any *\\\\n* (= next line) characters at the end of lines are stripped.
        If the ``stripLineFromBlanks`` boolean is set to ``True``, then also strip any blank spaces at the beginning and end of the lines.
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
                                        
                        fileContentClean [-1] = fileContentClean [-1].strip ()
                        
        
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
        
            

    # Reads data from a tabular text file and return a list of NumPy arrays, strings and the header.
    @staticmethod
    def readTable (textFileNameAndPath, endOfHeaderString = 'C_END', separatorString = ' '):
        '''
        :param textFileNameAndPath: file name (and path) of the table (text) file to read.
        :type textFileNameAndPath: str

        :param endOfHeaderString: string at the end of the header and the start of the tabular data, default = *C_END*
        :type endOfHeaderString: str

        :param separatorString: string that separates the elements in the table at each row, default is one empty space.
        :type separatorString: str

        :return: table content for each column as numbers (int or float, or NaN), as strings, and the content of a header.
        :rtype: list [ list [1-D NumPy array] ], one for each column in the list, list [ list [str] ], list [str]

        **Description:**
        With this method formatted human-readable tabular data from a text file can be loaded as a list of 1-Dimensional NumPy arrays.

        If the table file has header text, then the user can indicate the ``endOfHeaderString`` that signals the end of the header and the start of the data,
        in the example below ``C_END``, which is also the default value.
        
        .. code-block:: console
        
            header line 1
            header line 2
            
            C_END
            entry1 1 4.5 6
            entry2 2 5.3 7
            
        
        Alternatively there can be no header and just data. 
        
        .. code-block:: console

            entry1 1 4.5 6
            entry2 2 5.3 7
        
        
        .. attention:: 
        
            It is assumed there are no empty row amongst the data rows!

                 
        The tabular data is returned both as an array of NumPy arrays and as a list of strings. 
        Each array (list) in the list contains the content of one of the columns from the table.
        Note that the format of the data (int or float) can be different per column and is taken to be the format of the data entrees in the first row.
        In this manner integers and floats can be mixed as long as they are consistent throughout each column. 
        In the example above, the 2nd and  4th columns are int and the 3rd column a float. The first column is not numerical and is returned as a list of NaN.
        
        The ``HandyTools.getTextFileContent`` method is used to load the full content of the text file.
        '''
        
        tableContent = HandyTools.getTextFileContent (textFileNameAndPath, stripLineFromBlanks = True)

        if tableContent:
        
            if (endOfHeaderString + "\n").strip () in tableContent:
            
                iHeaderEnd = tableContent.index ( (endOfHeaderString + "\n").strip () )
            
            
            elif endOfHeaderString.strip () in tableContent:
            
                iHeaderEnd = tableContent.index ( endOfHeaderString.strip () )
            
            
            #---------- If file has no "C_END" marker, then assume the data start at the top.
            else:

                print ('')
                print ('---ATTENTION---')
                print (' From HandyTools.readTable: ')
                print('   File {} does not contain {} marker: assuming data only'.format ( textFileNameAndPath, endOfHeaderString.strip () ) )
                print('')
                
                iHeaderEnd = -1


            # Store each line of the header in the  tableContentHeader  string list.
            #  Do not store empty lines.
            tableContentHeader = []
            if iHeaderEnd >= 0:
            
                tableContentHeader = [ headerLine.strip ()  for headerLine in tableContent [0:iHeaderEnd]  if headerLine.strip () ]

            
            # Go through each data line and extract each element of the line as a string and if possible as a number.
            numberOfLines = len (tableContent) 
            numberOfDataLines = numberOfLines - iHeaderEnd - 1                                   
            if numberOfDataLines:

                numberOfEntriesPerLine = len ( tableContent [iHeaderEnd + 1].split (separatorString) )
                
                tableContentDataNumbers = [ []  for iLine in range (numberOfEntriesPerLine) ]
                tableContentDataStrings = [ []  for iLine in range (numberOfEntriesPerLine) ]
                       
                # Go through each element on the line: elements are based on the  separatorString
                for iTableContent in range (iHeaderEnd + 1, numberOfLines):
                
                    for iDataValueString, dataValuesString in enumerate ( tableContent [iTableContent].split (separatorString) ):

                        tableContentDataStrings [iDataValueString].append ( dataValuesString.strip () )  

                        try:
                            
                            if '.' in dataValuesString:
                            
                                tableContentDataNumbers [iDataValueString].append ( float (dataValuesString) ) 
                                
                            else:
                            
                                tableContentDataNumbers [iDataValueString].append ( int (dataValuesString) )

                        except:
                        
                            tableContentDataNumbers [iDataValueString].append (np.nan)

            
            tableContentDataNumbers =  [ np.asarray (tableContentDataNumbers [iColumn], dtype = type (tableContentDataNumbers [iColumn][0]) )  for iColumn in range (numberOfEntriesPerLine) ]                                                
            
            print ('')
            print ('-------------')
            print (' From HandyTools.readTable: ')
            print ( '  file {} has been loaded with'.format (textFileNameAndPath) )
            print ( '  {} data lines and {} columns'.format (numberOfDataLines, numberOfEntriesPerLine) )                

            
        else:
        
            print ('')
            print ('---WARNING---')
            print (' From HandyTools.readTable: ')
            print ( '  file {} does not contain any readable data.'.format (textFileNameAndPath) )                

            
        return tableContentDataNumbers, tableContentDataStrings, tableContentHeader
                

                   
    # Save content (list, dictionary, ...) to a numpy file with a custom extension.
    @staticmethod
    def saveContentToNumpyWithCustomExtension (contentToSave, fileName, extensionWithoutDot, overWrite = False):
        '''
        :param contentToSave: any content that the user wants to save to the numpy-format file.
        :type: anything!
        
        :param fileName: the file name of the file to be saved.
        :type fileName: str
 
        :param extensionWithoutDot: the extension of the file name to be saved.
        :type extensionWithoutDot: str
        
        :param overWrite: if set to ``True``, then overwrite any file with the same file name that might already exist.
        :type overWrite: bool; default = False
        
        :return: two booleans, the first for file saved successfully, the second for if the file already existed.
        :rtype: bool, bool
        
        **Description:**
        Save content (list, dictionary, ...) to a numpy file with a custom extension ``extensionWithoutDot``. If the file already exists, then it is overwritten only 
        if the used set the ``overWrite`` boolean to ``True``.
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
    @staticmethod
    def getDateAndTime ():
        '''
        :return: datetime object.
        :rtype: datetime object.

        **Description:**
        Obtain the date and time of the moment of the call of this function. The result is returned as a datetime object.
        This object has several attributes, for example *year*, *month*, *day*, *hour*, *minute* and *second*. Use ``dir`` for a full list of all attributes.
        The ``dataandtime`` module is used in this function.
        '''
        
        # This try - except loop is needed because apparently between python versions the date module has changed some of its structure
        try:

            dateAndTime = datetime.now ()

        except:

            dateAndTime = datetime.datetime.now ()

                
        return dateAndTime
 


    # Get the date and time now and return in a string format.
    @staticmethod
    def getDateAndTimeString (includeYMD = True, dateFormat = 'YMD', includeHMS = True):
        '''
        :param includeYMD: if set to ``True`` (default), then return the year, month and date.
        :type includeYMD: bool; default = True
        
        :param dateFormat: three formats are option: YMD (default), DMY and MDY.
        :type dateFormat: str
 
        :param includeHMS: if set to ``True`` (default), then return the hour, minute and second.
        :type includeHMS: bool; default = True
        
        :return: a string that contains the date and time in the format specified by the user.
        :rtype: str
        
        
        **Description:**
        Obtain the date and time at the moment of the function call and return it as a string. The year month and date order can be specified by the user.
        Per default the year month and date are returned with the hours, minutes and seconds. for example ``2024-02-05 at 10:44:23``.
        The user can opt to only return the year, month and date in three different formats (YMD, DMY, or MDY) or only the hours, minutes and seconds.
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
    @staticmethod
    def getHMSFromTotalNumberOfSeconds (numberOfSecondsTotal = 0, numberOfDigitsAccuracy = 4):
        '''
        :param numberOfSecondsTotal: the number of seconds to convert to hours, minutes and seconds.
        :type numberOfSecondsTotal: float; default = 0
 
        :param numberOfDigitsAccuracy: number of digits after the comma (default = 4).
        :type numberOfDigitsAccuracy: int; default = 4
        
        :return: two lists, the first has the number of hours, number of minutes contained in the ``numberOfSecondsTotal`` given by the used, the seconds has the strings of these three values.
        :rtype: list [float, float, float], list [str, str, str]
        
        
        **Description:**
        Calculate the number of hours, minutes and seconds that are in a given number of seconds. 
        '''
    
        numberOfHours = int (numberOfSecondsTotal / 3600)
        numberOfMinutes = int ( (numberOfSecondsTotal - numberOfHours * 3600) / 60 )
        
        accuracyFactor = 10**numberOfDigitsAccuracy
        numberOfSeconds = int ( ( numberOfSecondsTotal - numberOfHours * 3600 - numberOfMinutes * 60 ) * accuracyFactor ) / accuracyFactor

        return [numberOfHours, numberOfMinutes, numberOfSeconds], \
               [ '{:02d}'.format (numberOfHours), '{:02d}'.format (numberOfMinutes), '{:7.4f}'.format (numberOfSeconds) ]                 


    
    # Calculate the total number of seconds from an Hours Minutes Seconds list.
    @staticmethod
    def getTotalNumberOfSecondsfromHMS ( HMS = [0,0,0] ):
        '''
        :param HMS: list of the number of hours, minutes and seconds to convert to number of seconds.
        :type HMS: list [int or str, int or str, float or str]
        
        :return: number of seconds contained in the number of hours, minutes and seconds given by the user.
        :rtype: float
               
        **Description:**
        Calculate the total number of seconds from a list with the number of hours, minutes and seconds given by the user.
        Note that each of the elements of the list with the hours, minutes and seconds can be given as an int (float) or string.   
        '''
        
        for iElement, element in enumerate (HMS):
        
            if type (element) == str:
            
                HMS [iElement] = float (element)

            
        return HMS [0] * 3600 + HMS [1] * 60 + HMS [2]



    # Determine and print execution of a piece of code.
    @staticmethod
    def getRunTime (functionName, startTime = None, indent = 0, printResult = True):
        '''
        :param functionName: a string with the name of the function or part of the code that is being speed checked.
        :type functionName: str

        :param startTime: the start time with respect to which the time of execution is measured. 
        :type startTime: start time in seconds after start of epoch; default = None

        :param indent: number of indents for printing, each indent is three empty spaces.
        :type indent: int; default = 0

        :param printResult: if ``True`` (default), then print a statement with the result.
        :type printResult: bool; default = True
        
        :return: if ``startTime`` is ``None`` (default), then return the time of the call of this function in number of seconds since the start of the epoch (1 January 1970 at 00:00:00 UTC). If ``startTime`` is not ``None`` and  ``printResult`` is set to ``False``, then return the time passed between ``startTime`` and the time of the current call to this function in seconds.
        :rtype: float
        
        **Description:**
        With this function, the speed of execution of a piece of code can be evaluated. At the start of the (piece of) code, the start time can be retrieved by a call to the this function without ``startTime`` specified. The function returns the time at that moment (in number of seconds since the start of the epoch, which is defined as 1 January 1970 at 00:00:00 UTC). The user can specify the ``functionName`` string to print this information on the screen. After the execution of the (piece of) code, this function can be called again with the ``startTime`` retrieved earlier and the resulting number of seconds of the execution will be printed.
        '''
        
        # Get the start time of the function's execution. This is the default.
        if not startTime:
        
            if printResult:
            
                print ( ' ' * indent * 3 + ' --> Start time of {} retrieved'.format (functionName) )
        
            return time.time ()
        

        # If  printResult  is True, then print the end time, given the  startTime.        
        elif printResult:
        
            print ( ' ' * indent * 3 + ' --> Run time for {} = {:10.6f}s'.format (functionName, time.time () - startTime) )
        
        # If  printResult  is False, then return the end time, given the  startTime.
        else:
        
            return time.time () - startTime
        
