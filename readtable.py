
from pathlib import Path
import numpy as np
import os
import time
from datetime import datetime



class ReadTable:

     """ 
     Read a table with or without a header, followed by data.
     Header ends at C_END and data lines start immediately afterwards
      or data lines start immediately at the top.
     Data is returned in the "data" instance variable along with the 
      n_data_lines and n_data_columns instance variables.
     "data" is a Numpy array of dimension  n_data_lines,n_data_columns
      
     When the "silent" variable is set, no information is printed.
     
     The use of os.linesep will ensure this works across different platforms.

     Maarten Roos
     v.20180307     
     """
 
 
     
     def __init__(self, file_name, silent=0):
     
          self.file_name = file_name
          self.silent = silent


#--------- Do not print this information if silent = 1

          if self.silent != 1:
               print('')
               print('ReadTable v.20180307')
               print('')


#---------- Check if the file exists

          file_check = Path(self.file_name)
          
          if file_check.is_file():
                    
               file_open = open (self.file_name,"r")
               file_content = file_open.readlines ()
               file_open.close ()


#---------- If file has a "C_END" marker, then determine at what line it is
#----------  the data is assumed to start right after the marker

               if "C_END\n" in file_content:
          
                    n_comment_lines = file_content.index("C_END\n")


               elif "C_END" in file_content:
          
                    n_comment_lines = file_content.index("C_END")


#---------- If file has no "C_END" marker, then assume the data start at the top
          
               else:

                    if self.silent != 1:
                         print('The file {:s} does not contain "C_END" marker: assuming data only'.format(self.file_name))
                         print('')

                    n_comment_lines = -1



#---------- read the data in the "data" instance variable


               file_open = open(self.file_name,"r")
               file_content = file_open.readlines()
               file_open.close()


#---------- delete any empty lines at the end of the file

               while file_content[len(file_content)-1] == os.linesep:
                    del file_content[-1]


               self.n_data_lines = len(file_content) - n_comment_lines - 1

               data_in = file_content[n_comment_lines+1].split()
               self.n_data_columns = len(data_in)


               if self.silent != 1:          
                    print('The file {:s} has {:d} entries with {:d} values'.format(self.file_name,self.n_data_lines,self.n_data_columns) )
                    print('Data returned in "data" instance variable')
                    print(' as well as "n_data_lines" and "n_data_columns" ')

          
               self.data = np.empty([self.n_data_lines,self.n_data_columns])
               
               for i in range( n_comment_lines+1,len(file_content) ):             
               
                    data_in = file_content[i].split()
               
                    for j in range(self.n_data_columns):

                         self.data[i-n_comment_lines-1,j] = \
                              float(data_in[j])  



#---------- In case the file is not present, print a warning message                              
          
          else:
          
               print('The file {:s} does not exist'.format(self.file_name))
               print('')


          return             
          
          
