# distutils: language = c++

import numpy as np
import cython

cdef extern from "FilterToolsCPPCore.cpp":

    pass


cdef extern from "FilterToolsCPPCore.h":
    
    cdef cppclass FilterToolsCPPCore:
    
        FilterToolsCPPCore () except+
    
        void passAverageFilter ( 
            int *, #1
            float *, #2
            int, #3
            int #4
        )


def passAverageFilterPYtoCPP ( 
    listOfNumbers,
    widthOfWindow ):


    cdef FilterToolsCPPCore FilterToolsCPPCore


    # Make sure all arrays passed into the function from the Python code are stored contiguously and are integers.
    listOfNumbers = np.ascontiguousarray (listOfNumbers, dtype = np.int32)
    cdef int [::1] listOfNumbers_view = listOfNumbers
    
    cdef int numberOfElements = len (listOfNumbers)

    listOfNumbersFiltered = np.ascontiguousarray ( np.zeros (numberOfElements, dtype = np.single) )
    cdef float [::1] listOfNumbersFiltered_view = listOfNumbersFiltered

            
    # Call the C++ core function.
    FilterToolsCPPCore.passAverageFilter ( 
        &listOfNumbers_view [0], #1
        &listOfNumbersFiltered_view [0], #2
        numberOfElements, #3        
        widthOfWindow #4
    )


    # Return the filter list of numbers.
    return listOfNumbersFiltered
           
           