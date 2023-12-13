# distutils: language = c++

import numpy as np
import cython

cdef extern from "HandyToolsCPPCore.cpp":

    pass


cdef extern from "HandyToolsCPPCore.h":
    
    void getListOfAmplitudeSegments ( short int *, unsigned int, 
                                      unsigned int *, unsigned int&,
                                      short int *, float *, unsigned int *,
                                      unsigned int *, unsigned int&, unsigned int&, unsigned int&,
                                      unsigned int *, unsigned int&, unsigned int&, unsigned int& )


def getListOfAmplitudeSegmentsPYtoCPP (dataValues):


    # Make sure the dataValues list is a NumPy array.
    if type (dataValues) == list:
    
        dataValues = np.asarray (dataValues, dtype = np.short)

    # Make sure the array is stored contiguously.
    if not dataValues.flags ['C_CONTIGUOUS']:
    
        dataValues = np.ascontiguousarray (dataValues)
    
    
    cdef Py_ssize_t numberOfDataValues = dataValues.shape [0]
    cdef short int [::1] dataValues_view = dataValues


    cdef unsigned int numberOfSegments = 0
    cdef unsigned int numberOfSegmentsNegative = 0
    cdef unsigned int numberOfSegmentsPositive = 0
    cdef unsigned int iSteepestNegativeSlopeSegment = 0
    cdef unsigned int iSegmentStartIndicesSteepestNegativeSlope = 0
    cdef unsigned int iSteepestPositiveSlopeSegment = 0
    cdef unsigned int iSegmentStartIndicesSteepestPositiveSlope = 0
    

    # Initialise the arrays to the length of the the  dataValues  array.
    segmentAmplitudes = np.ascontiguousarray ( np.zeros (numberOfDataValues, dtype = np.short) )
    cdef short int [::1] segmentAmplitudes_view = segmentAmplitudes

    segmentSlopes = np.ascontiguousarray ( np.zeros (numberOfDataValues, dtype = np.single) )
    cdef float [::1] segmentSlopes_view = segmentSlopes

    segmentDurations = np.ascontiguousarray ( np.zeros (numberOfDataValues, dtype = np.uintc) )
    cdef unsigned int [::1] segmentDurations_view = segmentDurations

    segmentStartIndices = np.ascontiguousarray ( np.zeros (numberOfDataValues, dtype = np.uintc ) )
    cdef unsigned int [::1] segmentStartIndices_view = segmentStartIndices
    
    segmentStartIndicesNegative = np.ascontiguousarray ( np.zeros (numberOfDataValues, dtype = np.uintc ) )
    cdef unsigned int [::1] segmentStartIndicesNegative_view = segmentStartIndicesNegative
    
    segmentStartIndicesPositive = np.ascontiguousarray ( np.zeros (numberOfDataValues, dtype = np.uintc ) )
    cdef unsigned int [::1] segmentStartIndicesPositive_view = segmentStartIndicesPositive
    
    
    # Call the C++ core function.
    getListOfAmplitudeSegments ( &dataValues_view [0], numberOfDataValues, 
                                 &segmentStartIndices_view [0], numberOfSegments,
                                 &segmentAmplitudes_view [0], &segmentSlopes_view [0], &segmentDurations_view [0],
                                 &segmentStartIndicesNegative_view [0], numberOfSegmentsNegative, 
                                 iSteepestNegativeSlopeSegment, iSegmentStartIndicesSteepestNegativeSlope,
                                 &segmentStartIndicesPositive_view [0], numberOfSegmentsPositive, 
                                 iSteepestPositiveSlopeSegment, iSegmentStartIndicesSteepestPositiveSlope )


    # Reduce the lengths of the arrays to the actuals lengths.
    segmentAmplitudes = segmentAmplitudes [0:numberOfSegments]
    segmentSlopes = segmentSlopes [0:numberOfSegments]
    segmentDurations = segmentDurations [0:numberOfSegments]
    segmentStartIndices = segmentStartIndices [0:numberOfSegments]
    segmentStartIndicesNegative = segmentStartIndicesNegative [0:numberOfSegmentsNegative]
    segmentStartIndicesPositive = segmentStartIndicesPositive [0:numberOfSegmentsPositive]
    
    
    return segmentStartIndices, numberOfSegments, \
           segmentAmplitudes, segmentSlopes, segmentDurations, \
           segmentStartIndicesNegative, numberOfSegmentsNegative, \
           iSteepestNegativeSlopeSegment, iSegmentStartIndicesSteepestNegativeSlope, \
           segmentStartIndicesPositive, numberOfSegmentsPositive, \
           iSteepestPositiveSlopeSegment, iSegmentStartIndicesSteepestPositiveSlope

    


    