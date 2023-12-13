#ifndef HANDYTOOLSCPP_H
#define HANDYTOOLSCPP_H


class getListOfAmplitudeSegments {

    public:
    
        getListOfAmplitudeSegments ();
        getListOfAmplitudeSegments ( short dataValues [1], unsigned int numberOfDataValues, 
                                     unsigned int segmentStartIndices [1], unsigned int& numberOfSegments, 
                                     short segmentAmplitudes [1], float segmentSlopes [1], unsigned int segmentDurations [1],
                                     unsigned int segmentStartIndicesNegative [1], unsigned int& numberOfSegmentsNegative, 
                                     unsigned int& iSteepestNegativeSlopeSegment, unsigned int& iSegmentStartIndicesSteepestNegativeSlope,
                                     unsigned int segmentStartIndicesPositive [1], unsigned int& numberOfSegmentsPositive, 
                                     unsigned int& iSteepestPositiveSlopeSegment, unsigned int& iSegmentStartIndicesSteepestPositiveSlope );
        ~getListOfAmplitudeSegments ();


};



#endif

