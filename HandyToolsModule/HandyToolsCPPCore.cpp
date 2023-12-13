#include "HandyToolsCPPCore.h"


// Default constructor
getListOfAmplitudeSegments::getListOfAmplitudeSegments () {}

// Overloaded constructor
getListOfAmplitudeSegments::getListOfAmplitudeSegments ( short dataValues [1], unsigned int numberOfDataValues,
                                                         unsigned int segmentStartIndices [1], unsigned int& numberOfSegments,
                                                         short segmentAmplitudes [1], float segmentSlopes [1], unsigned int segmentDurations [1], 
                                                         unsigned int segmentStartIndicesNegative [1], unsigned int& numberOfSegmentsNegative, 
                                                         unsigned int& iSteepestNegativeSlopeSegment,  unsigned int& iSegmentStartIndicesSteepestNegativeSlope, 
                                                         unsigned int segmentStartIndicesPositive [1], unsigned int& numberOfSegmentsPositive, 
                                                         unsigned int& iSteepestPositiveSlopeSegment , unsigned int& iSegmentStartIndicesSteepestPositiveSlope)                                                        
{

    int deltaValue;


// ATTENTION: Use these variable as an index counter, until the end of this function.
    numberOfSegments = 0;
    numberOfSegmentsNegative = 0;
    numberOfSegmentsPositive = 0;

    // Initialise the  segmentAmplitudes  first value with the first difference in  dataValues .
    segmentAmplitudes [numberOfSegments] = dataValues [1] - dataValues [0];
    segmentStartIndices [numberOfSegments] = 0;
    int numberOfSamplesInSegment = 1;
    
    // Initialise the steepest segment variables.
    int steepestNegativeSlope = 0;
    iSteepestNegativeSlopeSegment = 0;
    iSegmentStartIndicesSteepestNegativeSlope = 0;

    int steepestPositiveSlope = 0;
    iSteepestPositiveSlopeSegment = 0;
    iSegmentStartIndicesSteepestPositiveSlope = 0;
      

    // Go through the  dataValues  list.
    for (unsigned int iSample = 1; iSample < numberOfDataValues - 1; iSample++){
 
        deltaValue =  dataValues [iSample + 1] - dataValues [iSample];
        
        // If the difference between the next pair of dataValues has the same sign as the last  segmentAmplitudes  value, then add it.
        if ( ( segmentAmplitudes [numberOfSegments] > 0 and deltaValue >= 0 ) or (segmentAmplitudes [numberOfSegments] < 0 and deltaValue <= 0) )
        {
        
            segmentAmplitudes [numberOfSegments] += deltaValue;
            numberOfSamplesInSegment++;         
        
        }

        // A new segment has started.  
        else
        {
                   
            segmentDurations [numberOfSegments] = numberOfSamplesInSegment;
            segmentSlopes [numberOfSegments] = static_cast<float> ( segmentAmplitudes [numberOfSegments] ) / segmentDurations [numberOfSegments];

            // Reset the number of samples in the new segment to 1.
            numberOfSamplesInSegment = 1;
                       
            if ( segmentAmplitudes [numberOfSegments] < 0 )
            {
            
                // Check if the new segment is the steepest of the negative segments.
                if ( segmentSlopes [numberOfSegments] < steepestNegativeSlope )
                {

                    steepestNegativeSlope = segmentSlopes [numberOfSegments];
                    iSteepestNegativeSlopeSegment = segmentStartIndices [numberOfSegments];
                    iSegmentStartIndicesSteepestNegativeSlope = numberOfSegments;      
                
                }
                
                segmentStartIndicesNegative [numberOfSegmentsNegative] = segmentStartIndices [numberOfSegments];
                numberOfSegmentsNegative++;
            
            }  


            if ( segmentAmplitudes [numberOfSegments] > 0 )
            {
            
                // Check if the new segment is the steepest of the positive segments.
                if ( segmentSlopes [numberOfSegments] > steepestPositiveSlope )
                {

                    steepestPositiveSlope = segmentSlopes [numberOfSegments];
                    iSteepestPositiveSlopeSegment = segmentStartIndices [numberOfSegments];
                    iSegmentStartIndicesSteepestPositiveSlope = numberOfSegments;

                } 

                segmentStartIndicesPositive [numberOfSegmentsPositive] = segmentStartIndices [numberOfSegments];
                numberOfSegmentsPositive++;
                          
            }  

            // Initialise the new segment's amplitude value and store its start index.                
            numberOfSegments++;
            segmentAmplitudes [numberOfSegments] = deltaValue;
            segmentStartIndices [numberOfSegments] = iSample; 
            
        };
    
    }
    

    segmentDurations [numberOfSegments] = numberOfSamplesInSegment;
    segmentSlopes [numberOfSegments] = static_cast<float> ( segmentAmplitudes [numberOfSegments] ) / segmentDurations [numberOfSegments];

    // First check if the new segment is the steepest of the negative or positive slopes.
    if ( segmentAmplitudes [numberOfSegments] < 0 )
    {
    
        // Check if the new segment is the steepest of the negative segments.
        if ( segmentSlopes [numberOfSegments] < steepestNegativeSlope )
        {

            steepestNegativeSlope = segmentSlopes [numberOfSegments];
            iSteepestNegativeSlopeSegment = segmentStartIndices [numberOfSegments];
            iSegmentStartIndicesSteepestNegativeSlope = numberOfSegments;      
        
        }
        
        segmentStartIndicesNegative [numberOfSegmentsNegative] = segmentStartIndices [numberOfSegments];
    
    }  


    if ( segmentAmplitudes [numberOfSegments] > 0 )
    {

        // Check if the new segment is the steepest of the positive segments.
        if ( segmentSlopes [numberOfSegments] > steepestPositiveSlope )
        {

            steepestPositiveSlope = segmentSlopes [numberOfSegments];
            iSteepestPositiveSlopeSegment = segmentStartIndices [numberOfSegments];
            iSegmentStartIndicesSteepestPositiveSlope = numberOfSegments;

        } 

        segmentStartIndicesPositive [numberOfSegmentsPositive] = segmentStartIndices [numberOfSegments];        
                   
    }  


    // The actual number of segments is one more, because the initial value was set to zero so that it could be used as an index counter.
    numberOfSegments++;
    numberOfSegmentsNegative++;
    numberOfSegmentsPositive++;
        
}


// Destructor
getListOfAmplitudeSegments::~getListOfAmplitudeSegments () {};