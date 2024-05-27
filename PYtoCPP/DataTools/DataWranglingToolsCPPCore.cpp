// Standard includes.
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

// Custom includes.
#include "DataWranglingToolsCPPCore.h"



// Default constructor
DataWranglingToolsCPPCore::DataWranglingToolsCPPCore () {}

// Destructor
DataWranglingToolsCPPCore::~DataWranglingToolsCPPCore () {};


void DataWranglingToolsCPPCore::getSegmentSpecsFromDataValues ( 
    short dataValues [1], //1
    unsigned int numberOfDataValues, //2
    unsigned int segmentStartIndices [1], //3
    unsigned int& numberOfSegments, //4
    int segmentAmplitudes [1], //5
    float segmentSlopes [1], //6
    unsigned int segmentDurations [1], //7
    unsigned int segmentStartIndicesNegative [1], //8
    unsigned int& numberOfSegmentsNegative, //9
    unsigned int& iSteepestNegativeSlopeSegment, //10  
    unsigned int& iSegmentStartIndicesSteepestNegativeSlope, //11 
    unsigned int segmentStartIndicesPositive [1], //12
    unsigned int& numberOfSegmentsPositive, //13
    unsigned int& iSteepestPositiveSlopeSegment, //14 
    unsigned int& iSegmentStartIndicesSteepestPositiveSlope //15
)                                                        
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
    float steepestNegativeSlope = 0;
    iSteepestNegativeSlopeSegment = 0;
    iSegmentStartIndicesSteepestNegativeSlope = 0;

    float steepestPositiveSlope = 0;
    iSteepestPositiveSlopeSegment = 0;
    iSegmentStartIndicesSteepestPositiveSlope = 0;
      

    int iSegmentBefore = numberOfSegments;
    // Go through the  dataValues  list.
    for (unsigned int iSample = 1; iSample < numberOfDataValues - 1; iSample++){
 
        deltaValue = dataValues [iSample + 1] - dataValues [iSample];
        
        // Determine the closest segment back in time that has a non-zero delta value to compare to the current delta value.
        iSegmentBefore = numberOfSegments;
        if (deltaValue == 0)
        {
        
          while ( iSegmentBefore >= 0 && segmentAmplitudes [numberOfSegments] == 0 )
          
              iSegmentBefore--;
        
        };
        
          
        // The "=" in "segmentAmplitudes [iSegmentBefore] >=" and "segmentAmplitudes [iSegmentBefore] <=" is to deal with the situation when the 
        //  list of data points starts as a flat line, i.e. the first segments have a delta of zero.
        if ( ( segmentAmplitudes [iSegmentBefore] >= 0 && deltaValue >= 0 ) || (segmentAmplitudes [iSegmentBefore] <= 0 && deltaValue <= 0) )
        {
        
            segmentAmplitudes [numberOfSegments] += deltaValue;
            numberOfSamplesInSegment++;         
        
        }

        // A new segment has started.  
        else
        {
                   
            segmentDurations [numberOfSegments] = numberOfSamplesInSegment;
            segmentSlopes [numberOfSegments] = static_cast <float> ( segmentAmplitudes [numberOfSegments] ) / segmentDurations [numberOfSegments];
  
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



void DataWranglingToolsCPPCore::getAverageVarAndSD (
    double dataValues [1], //1
    int numberOfValues, //2
    double& averageValue, //3
    double& standardDeviation, //4
    double& variance //5
)
{
    
    // Calculate the average value.
    double sumOfDataValues = 0.;
    for (int iValue = 0; iValue < numberOfValues; iValue++)
      
        sumOfDataValues += dataValues [iValue];
    
    averageValue = sumOfDataValues / numberOfValues;

    // Calculate the variance and the standard deviation.
    sumOfDataValues = 0.;
    double differenceSquared;
    for (int iValue = 0; iValue < numberOfValues; iValue++)
    {
    
        differenceSquared = dataValues [iValue] - averageValue;
        sumOfDataValues += differenceSquared * differenceSquared;
    
    }
            
    variance = sumOfDataValues / numberOfValues;
    standardDeviation = sqrt (variance);  

}



void DataWranglingToolsCPPCore::getMedianAndQuantiles (
    double dataValues [1], //1
    int numberOfValues, //2
    double& medianValue, //3
    float lowerQuantile, //4
    double& lowerQuantileValue, //5
    float upperQuantile, //6
    double& upperQuantileValue //7
)
{
 
    std::vector <double> dataValuesSorted (numberOfValues);    
    for (int iValue = 0; iValue < numberOfValues; iValue++)
    
        dataValuesSorted [iValue] = dataValues [iValue];
        

    std::sort ( dataValuesSorted.begin (), dataValuesSorted.end () );    
    if (numberOfValues % 2)
    
        medianValue = dataValuesSorted [numberOfValues / 2];
        
    else
    
        medianValue = ( dataValuesSorted [numberOfValues / 2] + dataValuesSorted [numberOfValues / 2 - 1] ) / 2;

    lowerQuantileValue = getQuantileValue (dataValuesSorted, numberOfValues, lowerQuantile);  
    upperQuantileValue = getQuantileValue (dataValuesSorted, numberOfValues, upperQuantile);  
            

}


double DataWranglingToolsCPPCore::getQuantileValue (
    std::vector <double> dataValuesSorted,
    int numberOfValues,
    float quantile
)
{

    float virtualIndex = quantile * (numberOfValues - 1);
    int iVirtualIndex = static_cast <int> (virtualIndex);
    float fractionVirtualIndex = virtualIndex - iVirtualIndex;


    if (iVirtualIndex == numberOfValues - 1)

        return dataValuesSorted [iVirtualIndex];
    
    
    return dataValuesSorted [iVirtualIndex] +  (dataValuesSorted [iVirtualIndex + 1] - dataValuesSorted [iVirtualIndex]) * fractionVirtualIndex;


}

