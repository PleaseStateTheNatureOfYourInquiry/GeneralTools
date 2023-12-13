import numpy as np
import HandyToolsPYtoCPP


testList = np.asarray ( [1, 2, 5, 6, 5, 3, 1, -1, -4, -10, -3, 0, 1, 5, 6, 11, 12, 11, 7, 1, 0, -4, -11, -3, 2, 8, 12, 5, 0], dtype = np.short )


segmentStartIndices, numberOfSegments, \
segmentAmplitudes, segmentSlopes, segmentDurations, \
segmentStartIndicesNegative, numberOfSegmentsNegative, \
iSteepestNegativeSlopeSegment, iSegmentStartIndicesSteepestNegativeSlope, \
segmentStartIndicesPositive, numberOfSegmentsPositive, \
iSteepestPositiveSlopeSegment, iSegmentStartIndicesSteepestPositiveSlope = HandyToolsPYtoCPP.getListOfAmplitudeSegmentsPYtoCPP (testList)

print (numberOfSegments)
print (*segmentStartIndices)
print (*segmentAmplitudes)
print (*segmentSlopes)
print (*segmentDurations)


print ('negative')
print (numberOfSegmentsNegative)
print (*segmentStartIndicesNegative)
print ('steepest negative slope =', segmentSlopes [ iSegmentStartIndicesSteepestNegativeSlope ])
print (iSteepestNegativeSlopeSegment, iSegmentStartIndicesSteepestNegativeSlope)

print ('positive')
print (numberOfSegmentsPositive)
print (*segmentStartIndicesPositive)
print ('steepest positive slope =', segmentSlopes [ iSegmentStartIndicesSteepestPositiveSlope ])
print (iSteepestPositiveSlopeSegment, iSegmentStartIndicesSteepestPositiveSlope)

