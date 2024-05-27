#include <iostream>
#include <cmath>
#include <vector>
#include "FilterToolsCPPCore.h"



// Default constructor
FilterToolsCPPCore::FilterToolsCPPCore () {};

// Destructor
FilterToolsCPPCore::~FilterToolsCPPCore () {}; 


void FilterToolsCPPCore::passAverageFilter (
    // Electrogram
    int listOfNumbers [1], //1
    float listOfNumbersFiltered [1], //2
    int numberOfElements, //3
    int widthOfWindow //4
)
{

    // Filter the first  widthOfWindow  elements of the  listOfNumbers .
    for (int iElement = 0; iElement < widthOfWindow; iElement++)
    {
  
        for (int iWindow = 0; iWindow < iElement + widthOfWindow + 1; iWindow++)  
                  
            listOfNumbersFiltered [iElement] += listOfNumbers [iWindow];          


        listOfNumbersFiltered [iElement] /= (iElement + widthOfWindow + 1);
 
    };
  
    // Run the filter over all the elements.
    for (int iElement = widthOfWindow; iElement < numberOfElements - widthOfWindow; iElement++)
    {
  
        for (int iWindow = -widthOfWindow; iWindow < widthOfWindow + 1; iWindow++)   
             
            listOfNumbersFiltered [iElement] += listOfNumbers [iElement + iWindow];


        listOfNumbersFiltered [iElement] /= (2 * widthOfWindow + 1);

    }


    // Filter the last  widthOfWindow  elements of the  listOfNumbers .
    for (int iElement = numberOfElements - widthOfWindow; iElement < numberOfElements; iElement++)
    {
  
        for (int iWindow = -widthOfWindow; iWindow < numberOfElements - iElement; iWindow++)
        
            listOfNumbersFiltered [iElement] += listOfNumbers [iElement + iWindow];          
  

        listOfNumbersFiltered [iElement] /= (numberOfElements - iElement + widthOfWindow);
 
    };

}

