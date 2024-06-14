

.. admonition:: GeneralTools has been developed by Maarten Roos-Serote

    | `ORCID <https://orcid.org/my-orcid?orcid=0000-0001-5001-1347>`_ (0000-0001-5001-1347) 
    | `Research Gate <https://www.researchgate.net/profile/Maarten-Roos-Serote>`_



GeneralTools for Scientists -  Documentation
=============================================

**GeneralTools** is a set of functions designed to help scientists working with any kind of digital data and (data) files. 
I have developed these functions over the years while working on a variety of different science projects. 
They are pieces of code that I kept reinventing and googling for, as I kept forgetting the exact details, while writing scripts to work with data and files.
By sharing I hope these functions can be useful to other people as well.
In the sections on these documentation pages, you can find the API's and descriptions of the functions.


There are currently two different pseudo classes in **GeneralTools**, which are stand-alone (they do not depend on each other): :ref:`HandyTools <handytoolssection>` and `DataTools <datatoolssection>`. 
The class **PyQtInformationDialog**  is only useful when developing GUIs with PyQt. Just ignore it if this is not your case.


I have tried hard to write them as general as possible, and to document their workings as clear as possible. 
Kindly let me know if there are things that are unclear, or can be done in a more efficient way.

On that last note, the more efficient way, I try to avoid difficult to read code. 
As an example, in general I prefer

.. code-block:: Python

    someList = []
    for iSomeOtherCounter in range ( len (aList) ):
    
        for iSomeCounter in range (100):
        
            if someCondition:
                
                someList.append ( someValue [iSomeCounter] )

as opposed to:

.. code-block:: Python

    someList = [ someValue [iSomeCounter]  for iSomeOtherCounter in range ( len (aList) )  for iSomeCounter in range (100)  if someCondition  ]


Also, I try to give good names to variables, so that the code can be read and understood *as a text*.
This sometimes leads to long variable names, but I prefer those over abbreviated names that I forget what they mean after a week.

If you have comments or suggestions, then please feel free to contact me. You are of course 100% free to make any chances to your clone as well.





.. toctree::
   :hidden:


   handytools
   datatools

   
