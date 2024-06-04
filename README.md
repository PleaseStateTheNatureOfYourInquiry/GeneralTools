**GeneralTools** is a set of functions that come in handy in many instances when doing analysis of any kind of digital data and when working with data files. 

There are currently two different pseudo classes in **GeneralTools**, which are stand-alone (they do not depend on each other): **HandyTools** and **DataTools**. 
The class **PyQtInformationDialog**  is only useful if you develop GUI with PyQt. You can happily ignore it if you do not.

In order to use these classes, the Python session must know their paths on your machine, so you need to tell Python via a Python startup file and/or in your .zprofile file (MAC) or other where to look. In Python language it is appending the path to the ```sys.path``` variable:

  ```
  sys.path.append ('/SomeWhereOnYourMachine/GeneralTools/DataTools')
  sys.path.append ('/SomeWhereOnYourMachine/GeneralTools/HandyTools')
  ```

In the .zprofile file on a MAC it means adding lines like these:

```
export PYTHONPATH="/SomePathOnYourMachine/GeneralTools/PYtoCPP/DataTools:$PYTHONPATH"
```

**HandyTools** is purely written in Python and  It is a collection of small functions, that come in ... **handy**, as the name suggests. I wrote and collected these over the years. This set is pretty stable (I think) and I do not make many updates or changes to it anymore. 

**DataTools** is written in Python, and some of the functions have with bindings to C++, via Cython. This is a more recent class and I am more actively adding and updating things here. 
The Cython and C++ code and the compiled libraries are in the `./PYtoCPP/DataTools` subfolder. The compilation has been done for Python 3.11, which I use on my machines. 
If you use a different version of Python, then you must compile the two classes **DataWranglingToolsPYtoCPP** and **FilterToolsPYtoCPP** using

  ```
  > cd /SomeWhereOnYourMachine/GeneralTools/PYtoCPP/DataTools
  > python DataWranglingToolsPYtoCPP_setup.py build_ext --inplace
  > python FilterToolsPYtoCPP_setup.py build_ext --inplace
  ```
Note that you might have to install Cython and/or some other compilers before you can compile, it depends on your Python distribution and packages.

Note that before instantiating **DataTools** the path to `/SomeWhereOnYourMachine/GeneralTools/PYtoCPP/DataTools` needs of course to be known by the Python session:

  ```
  sys.path.append ('/SomeWhereOnYourMachine/GeneralTools/PYtoCPP/DataTools')
  ```








