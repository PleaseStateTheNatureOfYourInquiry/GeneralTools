**GeneralTools** is a set of functions that come in handy in many instances when doing analysis of any kind of digital data and when working with data files. 

There are currently two different pseudo classes in **GeneralTools**, which are stand-alone (they do not depend on each other): **HandyTools** and **DataTools**. 
The class **PyQtInformationDialog**  is only useful if you develop GUI with PyQt. You can happily ignore it if you do not.

In order to use these classes, the Python session must know their paths on your machine, so you need to tell Python via a Python startup file and/or in your .zprofile file (MAC) or other where to look. In Python language it is appending the following path to the ```sys.path``` variable:

  ```
  sys.path.append ('/SomeWhereOnYourMachine/GeneralTools/DataTools')
  sys.path.append ('/SomeWhereOnYourMachine/GeneralTools/PYtoCPP/DataTools')
  sys.path.append ('/SomeWhereOnYourMachine/GeneralTools/HandyTools')
  sys.path.append ('/SomeWhereOnYourMachine/GeneralTools/PyQtInformationDialog')
  ```

In the .zprofile file on a MAC it means adding the following lines:

```
export PYTHONPATH="/SomeWhereOnYourMachine/GeneralTools/DataTools:$PYTHONPATH"
export PYTHONPATH="/SomeWhereOnYourMachine/GeneralTools/PYtoCPP/DataTools:$PYTHONPATH"
export PYTHONPATH="/SomeWhereOnYourMachine/GeneralTools/HandyTools:$PYTHONPATH"
export PYTHONPATH="/SomeWhereOnYourMachine/GeneralTools/PyQtInformationDialog:$PYTHONPATH"
```

**HandyTools** is purely written in Python. It is a collection of small functions, that are ... well **handy**, as the name suggests (at least they are for me). I developed these over the years, as I was working on different projects. This set is pretty stable and I do not make many updates or changes to it. 

**DataTools** is written in Python and some of the functions have C++ bindings, via Cython. This is a more recent class and I am more actively adding and updating code. 
The Cython and C++ code and the compiled libraries are in the `./PYtoCPP/DataTools` subfolder. The compilation has been done for Python 3.11 (both Mac OS and Windows), which I use on my machine. On Mac, it might be necessary to have XCode installed before you can run it. If you use a different version of Python, then you must compile the two classes **DataWranglingToolsPYtoCPP** and **FilterToolsPYtoCPP** using

  ```
  > cd /SomeWhereOnYourMachine/GeneralTools/PYtoCPP/DataTools
  > python DataWranglingToolsPYtoCPP_setup.py build_ext --inplace
  > python FilterToolsPYtoCPP_setup.py build_ext --inplace
  ```
Note that you might have to install Cython and/or some other compilers (and/or XCode on Mac) before you can compile, it depends on your commputer's setup, and Python distribution and packages.






