
from setuptools import setup, Extension
from Cython.Build import cythonize


setup ( ext_modules = cythonize ( "FilterToolsPYtoCPP.pyx", 
                                  compiler_directives = {'language_level' : "3"} ) )

                                  
                                  
# from distutils.core import setup
# from distutils.extension import Extension
# from Cython.Build import cythonize
# 
# setup(
# 
#     ext_modules = cythonize (
#         Extension(
#             "ProcessorModulePYtoCPP", ["ProcessorModulePYtoCPP.pyx"],
#             extra_compile_args = ["-O3"],
#             language="c++",
#         ),
#     ),
# )