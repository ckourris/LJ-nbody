"""
Use this code to compile the C++ acceleration library "accelerate_lib"
Compilation is done using the command:
python3 compilec.py build_ext --inplace
Note: This requires Cython3.
"""
from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup( cmdclass={'build_ext': build_ext},
       ext_modules=[Extension("accelerate_lib",
       sources=["accelerate_module.pyx", "accelerate.cpp"],
       include_dirs=[numpy.get_include()],
       language='c++')]
)
