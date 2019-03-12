from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup( cmdclass={'build_ext': build_ext},
       ext_modules=[Extension("accelerate_module",
       sources=["accelerate_lib.pyx", "accelerate.cpp"],
       include_dirs=[numpy.get_include()],
       language='c++')]
)
