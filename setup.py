from distutils.core import setup, Extension
import os
import numpy

os.environ["CC"] = "g++-6"
os.environ["CXX"] = "g++-6"
#os.environ["CFLAGS"] = "-fopenmp"
setup(name='rlist', version='1.0', ext_modules=[Extension('rlist', ['rlistmodule.c'], language="c++", include_dirs=[numpy.get_include()])], requires=['numpy'])
