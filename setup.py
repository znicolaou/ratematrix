from distutils.core import setup, Extension
import numpy
print(numpy.get_include())

setup(name='rlist', version='1.0', ext_modules=[Extension('rlist', ['rlistmodule.c'],include_dirs=[numpy.get_include()]
)], requires=['numpy'])
