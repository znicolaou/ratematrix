from distutils.core import setup, Extension
import os
import numpy

os.environ["CC"] = "g++"
setup(name='rlist', version='1.0', ext_modules=[Extension('rlist', ['rlistmodule.c'], language="c++", include_dirs=[numpy.get_include()], extra_compile_args=["-fopenmp"], extra_link_args=["-fopenmp"])], requires=['numpy'])
