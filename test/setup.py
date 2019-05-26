from distutils.core import setup, Extension
import os
import numpy

os.environ["CC"] = "g++"
# os.environ["CXX"] = "g++"
# os.environ["CFLAGS"] = "-fopenmp -lgomp"
setup(name='tlist', version='1.0', ext_modules=[Extension('tlist', ['tlistmodule.c'], language="c++", include_dirs=[numpy.get_include()], extra_compile_args=["-fopenmp"], extra_link_args=["-fopenmp"])], requires=['numpy'])
# setup(name='tlist', version='1.0', ext_modules=[Extension('tlist', ['tlistmodule.c'], language="c")])
