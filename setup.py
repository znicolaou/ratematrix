from distutils.core import setup, Extension
setup(name='rlist', version='1.0',  \
      ext_modules=[Extension('rlist', ['rlistmodule.c'])])
