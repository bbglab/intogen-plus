from distutils.core import setup
from Cython.Build import cythonize

setup(name='Strongest path cython',
      ext_modules=cythonize("schulze_strongest_path_cython.pyx"))