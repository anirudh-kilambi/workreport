
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system


# compile the fortran modules without linking
fortran_mod_comp = 'gfortran cstintrp_py.f90 -c -o cstintrp_py.o -O3 -fPIC'
print (fortran_mod_comp)
system(fortran_mod_comp)

ext_modules = [Extension(# module name:
                         'cstintrp_py',
                         # source file:
                         ['cstintrp_py.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['cstintrp_py.o'],
                         libraries = ['gfortran'])]

setup(name = 'cstintrp_py',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = cythonize(ext_modules, compiler_directives={'language_level' : "3"}))