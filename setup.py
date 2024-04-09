from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

examples_extension = Extension(
    name="pyRAPIDD",
    sources=["pyRAPIDD.pyx"],
    library_dirs=["lib/source/build"],
    libraries=["RAPIDD", "gsl", "gslcblas"],
    include_dirs=["lib/source"]
)
setup(
    name="pyphys_consts",
    ext_modules=cythonize([examples_extension])
)
