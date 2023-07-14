import numpy
from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("flow_class.pyx", compiler_directives={'linetrace': True}, annotate=True),
    include_dirs=[numpy.get_include()]
)
