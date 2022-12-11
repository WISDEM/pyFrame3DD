# setup.py
# only if building in place: ``python setup.py build_ext --inplace``

import os
import sys
import platform
import glob
from setuptools import find_packages
from skbuild import setup

#os.environ['NPY_DISTUTILS_APPEND_FLAGS'] = '1'

#if os.name == 'nt':  # Windows.
#    extra_compile_args = ['/TC', '/D', 'ANSI'] # for msvs
#     # TODO: Not with Anaconda MINGW
#else:
#extra_compile_args = ''

#froot = 'pyframe3dd' + os.sep + 'src' + os.sep
#pyframeExt = Extension('pyframe3dd._pyframe3dd', sources=[froot+'py_HPGmatrix.c',
#                                               froot+'HPGutil.c',
#                                               froot+'NRutil.c',
#                                               froot+'coordtrans.c',
#                                               froot+'preframe.c',
#                                               froot+'py_eig.c',
#                                               froot+'py_frame3dd.c',
#                                               froot+'py_io.c',
#                                               froot+'py_main.c'])

setup(
    name='pyFrame3DD',
    version='1.2',
    description='Python bindings to Frame3DD',
    author='NREL WISDEM Team',
    author_email='systems.engineering@nrel.gov',
    #package_dir={'': 'src'},
    #py_modules=['pyframe3dd'],
    #package_data={'pyframe3dd': []},
    packages=['pyframe3dd'],
    license='Apache License, Version 2.0',
#    ext_modules=[pyframeExt],
    zip_safe=False,
    #cmake_install_dir='pyframe3dd',
    #cmake_source_dir='pyframe3dd' + os.sep + 'src' + os.sep,
)
