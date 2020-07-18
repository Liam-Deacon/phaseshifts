#!/usr/bin/env python
# encoding: utf-8
import os

#from setuptools import setup, find_packages
#import fix_setuptools_chmod
try:
    from setuptools import find_packages
except ImportError:
    from distutils.core import find_packages

from numpy.distutils.core import Extension, setup

#import phaseshifts
import sys, os
try:
    import py2exe
except:
    pass

if len(sys.argv) == 1:
    sys.argv.append('install')

# build f2py extensions
f2py_exts = [Extension(name='phaseshifts.lib.libphsh',
                      extra_compile_args=['-fopenmp'],
                      extra_link_args=['-lgomp'],
                      sources=[os.path.join('phaseshifts','lib','libphsh.f')])
            ]
    
dist = setup(
        name = 'phaseshifts',
		#packages=['phaseshifts', 'phaseshifts.gui', 'phaseshifts.lib', 'phaseshifts.contrib'],
        packages = find_packages(),
        version='0.1.6-dev',
        author='Liam Deacon',
        author_email='liam.deacon@lightbytestechnology.co.uk',
        license='MIT License',
        url='https://pypi.python.org/pypi/phaseshifts',
        description='Python-based version of the Barbieri/Van Hove phase '
        'shift calculation package for LEED/XPD modelling',
        long_description=None if not os.path.exists('README.rst') else open('README.rst').read(),
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
			'Environment :: X11 Applications :: Qt',  # The end goal is to have Qt or other GUI frontend
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Physics',
            ],
        keywords='phaseshifts atomic scattering muffin-tin diffraction',
        #recursive-include phaseshifts *.py *.pyw
        include_package_data = True,
        package_data = {
			# If any package contains *.txt or *.rst files, include them:
			'' : ['*.txt', '*.rst', '*.pyw', 'ChangeLog'],
			'lib' : ['lib/*.f', 'lib/*.c', 'lib/*.h'],
            'gui' : ['gui/*.ui', 'gui/*.bat'],
            'gui/res' : ['gui/res/*.*']
			},
        scripts = [ "phaseshifts/phsh.py" ],
        #data_files = cython_exts,
        install_requires = ['scipy >= 0.7', 'numpy >= 1.3', 'periodictable'],
		ext_modules = f2py_exts,
        console=[os.path.join("phaseshifts", "phsh.py")],
        # options={
            # 'py2exe': { 
                        # 'skip_archive':1,
                        # 'compressed':0,  
                        # 'bundle_files': 2, 
                        # 'dist_dir': os.path.join("dist", "py2exe"),
                        # 'excludes':['tcl', 'bz2'],
                        # 'dll_excludes':['w9xpopen.exe', 'tk85.dll', 'tcl85.dll']
                       # }
               # },
        #zipfile = None
               
)

sys.argv.append('build_ext')
sys.argv.append('--inplace')


# End of file
