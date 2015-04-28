#!/usr/bin/env python
# encoding: utf-8
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2013-2015 Liam Deacon                             #
#                                                                            #
# License: MIT License                                                       #
#                                                                            #
# Permission is hereby granted, free of charge, to any person obtaining a    #
# copy of this software and associated documentation files (the "Software"), #
# to deal in the Software without restriction, including without limitation  #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,   #
# and/or sell copies of the Software, and to permit persons to whom the      #
# Software is furnished to do so, subject to the following conditions:       #
#                                                                            #
# The above copyright notice and this permission notice shall be included in #
# all copies or substantial portions of the Software.                        #
#                                                                            #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    #
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        #
# DEALINGS IN THE SOFTWARE.                                                  #
#                                                                            #
##############################################################################

try:
    from setuptools import find_packages
except ImportError:
    from distutils.core import find_packages

from numpy.distutils.core import Extension, setup
from numpy.distutils import fcompiler
from numpy.distutils import ccompiler
from numpy.distutils.fcompiler.intel import BaseIntelFCompiler
from numpy.distutils.fcompiler.gnu import GnuFCompiler
from tempfile import gettempdir
from abc import abstractmethod, ABCMeta
from glob import glob

try:
    from Cython.Build import BuildExecutable, cythonize, Cythonize
except ImportError:
    pass

import sys
import os

try:
    import py2exe
except:
    pass

if len(sys.argv) == 1:
    sys.argv.append('install')

phsh_lib = os.path.join('phaseshifts', 'lib')


class Builder(object):
    def __init__(self, name, sources,
                 include_dirs=None,
                 define_macros=None,
                 undef_macros=None,
                 library_dirs=None,
                 libraries=None,
                 runtime_library_dirs=None,
                 export_symbols=None,
                 extra_objects=None,
                 extra_compile_args=None,
                 extra_link_args=None,
                 depends=None,
                 language=None,
                 module_dirs=None,
                 debug=False,
                 ccompiler=ccompiler.new_compiler()):
        self.name = name
        self.sources = sources or []
        self.include_dirs = include_dirs or []
        self.define_macros = define_macros or []
        self.undef_macros = undef_macros or []
        self.library_dirs = library_dirs or []
        self.libraries = libraries or []
        self.runtime_library_dirs = runtime_library_dirs or []
        self.export_symbols = export_symbols or []
        self.extra_objects = extra_objects or []
        self.extra_compile_args = extra_compile_args or []
        self.extra_link_args = extra_link_args or []
        self.depends = depends or []
        self.language = language
        self.module_dirs = module_dirs or []
        self.debug = debug
        self.ccompiler = ccompiler

    def _clean(self, build_dir='.'):
        [os.remove(f) for f in glob(os.path.join(build_dir, '*.o'))]


class CBuilder(Builder):
    def __init__(self, name, sources,
                 include_dirs=None,
                 define_macros=None,
                 undef_macros=None,
                 library_dirs=None,
                 libraries=None,
                 runtime_library_dirs=None,
                 export_symbols=None,
                 extra_objects=None,
                 extra_compile_args=None,
                 extra_link_args=None,
                 depends=None,
                 language=None,
                 module_dirs=None,
                 debug=False,
                 ccompiler=ccompiler.new_compiler()):
        Builder.__init__(name, sources,
                         include_dirs=include_dirs, 
                         define_macros=define_macros, 
                         undef_macros=undef_macros, 
                         library_dirs=library_dirs, 
                         libraries=libraries, 
                         runtime_library_dirs=runtime_library_dirs, 
                         export_symbols=export_symbols, 
                         extra_objects=extra_objects, 
                         extra_compile_args=extra_compile_args, 
                         extra_link_args=extra_link_args, 
                         depends=depends, 
                         language=language, 
                         module_dirs=module_dirs, 
                         debug=debug, 
                         ccompiler=ccompiler)
        
        
class FortranBuilder(Builder):
    ''' 
    Class for building Fortran shared libraries and executables using NumPy.
    
    Notes
    -----
    
    Unlike the numpy.distutils.core.Extension() class, the FortranBuilder() 
    class allows pure Fortran libraries which are not wrapped with f2py 
    (useful for complicated sources/builds where f2py fails). 
    
    This class should only be used for experts as there are no/limited checks 
    for the integrity of the initialisation arguments.
    
    '''
    def __init__(self, name, sources,
                 include_dirs=None,
                 define_macros=None,
                 undef_macros=None,
                 library_dirs=None,
                 libraries=None,
                 runtime_library_dirs=None,
                 export_symbols=None,
                 extra_objects=None,
                 extra_compile_args=None,
                 extra_link_args=None,
                 depends=None,
                 language=None,
                 module_dirs=None,
                 extra_f77_compile_args=None,
                 extra_f90_compile_args=None,
                 debug=False,
                 fcompiler=fcompiler.new_fcompiler(requiref90=True)):
        Builder.__init__(self, name, sources, 
                         include_dirs=include_dirs, 
                         define_macros=define_macros, 
                         undef_macros=undef_macros, 
                         library_dirs=library_dirs, 
                         libraries=libraries, 
                         runtime_library_dirs=runtime_library_dirs, 
                         export_symbols=export_symbols, 
                         extra_objects=extra_objects, 
                         extra_compile_args=extra_compile_args, 
                         extra_link_args=extra_link_args, 
                         depends=depends, 
                         language=language, 
                         module_dirs=module_dirs, 
                         debug=debug)
        self.fcompiler = fcompiler
        self.extra_f77_compile_args = extra_f77_compile_args or []
        self.extra_f90_compile_args = extra_f90_compile_args or []
        self.fcompiler._is_customised = True
        
    def _clean(self, build_dir='.'):
        '''
        Cleans up compiled objects in the specified build directory
        '''
        Builder._clean(self, build_dir=build_dir)
        [os.remove(f) for f in glob(os.path.join(build_dir, '*.mod'))]
    
    def _compile_sources(self):
        ''' 
        Compiles source files and returns list of compiled object files 
        '''
        print("Compiling sources: {srcs}...".format(srcs=self.sources))
        cwd = os.path.abspath(os.path.curdir) 
        sources = [src if os.path.isabs(src) 
                   else os.path.join(os.path.abspath(cwd), src)
                   for src in self.sources]
        os.chdir(gettempdir())
        self._clean(gettempdir())
        objects = self.fcompiler.compile(sources, 
                                         output_dir=gettempdir(), 
                                         macros=self.define_macros, 
                                         include_dirs=self.include_dirs, 
                                         debug=self.debug, 
                                         extra_preargs=self.extra_compile_args, 
                                         depends=self.depends)
        os.chdir(cwd)
        return objects
        
    def make_lib(self, output_dir='.'):
        ''' 
        Creates a FORTRAN shared library for use with ctypes and places it 
        in output_dir. 
        '''
        objects = self._compile_sources() 
        linker_flags = set(self.extra_link_args + 
                           self.fcompiler.get_flags_linker_so())
        output_file = str('lib' + self.name + 
                          fcompiler.get_shared_lib_extension(False))
        
        print("Creating '{}' Fortran shared library...".format(output_file))
        self.fcompiler.link(self.fcompiler.SHARED_LIBRARY,
                            objects, 
                            output_filename=output_file,
                            output_dir=output_dir, 
                            libraries=self.libraries, 
                            library_dirs=self.library_dirs, 
                            runtime_library_dirs=self.runtime_library_dirs, 
                            export_symbols=self.export_symbols, 
                            debug=self.debug, 
                            extra_preargs=linker_flags,
                            target_lang=self.language)
        
    def make_exe(self, output_dir='.'):
        ''' 
        Creates an executable and places it into output_dir 
        '''
        objects = self._compile_sources()
        exe_flags = set(self.extra_link_args + 
                        self.fcompiler.get_flags_linker_exe())
        output_file = (self.name + '.exe' 
                       if str(sys.platform).startswith('win') else self.name)
        print("Creating '{}' executable...".format(output_file))
        self.fcompiler.link(self.fcompiler.EXECUTABLE,
                            objects, 
                            output_filename=output_file,
                            output_dir=output_dir, 
                            libraries=self.libraries, 
                            library_dirs=self.library_dirs, 
                            runtime_library_dirs=self.runtime_library_dirs, 
                            export_symbols=self.export_symbols, 
                            debug=self.debug, 
                            extra_preargs=exe_flags,
                            target_lang=self.language)
        return output_file

if 'install' in sys.argv or 'build' in sys.argv or 'build_ext' in sys.argv:
    default_fcompiler = fcompiler.new_fcompiler(requiref90=True)
    compile_args = []
    if isinstance(default_fcompiler, GnuFCompiler):
        compile_args = ['-O2', '-msse4.1', '-march=native', 
                        '-finline-functions', '-fbacktrace', '-fopenmp']
    elif isinstance(default_fcompiler, BaseIntelFCompiler):
        compile_args = ['-O2', '-xSSE4.1', '-finline-functions', '-traceback']
 
    builder = FortranBuilder(name='EEASiSSS',
                             sources=[os.path.join(phsh_lib, 'EEASiSSS', 
                                                   'EEASiSSS_2015_03_28.f90')], 
                             extra_compile_args=compile_args)

    # create shared library
    builder.make_lib(output_dir=phsh_lib)
    
    # now add executable
    builder.sources += [os.path.join(phsh_lib, 'EEASiSSS', 
                                     'EEASiSSS_main.f90')]
    # builder.make_exe(output_dir=phsh_lib)

# build f2py extensions
f2py_exts = [Extension(name='phaseshifts.lib.libphsh',
                       extra_compile_args=['-fopenmp'],
                       extra_link_args=['-lgomp'],
                       sources=[os.path.join(phsh_lib, 'libphsh.f')]), 
             
             Extension(name='phaseshifts.lib.libhartfock',
                       extra_compile_args=[],
                       sources=[os.path.join(phsh_lib, 'EEASiSSS', 'hf.f90')]),
             ]

readme = os.path.join('phaseshifts', 'README.rst')
    
dist = setup(name='phaseshifts', 
             packages=find_packages(),
             version='0.1.6-dev',
             author='Liam Deacon',
             author_email='liam.deacon@diamond.ac.uk',
             license='MIT License',
             url='https://pypi.python.org/pypi/phaseshifts',
             description='Python package for calculating phase shifts '
                         'for LEED/XPD modelling',
             long_description=(open(readme).read() 
                               if os.path.exists(readme) else None),
             classifiers=['Development Status :: 4 - Beta',
                          'Environment :: Console',
                          # The end goal is to have Qt or other GUI frontend
                          'Environment :: X11 Applications :: Qt',  
                          'Intended Audience :: Science/Research',
                          'License :: OSI Approved :: MIT License',
                          'Operating System :: OS Independent',
                          'Programming Language :: Python',
                          'Topic :: Scientific/Engineering :: Chemistry',
                          'Topic :: Scientific/Engineering :: Physics',
                          ],
             keywords='phaseshifts atomic scattering muffin-tin diffraction',
             # recursive-include phaseshifts *.py *.pyw
             include_package_data=True,
             # If any package contains *.txt or *.rst files, include them:
             package_data={'': ['*.txt', '*.rst', '*.pyw', 'ChangeLog'],
                           'lib': ['lib/*.f', 'lib/*.c', 'lib/*.h', 
                                   'lib/*.dll', 'lib/*.so'],
                           'gui': ['gui/*.ui', 'gui/*.bat'],
                           'gui/res': ['gui/res/*.*']
                           },
             scripts=[os.path.join("phaseshifts", "phsh.py"), 
                      os.path.join("phaseshifts", "lib", "EEASiSSS", "hf.py"),
                      os.path.join("phaseshifts", "lib", "EEASiSSS", 
                                   "eeasisss.py")
                      ],
             # data_files = cython_exts,
             install_requires=['scipy >= 0.7', 
                               'numpy >= 1.3', 
                               'periodictable'],
             ext_modules=f2py_exts,
             #console=[os.path.join("phaseshifts", "phsh.py")],
             options={'py2exe': {'skip_archive': 1,
                                 'compressed': 0,  
                                 'bundle_files': 2, 
                                 'dist_dir': os.path.join("dist", "py2exe"),
                                 'excludes': ['tcl', 'bz2'],
                                 'dll_excludes': ['w9xpopen.exe', 
                                                  'tk85.dll', 
                                                  'tcl85.dll']
                                 }
                      },
             # zipfile = None
               
             )

sys.argv.append('build_ext')
sys.argv.append('--inplace')


# End of file
