import os
import sys

from numpy.distutils.core import setup
from numpy.distutils.extension import Extension


# from distutils.core import setup
sys.argv.append('build_ext')
sys.argv.append('--inplace')

# we'd better have Cython installed, or it's a no-go
try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
except ImportError:
    sys.stderr.write("You don't seem to have Cython installed. Please get a "
                     "copy from www.cython.org and install it\n")
    sys.exit(1)

# build Cython extensions
for module in ['atorb.pyx', 'conphas.pyx', 'elements.pyx', 'leed.pyx', 'model.pyx']:
    setup(
        cmdclass={'build_ext': build_ext},
        ext_modules=cythonize([module]),
    )

# build f2py extensions
setup(
    ext_modules=[Extension(name='libphsh',
                           sources=[os.path.join('..', 'lib', 'libphsh.f')])],
)
