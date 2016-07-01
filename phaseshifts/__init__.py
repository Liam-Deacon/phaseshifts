from pkg_resources import get_distribution, DistributionNotFound
import os

__author__ = 'Liam Deacon'
__author_email__ = 'liam.m.deacon@gmail.com'
__project__ = 'phaseshifts'
__version__ = '0.1.6'  # required for initial installation

try:
    filename = os.path.join(os.path.dirname(__file__), 'LICENSE.txt')
    with open(filename, 'r') as f:
        __license__ = f.readlines()[0].strip('\n')
    del(filename)
except:
    __license__ = 'unknown'

try:
    dist = get_distribution(__project__)
    egg_info = os.path.join(dist.location, 'EGG-INFO', 'PKG-INFO')
    with open(egg_info, 'r') as f:
        pkg_info = f.read()
    __version__ = dist.version
except IOError:
    pass
except DistributionNotFound:
    VERSION = __project__ + '-' + '(local)'
else:
    VERSION = __project__ + '-' + __version__
finally:  
    # clean up
    try:
        del(dist)
        del(egg_info)
        del(pkg_info)
    except:
        pass
