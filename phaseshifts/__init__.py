from pkg_resources import get_distribution, DistributionNotFound
from os.path import dirname, join
import re

__project__ = 'phaseshifts'
__version__ = None  # required for initial installation
__author__ = 'Liam Deacon'

dictionary = {}

__author_email__ = 'liam.deacon@diamond.ac.uk'

try:
    filename = join(dirname(__file__), 'LICENSE.txt')
    with open(filename, 'r') as f:
        __license__ = f.read()
except:
    __license__ = 'unknown'

try:
    dist = get_distribution(__project__)
    egg_info = join(dist.location, 'EGG-INFO', 'PKG-INFO')
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


# prev = 0
# pos = 0
# while True:
#     match = pattern.search(text, pos)
#     if not match:
#         break
#     key = text[s:e-1]
#     s = match.start()
#     e = match.end()
#     
#     pkg_info[key] = s
#     # Move forward in text for the next search
#     pos = e
#     if prev > 0:
#         prev = s - 1