'''
testimports.py

Copyright (C) 2013 Liam Deacon
MIT License (see LICENSE file for details)

Test Imports - checks whether key modules and extensions can be imported
'''
from glob import glob1
from os.path import dirname

import importlib
import os
import sys
import unittest


class TestImports(unittest.TestCase):
    """ Test class for imports """

    imports = {'phaseshifts': ['atorb', 'conphas', 'elements', 'factories',
                               'leed', 'libphsh', 'model', 'phsh', 'utils',
                               'wrappers'],
               'phaseshifts.lib': ['libphsh'],
               'phaseshifts.lib.EEASiSSS': ['hf', 'EEASiSSS'],
               'phaseshifts.tests': [os.path.splitext(_module)[0] for _module
                                     in glob1(dirname(__file__), '*.py')
                                     if _module != '__init__.py']
               }

    @staticmethod
    def isimportable(package, _module):
        """
        Determine whether module is importable from given package

        Returns
        -------
        bool : Result is True if import is successful

        """
        try:
            sys.stderr.write('Testing import of %s... ' % _module)
            sys.stderr.flush()
            importlib.import_module(package, _module)
            sys.stderr.write('SUCCESS\n')
            sys.stderr.flush()
            return True
        except ImportError as exception:
            sys.stderr.write('FAILED\n%s' % exception)
            sys.stderr.flush()
            return False

    def testimports(self):
        """
        Function to determine if modules in packages are importable

        Returns
        -------
        FailsIf
            Any module cannot be imported from a given package
        """
        for packages in self.imports:
            sys.stderr.write('Inspecting package: %s\n' % packages)
            successes = [_module for _module in self.imports.get(packages)
                         if self.isimportable('phaseshifts', _module)]
            self.failIf(len(successes) < len(self.imports.get(packages)))
            sys.stderr.write('Failed to import %i out of %i modules\n\n'
                             % (len(self.imports.get(packages)) -
                                len(successes),
                                len(self.imports.get(packages))))
            sys.stderr.flush()

if __name__ == "__main__":
    sys.stderr.write('=' * 72 + '\n')
    sys.stderr.write('TESTING: %s\n' % os.path.basename(__file__))
    sys.stderr.write('=' * 72 + '\n')
    unittest.main()
