'''
testimports.py

Copyright (C) 2013 Liam Deacon
MIT License (see LICENSE file for details)

Test Imports - checks whether key modules and extensions can be imported
'''
import sys
import os
import unittest
import glob


def isimportable(package, module):
    """
    Determine whether module is importable from given package
    
    Returns
    -------
    bool : Result is True if import is successful
    
    """
    try:
        sys.stderr.write('Testing import of %s... ' % module)
        sys.stderr.flush()
        exec('from %s import %s' % (package, module))
        sys.stderr.write('SUCCESS\n')
        sys.stderr.flush()
        return True
    except ImportError as e:
        sys.stderr.write('FAILED\n')
        sys.stderr.flush()
        print(e)
        return False

imports = {'phaseshifts': ['atorb', 'conphas', 'elements', 'factories', 
                           'leed', 'libphsh', 'model', 'phsh', 'utils',
                           'wrappers'],
           'phaseshifts.lib': ['libphsh'],
           'phaseshifts.lib.EEASiSSS': ['hf', 'EEASiSSS'], 
           'phaseshifts.tests': [os.path.splitext(module)[0] for module in 
                                 glob.glob1(os.path.dirname(__file__), '*.py')
                                 if module != '__init__.py']
           }


class TestImports(unittest.TestCase):
    """ Test class for imports """
    def testimports(self):
        """ 
        Function to determine if modules in packages are importable 
        
        Returns
        -------
        FailsIf
            Any module cannot be imported from a given package 
        """
        for packages in imports:
            sys.stderr.write('Inspecting package: %s\n' % packages)
            successes = [module for module in imports.get(packages) 
                         if isimportable('phaseshifts', module)]
            self.failIf(len(successes) < len(imports.get(packages)))
            sys.stderr.write('Failed to import %i out of %i modules\n\n' 
                             % (len(imports.get(packages)) - len(successes), 
                                len(imports.get(packages))))
            sys.stderr.flush()

if __name__ == "__main__":
    sys.stderr.write('=' * 72 + '\n')
    sys.stderr.write('TESTING: %s\n' % os.path.basename(__file__))
    sys.stderr.write('=' * 72 + '\n')
    unittest.main()
    