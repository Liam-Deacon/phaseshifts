#!/usr/bin/env python
# encoding: utf-8
'''
Created on 26 Feb 2014

@author: Liam Deacon

@contact: liam.deacon@diamond.ac.uk

@copyright: 2014 Liam Deacon

@license: MIT License

'''

import sys
import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__version__ = 0.1
__date__ = '2014-02-23'
__updated__ = '2014-02-23'
__contact__ = 'liam.deacon@diamond.ac.uk'

class Debye_Waller(object):
    '''
    Debye_Waller is a class for calculating Debye-Waller factors
    '''
    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def debye_waller_factor(self):
        '''
        Calculate the Debye-Waller factor
        '''
        pass
    

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
        
    def __str__(self):
        return self.msg
    
    def __unicode__(self):
        return self.msg


def main(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, 
                                                     program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

      Created by Liam Deacon on %s.
      Copyright 2013-2014 Liam Deacon. All rights reserved.

      Licensed under the MIT license (see LICENSE file for details)

      Please send your feedback, including bugs notifications
      and fixes, to: %s

    usage:-
    ''' % (program_shortdesc, str(__date__), __contact__)

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, 
                                formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('-b', '--bulk', dest='bulk', metavar='<bulk_file>', 
                            help="path to MTZ bulk input file")
        
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    
    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2