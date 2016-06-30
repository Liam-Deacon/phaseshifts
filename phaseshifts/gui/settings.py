#!/usr/bin/env python
# encoding: utf-8
'''
Created on 4 Feb 2014

@author: Liam Deacon

@contact: liam.m.deacon@gmail.com

@copyright: Copyright 2014 Liam Deacon

@license: MIT License 

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to deal 
in the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to 
do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
'''

from os.path import join, expandvars
import sys


class Settings(object):
    '''
    Class to hold settings
    '''
    if sys.platform == 'win32':  # windows paths and environment variables
        DEFAULT_PATH = join(expandvars('%HOMEDRIVE%%HOMEPATH%'), 
                            '.phaseshifts')
    else:  # assume unix environment variables and paths
        DEFAULT_PATH = join(expandvars('${HOME}'), '.phaseshifts') 

    def __init__(self, grid=1000, mixing=0.5, tolerence=0.0005, ech=100, 
                 relic=0.0, atorb_exchange=0.0, relativistic_atorb=True, 
                 write_files=False, delete_files_after_use=False, 
                 filepath=DEFAULT_PATH, mufpot_exchange=0.0, mufpot_nh=10,  
                 phsh_method='Relativistic', phsh_format='CLEED', 
                 logfile=join(DEFAULT_PATH, 'log.txt'), logging=True,
                 **kwargs):
        '''
        Constructor
        '''
        self.members_dict = {
                        'grid': 'int',
                        'mixing': 'float',
                        'tolerence': 'float',
                        'ech': 'float',
                        'relic': 'float',
                        'atorb_exchange': 'float',
                        'relativistic_atorb': 'bool',
                        'write_files': 'bool',
                        'delete_files_after_use': 'bool',
                        'filepath': 'str',
                        'mufpot_exchange': 'float',
                        'mufpot_nh': 'int',
                        'phsh_method': 'str',
                        'phsh_format': 'str',
                        'logfile': 'str',
                        'logging': 'bool'
                        }
        
        for param in self.members_dict:
            exec('self.{0} = {0}'.format(param)) 
        
        # add extra arguments
        for kwarg in kwargs:
            self.members_dict.update({kwarg: str(type(kwargs.get(kwarg))
                        ).replace("<type ", '').replace("'", '').strip('>')})

        self.__dict__.update(kwargs)  # update class dictionary
        
    def update(self, **kwargs):
        '''
        Lazy update method
        '''
        for kwarg in kwargs:
            s = eval(self.members_dict.get(kwarg))
            if isinstance(kwargs.get(kwarg), s):
                self.__dict__.update({kwarg: kwargs.get(kwarg)})
