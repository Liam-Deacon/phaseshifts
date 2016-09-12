#!/usr/bin/env python
# encoding: utf-8
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.m.deacon@gmail.com                                           #                                                #
#                                                                            #
# Copyright: Copyright (C) 2015 Liam Deacon                                  #
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
"""
**utils.py**

Provides utility functions and classes for performing common low-level tasks.
"""
from __future__ import absolute_import, division, with_statement
from __future__ import print_function, unicode_literals

from shutil import copy
import os
import sys


def expand_filepath(path):
    """
    Expands `path` for environment and user variables

    Examples
    --------
    >>> from phaseshifts.utils import expand_filepath
    >>> expand_filepath('~/')
     'C:\\Users\\Liam'
    >>> expand_filepath(r'parent\dir/file.ext')
     'parent\\dir\\file.ext'

    """
    return os.path.normpath(
        os.path.expanduser(
            os.path.expandvars(
                os.path.expanduser(path))))


def fix_path(file_path, fill_char=''):
    """
    Fixes escaped characters in `file_path`. Implicitly calls
    :py:meth:`expand_filepath`

    .. note:: Offending path characters will be substituted with `fill_char`.

    Examples
    --------
    >>> from phaseshifts.utils import fix_path
    >>> fix_path(r"C:\test\a_file.txt")
     u'C:\\test\\a_file.txt'
    """
    if sys.platform.lower().startswith('win'):

        file_path = expand_filepath(file_path)
        fix_list = {'\a': '\\a', '\b': '\\b',
                    '\f': '\\f', '\n': '\\n',
                    '\r': '\\r', '\t': '\\t',
                    '\v': '\\v', '\\\\': '\\'}
        for fix in fix_list:
            file_path = file_path.replace(fix, fix_list[fix])

        for fix in fix_list:
            file_path = file_path.replace(fix, fix_list[fix])

    fill_char = '' if not isinstance(fill_char, str) else fill_char

    return "".join(x if x.isalnum() or x in list(':\\/-_.')
                   else fill_char for x in file_path)


def stringify(arg):
    """
    Returns string of `arg` or fancy string of items if `arg` is a
    :py:obj:`dict`, :py:obj:`list` or :py:obj:`tuple`.

    Raises
    ------
    TypeError
        If `arg` is unicode and cannot be coerced into a formatted string.

    Examples
    --------
    >>> from phaseshifts.utils import stringify
    >>> stringify(None)
     u'None'
    >>> stringify([1, 'two', None])
     u"'1', 'two', or 'None'"
    >>> stringify({3: 'three', 'four': 4, False: None})
     u"'3', 'four', or 'False'"

    """
    try:
        if (isinstance(arg, list) or
                isinstance(arg, tuple) or
                isinstance(arg, dict)):
            arr = ("'" + "', '".join(["{}".format(item)
                                      for item in arg]) + "'")
            var = arr.split()
            var.insert(len(arg) - 1, 'or')
            return " ".join(var)
        else:
            return "{}".format(arg)
    except TypeError:
        # !TODO: error handling can go here
        raise TypeError


class FileUtils(object):
    """
    Class for performing phase shift related file operations
    """

    def __init__(self, params):
        """
        Constructor
        """
        pass

    @staticmethod
    def copy_files(files, dst, verbose=False):
        """copy list of files into destination directory"""
        # check if using native Windows Python with cygwin
        env = ''
        if str(sys.platform).startswith('win') and dst.startswith('/cygdrive'):
            if os.environ['CLEED_PHASE'] == dst:
                env = 'CLEED_PHASE='
            dst = '"%s"' % (dst.split('/')[2] + ':' +
                            os.path.sep.join(dst.split('/')[3:]))

        # do check and create directory if needed
        if os.path.isfile(dst):
            dst = os.path.dirname(dst)
        if not os.path.exists(dst):
            try:
                os.makedirs(dst)
            except WindowsError:
                pass

        # copy each phase shift file to directory
        if verbose:
            print("\nCopying files to %s'%s'" % (env, dst))
        for filename in files:
            try:
                copy(filename, dst)
                if verbose:
                    print(os.path.basename(filename))
            except IOError:
                sys.stderr.write("Cannot copy file '%s'\n" % filename)
                sys.stderr.flush()
