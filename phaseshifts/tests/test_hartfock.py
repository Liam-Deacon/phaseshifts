'''
Created on 27 Feb 2015

@author: kss07698
'''
import unittest
from phaseshifts.lib.hartfock import HartFock as hartfock 

class Test(unittest.TestCase):


    def testHartFockInit(self):
        ''' Test HartFock class initialization '''
        hf = hartfock()
        assert(isinstance(hf, hartfock))

    def testHartFockReadInputFile(self):
        ''' Test HartFock input '''
        hf = hartfock()
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()