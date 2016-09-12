'''
Created on 10 Jul 2016

@author: Liam Deacon
'''
import unittest


class TestImportDialog(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestImportDialog, cls).setUpClass()
        from qtsix.Qt import QApplication
        cls.app = QApplication()

    @classmethod
    def tearDownClass(cls):
        super(TestImportDialog, cls).tearDownClass()
        cls.app.exec_()

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testName(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
