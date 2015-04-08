'''
Created on 2 Apr 2015

@author: dave
'''
import unittest
from felt.tests.lib.PyStubTemplate import fn

class test_PyStubTemplate(unittest.TestCase):

	def test(self):
		absval = fn(-1)
		self.assertEqual(absval, 1)


if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()