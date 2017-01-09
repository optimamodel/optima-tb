"""
Example to illustrate how to use the unit tests module

Required: 
    The following folders need to exist:
          
Generates:


@date:   06-Jan-2017
@author: hussainazfar

"""
######################################################################
#Will be used to create a suite of unit tests once all tests are ready
######################################################################
from test_model import ModelTest
import unittest

databook = './tests/databooks/databook_model_simple.xlsx'
cascade =  './tests/cascade_spreadsheet/cascade_model_simple.xlsx'

suite = unittest.TestSuite()
suite.addTest(ModelTest())
result = unittest.TestResult()
suite.run(result)