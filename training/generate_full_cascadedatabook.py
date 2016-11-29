from project import Project
import pylab 

"""
Example to illustrate how to run and modify for simple cascade structure.

Required: full-cascade.xlsx
          databook-full-cascade.xlsx
          
Generates: sample-output.xlsx, which will be the same format as databook-simple-cascade.xlsx

Note: this script should be run from the home folder for this project i.e. ~/git/tb-ucl/


@date:   14-Nov-2016
@author: sjarvis

"""

"""
Populations are:
    - 0-4
    - 5-14
    - 15-64
    - 65+
    - Prisoners
    - 15-64 HIV
    - 65+ HIV
"""
num_pop = 7 
#this is a test

databook = './training/Belarus-databook-cascade_161122.xlsx'
cascade = './training/cascade_161122.xlsx'


proj= Project(name = 'test-full-cascade', cascade_path = cascade)
# Set the year range we wish to enter data points for: from 2000 to 2016 inclusive
proj.setYear([2000.,2016.])


# ----------------------------------------------------------------
# 1. make spreadsheet after implementing number format
proj.makeSpreadsheet(databook_path=databook, num_pops = num_pop)

