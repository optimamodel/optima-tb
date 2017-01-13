from optima_tb.project import Project

import os

import pylab 


"""
Example to illustrate how to run and modify for simple cascade structure.

Required: simple-cascade.xlsx
          databook-simple-cascade.xlsx
          
Generates: sample-output.xlsx, which will be the same format as databook-simple-cascade.xlsx

Note: this script should be run from the home folder for this project i.e. ~/git/tb-ucl/


@date:   14-Nov-2016
@author: sjarvis

"""

num_pop = 2

databook = os.path.abspath('tests/databooks/databook_model_simple(counterfactual).xlsx')
cascade = os.path.abspath('tests/cascade_spreadsheet/cascade_model_simple.xlsx')


proj= Project(name = 'test-simple-birth', cascade_path = cascade)
# Set the year range we wish to enter data points for: from 2000 to 2016 inclusive
proj.setYear([2000.,2016.])
proj.loadSpreadsheet(databook_path = databook)
proj.makeParset()
r1 = proj.runSim(plot = False)

"""
# ----------------------------------------------------------------
# 1. make spreadsheet after implementing number format
proj.makeSpreadsheet(databook_path="sample-output.xlsx", num_pops = num_pop)

# Note that we'll output the data sheet to another name to prevent accidentally overwriting it.
# Normally, we would probably set the databook path to be:
#proj.makeSpreadsheet(databook_path=databook num_pops = num_pop)

# ----------------------------------------------------------------
# ----------------------------------------------------------------
# 2, load spreadsheet with number format and run simulation

proj.loadSpreadsheet(databook_path = databook)
 
proj.makeParset()
r1, o1 = proj.runSim(plot=True)
pylab.show()
# ----------------------------------------------------------------
"""

