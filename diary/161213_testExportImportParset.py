
from optima.tb.project import Project
import numpy as np

databook = './data/databook-simple-cascade-calibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple-calibration.xlsx')

proj.loadSpreadsheet(databook_path = databook)
# make initial parameter set 
proj.makeParset()
# export to file
proj.exportParset("default")
# load in from file
proj.importParset("default.csv",new_name='def2')

"""
# compare parameter sets, which should be identical
ps_orig = proj.parsets['default']
ps_new  = proj.parsets['def2']
for (i,cid) in enumerate(ps_orig.pars['cascade']):
    print "Original parameter set:"
    print cid
    print "Imported parameter set:"
    print ps_new.pars['cascade'][i]
    print "-----------------------"
"""