
from project import Project
import numpy as np
import pylab

num_pop = 2
plot = True
#plot=False


databook = './data/databook-simple-cascade-calibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple-calibration.xlsx')
proj.setYear([2000.,2030.],False) 

"""
# 1. make spreadsheet after implementing number format
proj.makeSpreadsheet(databook_path=databook, num_pops = num_pop)

"""
# 2, load spreadsheet with number format. 
# Test by setting aging --> n
#print proj.settings.linkpar_specs
proj.loadSpreadsheet(databook_path = databook)

 
proj.makeParset()
#print proj.data
results1 = proj.runSim(plot=plot)
pylab.show()
