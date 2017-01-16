from optima_tb.project import Project

import pylab 

"""
Aim: test and validation that implementation of births

"""

num_pop = 2
plot = False

    
databook = './data/databook-simple-cascade.xlsx'
proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple.xlsx')
proj.setYear([2000.,2015.],False) # simulate to 2031
"""
# 1. make spreadsheet after implementing number format
proj.makeSpreadsheet(databook_path=databook, num_pops = num_pop)

"""
# 2, load spreadsheet with number format. 
# Test by setting aging --> n
proj.loadSpreadsheet(databook_path = databook)
 
 
proj.makeParset()
print proj.settings.charac_specs
#print proj.data
plot=True
r1 = proj.runSim(plot=plot)
pylab.show()
