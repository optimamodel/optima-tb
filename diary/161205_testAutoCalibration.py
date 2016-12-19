import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from project import Project
import numpy as np
import pylab 

"""
Aim: test autocalibration

"""

num_pop = 2
plot = True
#plot=False


databook = './data/databook-simple-cascade-autocalibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = './data/cascade-simple-calibration.xlsx')
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

#r1, o1,s1,results1 = proj.runSim(plot=plot)

# Autofit over *all* characteristics:
proj.runAutofitCalibration(new_parset_name='bob')
# Autofit over a subset of characteristics:
proj.runAutofitCalibration(new_parset_name='eve',target_characs=['vaccin','lt_treat','at_treat','lt_inf','ac_inf'])
# Run using the second autofitted parameters
r2, o2,s2,results2 = proj.runSim(parset_name='eve',plot=plot)

pylab.show()
