#%% DJK System Hack

import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

"""
Run model simulation and plot output:

Loads and runs model data

"""

from optima_tb.project import Project
import pylab

proj = Project(name = 'Belarus', cascade_path = './cascade-belarus.xlsx', validation_level = 'avert')

#set the year range we simulate over as starting in 1995:
proj.setYear([1995,2030],False)


proj.loadSpreadsheet(databook_path = './databook-belarus-template.xlsx')
proj.makeParset(name = 'default')

# run and plot simulations

results = proj.runSim(plot = True)
# calculate a score for how good a fit the model is to the observed data 
proj.calculateFit(results)

pylab.show()