#%% DJK System Hack

import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

#%% Actual script

from project import Project
import pylab

proj = Project(name = 'Belarus', cascade_path = './cascade-belarus.xlsx')
proj.makeSpreadsheet(databook_path = './databook-belarus-template.xlsx', num_pops = 6, num_migrations = 2)



proj.loadSpreadsheet(databook_path = './databook-belarus.xlsx')
proj.makeParset(name = 'default')

# run and plot simulations

results = proj.runSim(plot = True)
# calculate a score for how good a fit the model is to the observed data 
proj.calculateFit(results)

pylab.show()