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



proj.loadSpreadsheet(databook_path = './databook-belarus-template.xlsx')
proj.makeParset(name = 'default')

# run and plot simulations
r1, o1, s1, results1 = proj.runSim(parset_name = 'default', plot = True)
#proj.runAutofitCalibration(new_parset_name='autocalibrated')
#r2, o2, s2, results2 = proj.runSim(parset_name = 'autocalibrated', plot = True)

pylab.show()