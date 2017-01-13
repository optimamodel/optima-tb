#%% Imports

import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from utils import tic, toc
from optima_tb.project import Project
import pylab as pl



#%% Test project functionality

tt = tic()

proj = Project(name = 'test', cascade_path = '../data/cascade.xlsx')
proj.makeSpreadsheet(databook_path = '../data/databook-test.xlsx', num_pops = 3, num_migrations = 2)


proj.loadSpreadsheet(databook_path = '../data/databook-test.xlsx')
proj.makeParset()

# run and plot simulations
results = proj.runSim(plot = True, debug = True)

toc(tt, label = 'entire process')

pl.show()