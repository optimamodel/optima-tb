import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from utils import odict
from project import Project
from plotting import Plotter
import numpy as np
import pylab

num_pop = 2
plot = True

databook = '../data/databook-simple-cascade-autocalibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = '../data/cascade-simple-calibration.xlsx')
proj.setYear([2000.,2030.],False) 

"""

Purpose: create automatic reporting methods for disease progression and metrics.

Reqs: 
- this has to be for each population (or list of populations), as different populations
can have different user-specified rates.
- this should accept different starting compartments and reporting populations. 
    Each starting compartment should be populated at 100% for some fixed value, and the 
    reporting populations should have the number of people 
- output should be of form "from latent, 45% moved into active after one year, 
    and 54% moved into active after two years"
- plots should also be generated and saved. 

Bonus:
- plots have lines and intercepts indicating yearly progression 

"""

# setup and read-in
proj.loadSpreadsheet(databook_path = databook)
proj.makeParset()

parset = proj.parsets['default']

# define populations and progressions
specified_populations = parset.pop_labels 
specified_progressions = {'sus': ['ltu','acu'],
                          'ltu': ['ltt','acu'],
                          'acu': ['act','rec']}

year_track = np.arange(1.,5.1) # times that we report on 

# turn off aging and other transfers
for transfer_type in parset.transfers.keys():   
    parset.transfers[transfer_type] = odict() 
    
# remove all references to entry points in project.data characteristics so that we can add our own ...
for charac in proj.settings.charac_specs:
    if 'entry_point' in charac.keys():
        del charac['entry_point']

# ... and now we add our own entry points and characteristics for reporting



# reset based on populations and compartments specified
for prog in specified_progressions:
    for (pop, reporters) in specified_populations:
        
        pass

        # set up populations and init values
        
        # run simulations
        
        # get outputs
        
        # save plots 



