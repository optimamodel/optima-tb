import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from optima_tb.utils import odict
from optima_tb.project import Project
from optima_tb.plotting import getCategoryColors, plotPopulation


import numpy as np
import pylab

num_pop = 2
plot = True
#plot=False


databook = '../data/databook-simple-cascade-autocalibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = '../data/cascade-simple-calibration.xlsx')
proj.setYear([2000.,2030.],False) 

# setup: define category color list. Note that using colormappings will index by compartment keys, so
# the color maps do not have to be defined in sequence of their appearance within a population's compartment list.
# The only constraint is that all plottable compartments must be assigned a colour. If a compartment is included
# in more than one colormap list, its cmap will be set as last list it appeared on. 
cat_list = odict()
cat_list['Blues'] =    ['sus','vac']
cat_list['Purples'] =  ['ltu','ltt']
cat_list['Reds'] =     ['acu','act']
cat_list['Greens'] =   ['rec']
cat_list['Greys'] =    ['dead']


"""
# 1. Test category color list
print cat_list
col_list = getCategoryColors(cat_list)
print col_list
"""



# 2, load spreadsheet with number format. 
# Test by setting aging --> n
#print proj.settings.linkpar_specs
proj.loadSpreadsheet(databook_path = databook)
proj.makeParset()
results = proj.runSim()
proj.plotResults(results,debug=False,colormappings=cat_list)



# 3. example of plotting subset of compartments
plot_comp_labels = ['acu','act']
plotPopulation(results=results,
               data=proj.data,
               pop_labels=results.pop_labels,
               title=proj.name.title()+" Infectious", 
               colormappings=cat_list,
               plot_observed_label="ac_inf",
               plot_comp_labels=plot_comp_labels,
               plotdict = proj.settings.plot_settings,
               save_fig=True,
               fig_name="Active Infections")

plot_comp_labels = ['ltu','ltt','acu','act','rec']
plotPopulation(results=results,
               data=proj.data,
               pop_labels=results.pop_labels,
               title=proj.name.title()+" Infectious", 
               colormappings=cat_list,
               plot_observed_data=False,
               plot_comp_labels=plot_comp_labels,
               plotdict = proj.settings.plot_settings,
               save_fig=True,
               fig_name="Example Infections")

pylab.show()