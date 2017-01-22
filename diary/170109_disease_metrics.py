import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from optima_tb.utils import odict
from optima_tb.project import Project
import numpy as np
import pylab

from optima_tb.cascade import __addCharacteristic
from optima_tb.databook import __addCharacteristicData

start_year = 2000.
starting_pop = 1e6
dt = 1./4

output_file = "ProgressData.csv"
databook = 'data/databook-simple-cascade-autocalibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple-calibration.xlsx')
proj.setYear([start_year,2030.],False) 
proj.settings.tvec_dt = dt

# setup and read-in
proj.loadSpreadsheet(databook_path = databook)
proj.makeParset() # make a default parset for population labels etc.
data = proj.data

# ================================================
# 2. Define populations and progressions:
# For the moment, we will output for all populations. 
specified_populations = proj.parsets[0].pop_labels 
# To evaluate given transitions between compartments, update 'specified_progressions':
specified_progressions = {'sus': ['ltu','acu'],
                          'ltu': ['ltt','acu','ltu'],
                          'acu': ['act','rec']}
# Finally, specify over which time steps in the future we wish to examine
year_track = [1.,2.,5.] # times that we report on. Note that these can only be multiple of dt
# ================================================




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

def saveResults(filehandle,results,pop,progression_from,progression_to,start_year,year_report):
    """
    Saves the results into a formatted csv file: 
    
    """
    indices = np.in1d(results.t_step-start_year, year_track)
    
    for prog in progression_to:
        filehandle.write("%s --> %s,"%(progression_from,prog))
    
        filehandle.write(",".join(map(str,results.outputs["entry_%s"%prog][pop][indices]/starting_pop)))
        filehandle.write("\n")
        
        
        
        
    
# turn off aging and other transfers
for transfer_type in data['transfers']:
    data['transfers'][transfer_type] = odict() 

# turn off all plotting except for population
for (charac,val) in proj.settings.charac_specs.iteritems():
    val['plot_characteristic'] = 'n'
    val['databook_order'] = -1

# set births to 0: here, we set the rate to 0
for (prog, reporters) in specified_progressions.iteritems():
    for pop in specified_populations:
        data['linkpars']['birth_transit'][pop]['t'] = [start_year]
        data['linkpars']['birth_transit'][pop]['y'] = [0.]
    
# remove all references to entry points in project.data characteristics so that we can add our own ...
for charac,spec in proj.settings.charac_specs.iteritems():
    if 'entry_point' in spec.keys():
        del spec['entry_point']

# ... and now we add our own entry points and characteristics for reporting.
# This is where it gets hacky. We have already loaded settings and data, so we have to add
# to both.
full_compartment_list = list(set([x for v in specified_progressions.itervalues() for x in v]+specified_progressions.keys()))
for prog in full_compartment_list: #specified_progressions.iteritems():
    charac_label='entry_%s'%prog
    __addCharacteristic(proj.settings,charac_label=charac_label,full_name=charac_label,entry_point=prog,includes=[prog])
    proj.settings.charac_specs[charac_label]['plot_characteristic'] = 'n'
    for pop in specified_populations:
        data = __addCharacteristicData(data,charac_label,pop,ts=[start_year],ys=[0.],y_format='number')
        


# Setup simulations, save output and plot:

output_file_handle = open(output_file,'w')

for pop in specified_populations:
    
    output_file_handle.write("%s,"%pop)
    output_file_handle.write(",".join(map(str,year_track)))
    output_file_handle.write("\n")
    
    for (prog, reporters) in specified_progressions.iteritems():

        parset_name = "prog%s_pop%s"%(prog,pop)
        charac_label = 'entry_%s'%prog
        proj.makeParset(parset_name)
        parset = proj.parsets[parset_name]

        # set up populations and init values
        par = parset.pars['characs'][parset.par_ids['characs'][charac_label]]
        par.y[pop][0] = starting_pop
        
        # run simulations
        results = proj.runSim(parset_name=parset_name)

        # get outputs
        saveResults(output_file_handle,results,pop,prog,reporters,start_year,year_track)
        
        # save plots 
        proj.plotResults(results,plot_observed_data=False,savePlot=True,figName=parset_name)
        
        # reset for the next loop
        par.y[pop][0] = 0.
    
    output_file_handle.write("\n")
    
output_file_handle.close()


