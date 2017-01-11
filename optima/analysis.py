import logging
logger = logging.getLogger(__name__)

import numpy as np

from utils import OptimaException, odict
from cascade import __addCharacteristic
from databook import __addCharacteristicData



def __saveProgressionResults(filehandle,results,pop,progression_from,progression_to,start_year,year_report,starting_pop):
    """
    Saves the results into a formatted csv file: 
    
    """
    indices = np.in1d(results.t_step-start_year, year_report)
    
    for prog in progression_to:
        filehandle.write("%s --> %s,"%(progression_from,prog))
    
        filehandle.write(",".join(map(str,results.outputs["entry_%s"%prog][pop][indices]/starting_pop)))
        filehandle.write("\n")
        
        
def evaluateDiseaseProgression(proj, specified_progressions, specified_populations=None,starting_pop = 1e6,
                               year_track=[1.,2.,3.,4.,5.], birth_transit='birth_transit',output_file="ProgressionData.csv"):
    """
    
    
    specified_progressions = {'sus': ['ltu','acu'],
                              'ltu': ['ltt','acu','ltu'],
                              'acu': ['act','rec']}
    
    
    """   
    if specified_populations is None:
        specified_populations = proj.parsets[0].pop_labels 
        logger.info("Evaluating disease progression for all populations")
    
    data = proj.data
    dt = proj.settings.tvec_dt
    start_year = proj.settings.tvec_start
    
        
    # turn off aging and other transfers
    for transfer_type in data['transfers']:
        data['transfers'][transfer_type] = odict() 
    
    # turn off all plotting except for population
    for (charac,val) in proj.settings.charac_specs.iteritems():
        val['databook_order'] = -1
    
    # set births to 0: here, we set the rate to 0
    for (prog, reporters) in specified_progressions.iteritems():
        for pop in specified_populations:
            data['linkpars'][birth_transit][pop]['t'] = [start_year]
            data['linkpars'][birth_transit][pop]['y'] = [0.]
        
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
        for pop in specified_populations:
            data = __addCharacteristicData(data,charac_label,pop,[start_year],[0.],'number')
            
    
    
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
            __saveProgressionResults(output_file_handle,results,pop,prog,reporters,start_year,year_track,starting_pop)
            
            # save plots 
            proj.plotResults(results,plotObservedData=False,savePlot=True,figName=parset_name)
            
            # reset for the next loop
            par.y[pop][0] = 0.
        
        output_file_handle.write("\n"*2)
        
    output_file_handle.close()
