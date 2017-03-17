import optima_tb.asd as asd
from optima_tb.model import runModel

import logging
import logging.config
logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()



def optimizeFunc(settings, parset, progset, options = None):
    
    if options is None:
        logger.warn("An options dictionary was not supplied for optimisation. A default one will be constructed.")
        options = dict()
        
    if 'progs_start' not in options:
        options['progs_start'] = 2015.0
        logger.warn("Programs will start overwriting calibrated parameter values in the year: %f" % options['progs_start'])

    if 'init_alloc' not in options:
        options['init_alloc'] = {}
        
    total_budget = 0
    for prog in progset.progs:
        if prog.label not in options['init_alloc']:
            options['init_alloc'][prog.label] = prog.getDefaultBudget()
        total_budget += options['init_alloc'][prog.label]
            
    print options
    print total_budget

    results = runModel(settings = settings, parset = parset, progset = progset, options = options)
    
    return results 
        
#options['init_alloc'] = {'HT-DS':1e6,'SAT-DS':1e6,'HT-MDR':1e6}     # Which programs will have user-specified initial budgets?
      
        