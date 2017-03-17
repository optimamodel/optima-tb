import optima_tb.asd as asd
from optima_tb.model import runModel

import logging
import logging.config
logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

import numpy as np

#def calculateOutcome()

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
        
    if 'outcome_weight' not in options:
        options['outcome_weight'] = {settings.charac_pop_count:'-1'}
            
    print options
    print total_budget

    results = runModel(settings = settings, parset = parset, progset = progset, options = options)
    
    index_start = np.where(results.sim_settings['tvec'] >= options['progs_start'])[0][0]
    print index_start
    print results.sim_settings['tvec'][index_start:]
    charac_label = options['outcome_weight'].keys()[0]
    outcome = 0
    for pop_label in results.outputs[charac_label].keys():
    print outcome
    
    return results
    