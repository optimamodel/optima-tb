#%% Imports

from optima_tb.utils import odict, OptimaException
import optima_tb.settings as project_settings

import logging
logger = logging.getLogger(__name__)

import numpy as np

#%% Functions to generate default structures

def defaultOptimOptions(settings, progset = None):
    
    options = dict()
    
    options['progs_start'] = 2015.0
    options['init_alloc'] = odict()
    options['constraints'] = {'limits':odict()}
    
    if not progset is None:
        for prog in progset.progs:
            options['init_alloc'][prog.label] = prog.getDefaultBudget()
            options['constraints']['limits'][prog.label] = {'vals':[0.0,np.inf],'rel':True}
            if prog.func_specs['type'] == 'cost_only':
                options['constraints']['limits'][prog.label]['vals'] = [1.0,1.0]

    options['constraints']['total'] = sum(options['init_alloc'].values())
    options['objectives'] = {settings.charac_pop_count : {'weight':-1,'year':2030.0}}
    
    return options