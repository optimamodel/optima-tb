#%% Imports

from optima_tb.utils import odict, OptimaException
import optima_tb.settings as project_settings

import logging
logger = logging.getLogger(__name__)

import numpy as np
from copy import deepcopy as dcp

#%% Functions to generate default structures

def defaultOptimOptions(settings, progset = None, progs_start = None):
    
    options = dict()
    
    if progs_start is None: options['progs_start'] = 2015.0
    else: options['progs_start'] = progs_start
    options['progs_end'] = np.inf
    options['init_alloc'] = odict()
    options['constraints'] = {'limits':odict(), 'max_yearly_change':odict(), 'impacts':odict()}
    
    if not progset is None:
        for prog in progset.progs:
            options['init_alloc'][prog.label] = prog.getDefaultBudget(year=progs_start)
            options['constraints']['limits'][prog.label] = {'vals':[0.0,np.inf],'rel':True}
            if prog.func_specs['type'] == 'cost_only':
                options['constraints']['limits'][prog.label]['vals'] = [1.0,1.0]
            else:
                # All programs that are not fixed-cost can have a default ramp constraint.
                # This should be fine for fixed-cost programs too, but is redundant and does not need to be explicit.
                options['constraints']['max_yearly_change'][prog.label] = {'val':np.inf, 'rel':True}
        for impact in progset.impacts.keys():
            options['constraints']['impacts'][impact] = {'vals':[0.0,np.inf]}

    options['constraints']['total'] = sum(options['init_alloc'].values())
    options['objectives'] = {settings.charac_pop_count : {'weight':-1,'year':2030.0}}
    options['saturate_with_default_budgets'] = True     # Set True so that optimization redistributes funds across entire progset.
    
    return options