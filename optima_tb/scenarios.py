import logging
logger = logging.getLogger(__name__)


from uuid import uuid4 as uuid
from copy import deepcopy as dcp

from optima_tb.utils import odict
from optima_tb.parameters import ParameterSet
from optima_tb.databook import getEmptyData


class Scenario(object):

    def __init__(self, name, settings, run_scenario=False, overwrite=True):

        self.name = name
        self.uid = uuid()
        self.run_scenario = run_scenario
        self.overwrite = overwrite
        self.settings = settings

    def getScenarioParset(self, parset):
        return parset

    def getScenarioProgset(self, progset, options):
        return progset, options

    def __repr__(self):
        return '%s "%s" (%s)' % (self.__class__.__name__,self.name,self.uid)

class ParameterScenario(Scenario):

    def __init__(self, name, settings, run_scenario=False, overwrite=True, scenario_values=None):
        """
        Given some data that describes a parameter scenario, creates the corresponding parameterSet 
        which can then be combined with a ParameterSet when running a model.
        
        Params:
            scenario_values:     list of values, organized such to reflect the structure of a linkpars structure in data
                                 data['linkpars'] = {parameter_label : {pop_label : odict o }  where
                                     o = odict with keys:
                                         t : np.array or list with year values
                                         y : np.array or list with corresponding parameter values

                                         
        Example:
            scvalues = odict()
            param = 'birth_transit'
            scvalues[param] = odict()
            scvalues[param]['Pop1'] = odict()
            scvalues[param]['Pop1']['y'] = [3e6, 1e4, 1e4, 2e6]
            scvalues[param]['Pop1']['t'] = [2003.,2004.,2014.,2015.]

            pscenario = ParameterScenario(name="examplePS",scenario_values=scvalues)
    
        """
        super(ParameterScenario, self).__init__(name, settings, run_scenario, overwrite)
        self.scenario_values = scenario_values

    def getScenarioParset(self, parset):
        """
        Get the corresponding parameterSet for this scenario, given an input parameterSet for the default baseline 
        activity. 
        
        The output depends on whether to overwrite (replace) or add values that appear in both the 
        parameterScenario's parameterSet to the baseline parameterSet. 
        """

        # Note - the parset will be overwritten between the first and last year specified in scvalues
        # on a per-parameter+pop basis. Within the scenario, only the data points in scvalues will be used

        import numpy as np
        tvec = np.arange(self.settings.tvec_start, self.settings.tvec_end + self.settings.tvec_dt / 2, self.settings.tvec_dt)

        if self.overwrite:  # update values in parset with those in scenario_parset
            new_parset = dcp(parset)

            for par_label in self.scenario_values.keys():
                par = new_parset.getPar(par_label) # This is the parameter we are updating
                for pop_label,overwrite in self.scenario_values[par_label].items():
                    start_year = min(overwrite['t'])
                    stop_year = max(overwrite['t'])

                    if 'smooth_onset' not in overwrite:
                        onset = self.settings.tvec_dt
                    else:
                        onset = overwrite['smooth_onset']

                    # First, remove values during the overwrite, and add back the onset value
                    # The par values are interpolated first to ensure they exactly match
                    par_dt_vals = par.interpolate(tvec,pop_label)
                    par.t[pop_label] = tvec
                    par.y[pop_label] = par_dt_vals/np.abs(par.y_factor[pop_label]) # Don't forget to divide by the y-factor because this value will get passed through interpolate() later
                    par.removeBetween([start_year-onset, stop_year],pop_label)

                    # Now, insert all of the program overwrites
                    for i in xrange(0,len(overwrite['t'])):
                        par.insertValuePair(overwrite['t'][i], overwrite['y'][i]/par.y_factor[pop_label], pop_label)

            return new_parset

        else:  # add the two together
            # inflate both the two parameter sets first
            parset.inflate(tvec)
            self.scenario_parset.inflate(tvec)
            return parset + self.scenario_parset

class BudgetScenario(Scenario):

    def __init__(self, name, run_scenario=False, overwrite=True, scenario_values=None, pop_labels=None, **kwargs):
        super(BudgetScenario, self).__init__(name, run_scenario, overwrite)
        self.makeScenarioProgset(budget_allocation=scenario_values)
        self.budget_allocation = budget_allocation

    def getScenarioProgset(self, progset, budget_options):
        """
        Get the updated program set and budget allocation for this scenario. 
        This combines the values in the budget allocation with the values for the scenario. 
        
        Note that this assumes that all other budget allocations that are NOT
        specified in budget_options are left as they are. 
        
        Params:
            progset            program set object
            budget_options     budget_options dictionary
        """
        new_budget_options = dcp(budget_options)
        if self.overwrite:
            for prog in self.budget_allocation.keys():
                new_budget_options['init_alloc'][prog] = self.budget_allocation[prog]

        else:  # we add the amount as additional funding
            for prog in self.budget_allocation.keys():

                if new_budget_options['init_alloc'].has_key(prog):
                    new_budget_options['init_alloc'][prog] += self.budget_allocation[prog]
                else:
                    new_budget_options['init_alloc'][prog] = self.budget_allocation[prog]

        return progset, new_budget_options

class CoverageScenario(BudgetScenario):

    def __init__(self, name, run_scenario=False, overwrite=True, scenario_values=None, pop_labels=None, **kwargs):
        super(CoverageScenario, self).__init__(name, run_scenario=run_scenario, overwrite=overwrite, scenario_values=scenario_values, pop_labels=pop_labels, **kwargs)

    def getScenarioProgset(self, progset, options):
        progset, options = super(CoverageScenario, self).getScenarioProgset(progset, options)
        options['alloc_is_coverage'] = True
        return progset, options
