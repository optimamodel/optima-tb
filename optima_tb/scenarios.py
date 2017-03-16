

from uuid import uuid4 as uuid
from copy import deepcopy as dcp

from optima_tb.utils import odict
from optima_tb.parameters import ParameterSet
from optima_tb.databook import getEmptyData


class Scenario(object):
    
    def __init__(self,name,settings,run_scenario=False,overwrite=True):
        
        self.name = name
        self.uid  = uuid()
        self.run_scenario = run_scenario
        self.overwrite = overwrite
        self.scenario_parset= None # Placeholder for scenario values
        self.settings = settings
    
    
    def makeScenarioParset(self):
        raise NotImplementedError
    
    def getScenarioParset(self, parset):
        raise NotImplementedError
    
    def makeScenarioProgset(self):
        raise NotImplementedError
    
    def getScenarioProgset(self, progset,options):
        raise NotImplementedError
    
class ParameterScenario(Scenario):
    
    def __init__(self,name,settings,run_scenario=False,overwrite=True,scenario_values=None,pop_labels=None,**kwargs):
        super(ParameterScenario,self).__init__(name,settings,run_scenario,overwrite)
        self.makeScenarioParset(scenario_values,pop_labels=pop_labels)
        
        
    def makeScenarioParset(self,scenario_values,pop_labels):
        """
        Given some data that describes a parameter scenario, creates the corresponding parameterSet 
        which can then be combined with a ParameterSet when running a model.
        
        Params:
            scenario_values:     list of values, organized such to reflect the structure of a linkpars structure in data
                                 data['linkpars'] = {parameter_label : {pop_label : odict o }  where
                                     o = odict with keys:
                                         t : np.array or list with year values
                                         y : np.array or list with corresponding parameter values
                                         y_format : string, describing format of y value. Possible values = 'number', 'fraction' or 'proportion'
                                         y_factor : float, preferably chosen as settings.DO_NOT_SCALE or settings.DEFAULT_YFACTOR
                                         
                                         
        Example:
            from optima_tb.settings import DO_NOT_SCALE
            scvalues = odict()
            param = 'birth_transit'
            scvalues[param] = odict()
            scvalues[param]['Pop1'] = odict()
            scvalues[param]['Pop1']['y'] = [3e6, 1e4, 1e4, 2e6]
            scvalues[param]['Pop1']['t'] = [2003.,2004.,2014.,2015.]
            scvalues[param]['Pop1']['y_format'] = 'number'
            scvalues[param]['Pop1']['y_factor'] = DO_NOT_SCALE
            pops = {'Pop1':'Pop1','Pop2':'Pop2'}
                                     
            pscenario = ParameterScenario(name="examplePS",scenario_values=scvalues,pop_labels=pops)
    
        """
        data = getEmptyData()
        
        if scenario_values is None:
            scenario_values = odict()
        data['linkpars'] = scenario_values
        # values required when creating a parameter set
        data['pops']['label_names'] = pop_labels
        
        ps = ParameterSet(self.name)
        ps.makePars(data)
        self.scenario_parset = ps
    
    def getScenarioParset(self, parset):
        """
        Get the corresponding parameterSet for this scenario, given an input parameterSet for the default baseline 
        activity. 
        
        The output depends on whether to overwrite (replace) or add values that appear in both the 
        parameterScenario's parameterSet to the baseline parameterSet. 
        """
        import numpy as np
        tvec = np.arange(self.settings.tvec_start, self.settings.tvec_end + self.settings.tvec_dt/2, self.settings.tvec_dt)
        
        # inflate both the two parameter sets first
        parset.inflate(tvec)
        self.scenario_parset.inflate(tvec)
             
        if self.overwrite: # update values in parset with those in scenario_parset
            return parset << self.scenario_parset
        else: # add the two together
            return parset + self.scenario_parset
    
    

    def getScenarioProgset(self, progset,options):
        return progset, options
    
    def __repr__(self, *args, **kwargs):
        return "ParameterScenario: \n%s"%self.scenario_parset

    
class BudgetScenario(Scenario):
    
    def __init__(self,name,run_scenario=False,overwrite=True,scenario_values=None,pop_labels=None,**kwargs):
        super(BudgetScenario,self).__init__(name,run_scenario,overwrite)
        self.makeScenarioProgset(budget_allocation=scenario_values)
        
    
    def getScenarioParset(self, parset):
        """

        """
        return parset 

    def makeScenarioProgset(self, budget_allocation):
        """
        Sets up the program set budgetary allocation.
        
        Params:
            budget_allocation     A dict of program label: budget allocation pairs
            
        Example:
            budget_allocation = {'HT-DS': 3.14e6}
            makeScenarioProgset(budget_allocation)
            
        """
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
                
        else: # we add the amount as additional funding
            for prog in self.budget_allocation.keys():
                
                if new_budget_options['init_alloc'].has_key(prog):
                    new_budget_options['init_alloc'][prog] += self.budget_allocation[prog]
                else:
                    new_budget_options['init_alloc'][prog] = self.budget_allocation[prog]
                    
        return progset, new_budget_options
    
    def __repr__(self, *args, **kwargs):
        return "BudgetScenario: \n"+''.join('{}={}\n'.format(key, val) for key, val in sorted(self.budget_allocation.items()))



class CoverageScenario(Scenario):
    
    def __init__(self):
        super(CoverageScenario,self).__init__()
        
        
    def makeScenarioParset(self, parset):
        """
        
        Coverage Scenarios do not make any changes to the ParameterScenario to be used
        """
        return parset
    
    def getScenarioParset(self, parset):
        raise NotImplementedError
    
    
    def getScenarioProgset(self, progset,options):
        raise NotImplementedError
    
    
    
    