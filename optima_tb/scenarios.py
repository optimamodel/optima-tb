

from uuid import uuid4 as uuid

from optima_tb.utils import odict
from optima_tb.parameters import ParameterSet


class Scenario(object):
    
    def __init__(self,name,run_scenario=False,overwrite=True):
        
        self.name = name
        self.uid  = uuid()
        self.run_scenario = run_scenario
        self.overwrite = overwrite
        self.scenario_parset= None # Placeholder for scenario values
    
    
    def makeScenarioParset(self):
        raise NotImplementedError
    
    def getScenarioParset(self):
        raise NotImplementedError
    
class ParameterScenario(Scenario):
    
    def __init__(self,name,run_scenario=False,overwrite=True,scenario_values=None,pop_labels=None,**kwargs):
        super(ParameterScenario,self).__init__(name,run_scenario,overwrite)
        self.makeScenarioParset(scenario_values,pop_labels=pop_labels)
        
        
    def makeScenarioParset(self,scenario_values,pop_labels):
        """
        
        
        Params:
            scenario_values:     list of values, organized such to reflect the structure of a linkpars structure in data
                                 data['linkpars'] = {parameter_label : {pop_label : odict o }  where
                                     o = odict with keys:
                                         t
                                         y
                                         y_format
                                         y_factor
                                     
            
            
            
            
            
            parameter set:
                                 [ p ] where
                                     p = Parameter with the fields
                                         .label
                                         .y
                                         .t
                                         .y_format    
                                         .y_factor    Assumed to be settings.DO_NOT_SCALE
                                
        """
        data = odict()
        if scenario_values is None:
            scenario_values = odict()
        data['linkpars'] = scenario_values
        # values required when creating a parameter set
        data['characs'] = odict()
        data['transfers'] = odict()
        data['pops'] = odict()
        data['pops']['name_labels'] = pop_labels
        
        ps = ParameterSet(self.name)
        ps.makePars(data)
        self.scenario_parset = ps
    
    def getScenarioParset(self, parset):
        """
        
        """
        if self.overwrite: # update values in parset with those in scenario_parset
            return parset << self.scenario_parset
        else: # add the two together
            return parset + self.scenario_parset
    
    
    
    
    
class BudgetScenario(Scenario):
    
    def __init__(self):
        super(BudgetScenario,self).__init__()
        
        
    def makeScenarioParset(self):
        """
        
        @TODO: implement make scenario for BudgetScenario
        """
        pass
    
class CoverageScenario(Scenario):
    
    def __init__(self):
        super(CoverageScenario,self).__init__()
        
        
    def makeScenarioParset(self):
        """
        
        @TODO: implement make scenario for CoverageScenario
        """
        pass