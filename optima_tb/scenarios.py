

from uuid import uuid4 as uuid


class Scenario:
    
    def __init__(self,name,parset_id):
        
        self.name = name
        self.uid  = uuid()
        self.parset_id = parset_id
    
    
    def makeScenario(self):
        raise NotImplementedError
    
    def runScenario(self, project):
        pass
    
    
class ParameterScenario(Scenario):
    
    def __init__(self,parset):
        super(ParameterScenario,self).__init__()
        self.parameter_set = parset
        
        
    def makeScenario(self):
        """
        
        @TODO: implement make scenario for ParameterScenario
        """
        pass
    
    
class BudgetScenario(Scenario):
    
    def __init__(self):
        super(BudgetScenario,self).__init__()
        
        
    def makeScenario(self):
        """
        
        @TODO: implement make scenario for BudgetScenario
        """
        pass
    
class CoverageScenario(Scenario):
    
    def __init__(self):
        super(CoverageScenario,self).__init__()
        
        
    def makeScenario(self):
        """
        
        @TODO: implement make scenario for CoverageScenario
        """
        pass