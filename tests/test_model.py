import unittest
from project import Project


class ModelTest(unittest.TestCase):
    """
    
    
    
    """
    def setUp(self,cascade,databook):
        self.proj= Project(name = 'unittest', cascade_path = cascade)
        self.proj.loadSpreadsheet(databook_path = databook)
        self.proj.makeParset()

    def tearDown(self):
        self.proj = None
        

class SimpleModel(ModelTest):
    """
    Minimal model: assumes
        - 2 populations
        - no aging or transfers between populations
        - no deaths or births
    """
    def test_simple_model(self):
        results = self.proj.runSim()


class DeathModel(ModelTest):
    """
    Minimal model: assumes
        - as for SimpleModel but with deaths included
    """
    def test_simple_model(self):
        results = self.proj.runSim()

      
class BirthModel(ModelTest):
    """
    Minimal model: assumes
        - as for SimpleModel but with births included
    """
    def test_simple_model(self):
        results = self.proj.runSim()
        
                
class AgingModel(ModelTest):
    """
    Minimal model: assumes
        - as for SimpleModel but with aging between the two populations
    """
    
    def test_simple_model(self):
        pass


class TransfersModel(ModelTest):
    """
    Minimal model: assumes
        - as for SimpleModel but with transfers between the two populations
    """
    
    def test_simple_model(self):
        pass
    
    
class FullModel(ModelTest):
    """
    This is the model with the works ... 
        - 4 populations
        - aging between 3 of the populations 
        - transfers between 2 of the populations
        - deaths 
        - births
    """
    def test_full_model(self):
        pass     
     

class EvilModels(ModelTest):
    """
    Models that we know should fail:
        - models that are all junctions
        - ... 
    
    @TODO: implement
    """
    pass    
     
          
        
if __name__ == '__main__':
    unittest.main()
    
    """
    suite = unittest.TestSuite()
    suite.addTest(SimpleModel(),)
    result = unittest.TestResult()
    suite.run(result)
    
    """