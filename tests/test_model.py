import unittest
from project import Project


class ModelTest(unittest.TestCase):
    """
    
    """
    cascade = None
    databook = None
    
    def setUp(self):
        self.proj= Project(name = 'unittest', cascade_path = cascade)
        self.proj.loadSpreadsheet(databook_path = databook)
        self.proj.makeParset()

    def tearDown(self):
        self.proj = None
        

class SimpleModel(ModelTest):
    """
    Minimal model: assumes
        - 2 populations
    """
    def test_simple_model(self):
        """
        Assumptions:
            - no aging or transfers between populations
            - no deaths or births
        """
        results = self.proj.runSim()
        #self.assertEqual(results.outputs[])
        return None
        
    def test_birth_model(self):
        """
        Assumptions:
            - as for SimpleModel but with births included
        """
        results = self.proj.runSim()
        return None
    
    def test_death_model(self):
        """
        Assumptions:
            - as for SimpleModel but with deaths included
        """
        results = self.proj.runSim()
        return None
        
                
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
     
#databook = './tests/databooks/simple-cascade.xlsx'
#cascade =  './tests/cascade_spreadsheet/model-simple.xlsx'      
databook = 'C:\\Users\\Azfar\\Documents\\tb-ucl\\tests\\databooks\\databook_simple.xlsx'
cascade =  'C:\\Users\\Azfar\\Documents\\tb-ucl\\tests\\cascade_spreadsheet\\cascade_simple.xlsx'      

if __name__ == '__main__':
    ModelTest.cascade = cascade
    ModelTest.databook = databook
    unittest.main()
    
    """
    suite = unittest.TestSuite()
    suite.addTest(SimpleModel(),)
    result = unittest.TestResult()
    suite.run(result)
    
    """