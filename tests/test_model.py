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
        self.assertEqual(200000, results.outputs['alive']['SAC'][0], 'The Initialised Children Population is not 200000')
        self.assertEqual(200000, results.outputs['alive']['GEN'][0], 'The Initialised Adult Population is not 200000')
        self.assertEqual(200000, results.outputs['alive']['SAC'][-1], 'The Children Population is different from initial population of 200000')
        self.assertEqual(200000, results.outputs['alive']['GEN'][-1], 'The Adult Population is different from initial population of 200000')
        return None
        
#    def test_birth_model(self):
#        """
#        Assumptions:
#            - as for SimpleModel but with births included
#        """
#        results = self.proj.runSim()
#        return None
#    
#    def test_death_model(self):
#        """
#        Assumptions:
#            - as for SimpleModel but with deaths included
#        """
#        results = self.proj.runSim()
#        return None
#    
#    def test_aging_model(self):
#        """
#        Assumptions:
#            - as for SimpleModel but with aging between the two populations
#        """
#        pass
#    
#    def test_simple_model(self):
#        """
#        Assumptions:
#            - as for SimpleModel but with transfers between the two populations
#            - One transfer will be as numbers
#            - The reverse transfer will be as fractions
#        """
#        pass    
    
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
     
#databook = './tests/databooks/databook_model_simple.xlsx'
#cascade =  './tests/cascade_spreadsheet/cascade_model_simple.xlsx'      
databook = 'C:\\Users\\Azfar\\Documents\\tb-ucl\\tests\\databooks\\databook_model_simple.xlsx'
cascade =  'C:\\Users\\Azfar\\Documents\\tb-ucl\\tests\\cascade_spreadsheet\\cascade_model_simple.xlsx'      

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