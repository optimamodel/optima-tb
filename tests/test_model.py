import unittest
from project import Project
from copy import deepcopy as dcp


class ModelTest(unittest.TestCase):
    """
    
    """
    cascade = None
    databook = None
    
    def setUp(self):
        self.proj= Project(name = 'unittest', cascade_path = cascade)
        self.proj.loadSpreadsheet(databook_path = databook)
        self.proj.settings.tvec_start = 2000.0
        self.proj.settings.tvec_end = 2030.0
        self.proj.settings.tvec_observed_end = 2015.0
        self.proj.settings.tvec_dt = 1.0/4
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
        self.assertEqual(200000, int(results.outputs['alive']['SAC'][-1]), 'Final Children Population is different from initial population of 200000')
        self.assertEqual(200000, int(results.outputs['alive']['GEN'][-1]), 'Final Adult Population is different from initial population of 200000')
        return None
        
    def test_birth_model(self):
        """
        Assumptions:
            - as for SimpleModel but with births included
        """
        self.proj.data['linkpars']['birth_transit']['SAC']['y'][-1] = 100.
        results = self.proj.runSim()
        self.assertEqual(203000, int(results.outputs['alive']['SAC'][-1]), 'The Children Population number of births is not increasing by 100 per year')
        self.assertEqual(200000, int(results.outputs['alive']['GEN'][-1]), 'The Adult Population is including births')
        self.proj.data['linkpars']['birth_transit']['SAC']['y'][-1] = 0.
        return None
    
    def test_death_model(self):
        """
        Assumptions:
            - as for SimpleModel but with deaths included
        """
        #Test to see whether 10% of untreated population die annually in children
        self.proj.data['linkpars']['mort_u']['SAC']['y'][-1] = 0.1
        results = self.proj.runSim()
        self.assertAlmostEqual(190423, int(results.outputs['alive']['SAC'][-1]), 6, 'Children should be dying at an annual rate of 10% using initial value of 10,000(i.e approximately 9,576 deaths)')
        self.assertEqual(200000, int(results.outputs['alive']['GEN'][-1]), 'Adult Population is dying even though death rate for adults is 0%')
        self.proj.data['linkpars']['mort_u']['SAC']['y'][-1] = 0.
        
        #Test to see whether 10 people on treatment die annually
        self.proj.data['linkpars']['mort_t']['GEN']['y'][-1] = 10.
        results = self.proj.runSim()
        self.assertEqual(200000, int(results.outputs['alive']['SAC'][-1]), 'Children Population is dying even though death rate for adults is 0%')
        self.assertEqual(199700, int(results.outputs['alive']['GEN'][-1]), 'Adult Population is dying even though death rate for adults is 0%')
        self.proj.data['linkpars']['mort_t']['GEN']['y'][-1] = 0.
        return None
    
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
#databook = 'C:\\Users\\Azfar\\Documents\\tb-ucl\\tests\\databooks\\databook_model_simple.xlsx'
databook = 'C:\\Users\\Azfar\\Documents\\tb-ucl\\tests\\databooks\\databook_model_simple(counterfactual).xlsx'
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