import unittest
from project import Project
from copy import deepcopy as dcp
import numpy as np
from utils import odict


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
        self.assertEqual(200000, int(results.outputs['alive']['SAC'][-1]), 'Final Children Population is different from initial population of 200,000')
        self.assertEqual(200000, int(results.outputs['alive']['GEN'][-1]), 'Final Adult Population is different from initial population of 200,000')
        return None
        
    def test_birth_model(self):
        """
        Assumptions:
            - as for SimpleModel but with births included
        """
        self.proj.data['linkpars']['birth_transit']['SAC']['y'][-1] = 100.
        results = self.proj.runSim()
        self.assertEqual(-3000, int(results.outputs['birth_label']['SAC'][-1]), 'Cummulative number of child births for children at end of simulation period is incorrect')
        self.assertEqual(0, int(results.outputs['birth_label']['GEN'][-1]), 'Cummulative number of adult births at end of simulation period is non-zero')
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
        self.assertAlmostEqual(9576, int(results.outputs['death_label']['SAC'][-1]), 6, 'Total number of children deaths should be approximately 9,576 deaths at a rate of 10% per annum without births')
        self.assertEqual(0, int(results.outputs['death_label']['GEN'][-1]), 'Total number of adult deaths should be 0 deaths at a rate of 0% per annum without aging')
        self.assertAlmostEqual(190423, int(results.outputs['alive']['SAC'][-1]), 6, 'Children should be dying at an annual rate of 10% using initial value of 10,000(i.e approximately 9,576 deaths) without births')
        self.assertEqual(200000, int(results.outputs['alive']['GEN'][-1]), 'Adult Population is dying even though death rate for adults is 0% without aging')
        self.proj.data['linkpars']['mort_u']['SAC']['y'][-1] = 0.
        
        #Test to see whether 10 people on treatment die annually
        self.proj.data['linkpars']['mort_t']['GEN']['y'][-1] = 10.
        results = self.proj.runSim()
        self.assertEqual(0, int(results.outputs['death_label']['SAC'][-1]), 'Total number of children deaths should be 0 deaths since number of deaths per year is equal to number 0')
        self.assertEqual(300, int(results.outputs['death_label']['GEN'][-1]), 'Total number of adult deaths should be 300 deaths since number o deaths per year is equal to number 10')
        self.assertEqual(200000, int(results.outputs['alive']['SAC'][-1]), 'Children Population is dying even though death rate for adults is 0%')
        self.assertEqual(199700, int(results.outputs['alive']['GEN'][-1]), 'Adult Population is dying even though death rate for adults is 0%')
        self.proj.data['linkpars']['mort_t']['GEN']['y'][-1] = 0.
        return None
    
    def test_aging_model(self):
        """
        Assumptions:
            - as for SimpleModel but with aging between the two populations
        """
        self.proj.data['transfers']['aging']['SAC'] = odict()
        self.proj.data['transfers']['aging']['SAC']['GEN'] = odict()
        self.proj.data['transfers']['aging']['SAC']['GEN']['t'] = np.array([2000.])
        self.proj.data['transfers']['aging']['SAC']['GEN']['y'] = np.array([1./15.])
        self.proj.data['transfers']['aging']['SAC']['GEN']['y_format'] = 'fraction'
        self.proj.data['transfers']['aging']['SAC']['GEN']['y_factor'] = 1.0
        self.proj.makeParset(name='test_aging')
        results = self.proj.runSim(parset_name='test_aging')
        self.assertAlmostEqual(25242, int(results.outputs['alive']['SAC'][-1]), 5, 'Children(0-14) population size at end of simulation period should be approximately 25,242 without deaths')
        self.assertAlmostEqual(374757, int(results.outputs['alive']['GEN'][-1]), 6, 'Adult population size at end of simulation period should be approximately 374,757 without deaths')
        self.proj.data['transfers']['aging'] = odict()
        return None
    
    def test_transfer_model(self):
        """
        Assumptions:
            - as for SimpleModel but with transfers between the two populations
            - One transfer will be as numbers
            - The reverse transfer will be as fractions
        """
        #Transfer from GEN -> SAC (i.e. migration type 2) at a fraction of 10% per annum
        self.proj.data['transfers']['migration_type_2']['GEN'] = odict()
        self.proj.data['transfers']['migration_type_2']['GEN']['SAC'] = odict()
        self.proj.data['transfers']['migration_type_2']['GEN']['SAC']['t'] = np.array([2000.])
        self.proj.data['transfers']['migration_type_2']['GEN']['SAC']['y'] = np.array([0.1])
        self.proj.data['transfers']['migration_type_2']['GEN']['SAC']['y_format'] = 'fraction'
        self.proj.data['transfers']['migration_type_2']['GEN']['SAC']['y_factor'] = 1.0
        self.proj.makeParset(name='test_transferfraction')
        results = self.proj.runSim(parset_name='test_transferfraction')
        self.assertAlmostEqual(391521, int(results.outputs['alive']['SAC'][-1]), 6, 'Children(0-14) population size at end of simulation period should be approximately 391,521 without deaths, births or aging')
        self.assertAlmostEqual(8478, int(results.outputs['alive']['GEN'][-1]), 4, 'Adult population size at end of simulation period should be approximately 8,478 without deaths or aging')
        self.proj.data['transfers']['migration_type_2'] = odict()
        
        #Transfer from SAC -> GEN (i.e. migration type 1) at an annual figure of 10
        self.proj.data['transfers']['migration_type_1']['SAC'] = odict()
        self.proj.data['transfers']['migration_type_1']['SAC']['GEN'] = odict()
        self.proj.data['transfers']['migration_type_1']['SAC']['GEN']['t'] = np.array([2000.])
        self.proj.data['transfers']['migration_type_1']['SAC']['GEN']['y'] = np.array([10.])
        self.proj.data['transfers']['migration_type_1']['SAC']['GEN']['y_format'] = 'number'
        self.proj.data['transfers']['migration_type_1']['SAC']['GEN']['y_factor'] = 1.0
        self.proj.makeParset(name='test_transfernumber')
        results = self.proj.runSim(parset_name='test_transfernumber')
        self.assertEqual(199700, int(results.outputs['alive']['SAC'][-1]), 'Children(0-14) population size at end of simulation period should be 199,700 without deaths, births or aging')
        self.assertEqual(200300, int(results.outputs['alive']['GEN'][-1]), 'Adult population size at end of simulation period should be 200,300 without deaths or aging')
        self.proj.data['transfers']['migration_type_1'] = odict()
        
        
        
        return None
    
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
     
databook = './tests/databooks/databook_model_simple.xlsx'
cascade =  './tests/cascade_spreadsheet/cascade_model_simple.xlsx'

if __name__ == '__main__':
    ModelTest.cascade = cascade
    ModelTest.databook = databook
    unittest.main()
"""
databook = './tests/databooks/databook_model_simple.xlsx'
cascade =  './tests/cascade_spreadsheet/cascade_model_simple.xlsx'

suite = unittest.TestSuite()
suite.addTest(ModelTest())
result = unittest.TestResult()
suite.run(result)
"""