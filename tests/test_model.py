from optima_tb.project import Project
from optima_tb.utils import odict
from copy import deepcopy as dcp
import unittest
import numpy as np


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
        dt = self.proj.settings.tvec_dt
        plot_over = [2000,2030]
        self.proj.settings.plot_settings['x_ticks'] = [np.arange(plot_over[0],plot_over[1]+dt,5,dtype=int),np.arange(plot_over[0],plot_over[1]+dt,5,dtype=int)]
        self.proj.settings.plot_settings['xlim'] = (plot_over[0]-0.5,plot_over[1]+0.5)
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
        self.assertEqual(self.proj.data['characs']['alive']['SAC']['y'][0], results.outputs['alive']['SAC'][-1], 'Final Children Population is different from initial population of 200,000')
        print('############################Databook value: %s, Simulation value: %s' %(self.proj.data['characs']['alive']['GEN']['y'][0], results.outputs['alive']['GEN'][-1]))
        self.assertEqual(self.proj.data['characs']['alive']['GEN']['y'][0], results.outputs['alive']['GEN'][-1], 'Final Adult Population is different from initial population of 200,000')
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
        print('############################Databook value: %s, Simulation value: %s' %(self.proj.data['characs']['alive']['GEN']['y'][0], results.outputs['alive']['GEN'][-1]))
        self.assertEqual(self.proj.data['characs']['alive']['GEN']['y'][0], results.outputs['alive']['GEN'][-1], 'The Adult Population is including births')
        self.proj.data['linkpars']['birth_transit']['SAC']['y'][-1] = 0.
        return None
    
    def test_death_model(self):
        """
        Assumptions:
            - as for SimpleModel but with deaths included
        """
        #Test to see whether 10% of untreated population die annually in children
        self.proj.data['linkpars']['mort_u']['GEN']['y'][-1] = 0.1
        results = self.proj.runSim()
        self.assertAlmostEqual(6224, int(results.outputs['death_label']['GEN'][-1]), 4, 'Total number of adult deaths should be approximately 6224 deaths at a rate of 10% per annum without births')
        self.assertEqual(0, int(results.outputs['death_label']['SAC'][-1]), 'Total number of children deaths should be 0 deaths at a rate of 0% per annum without aging')
        self.assertAlmostEqual(193775, int(results.outputs['alive']['GEN'][-1]), 6, 'Adult population should be dying at an annual rate of 10% for untreated cases without births')
        self.assertEqual(200000, int(results.outputs['alive']['SAC'][-1]), 'Children Population is dying even though death rate for adults is 0% without aging')
        self.proj.data['linkpars']['mort_u']['SAC']['y'][-1] = 0.
        
        #Test to see whether 10 people on treatment die annually
        self.proj.data['linkpars']['mort_t']['GEN']['y'][-1] = 10.
        results = self.proj.runSim()        
        self.assertEqual(0, int(results.outputs['death_label']['SAC'][-1]), 'Total number of children deaths should be 0 deaths since number of deaths per year is equal to number 0')
        self.assertEqual(6524, int(results.outputs['death_label']['GEN'][-1]), 'Total number of adult deaths should be 6524 deaths since number of deaths per year is equal to number 10')
        self.assertEqual(200000, int(results.outputs['alive']['SAC'][-1]), 'Children Population is dying even though death rate for adults is 0%')
        self.assertEqual(193475, int(results.outputs['alive']['GEN'][-1]), 'Adult Population is deaths could not be mapped properly')
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
        print('############################Expected value: 200300, Simulation value: %s' %(results.outputs['alive']['GEN'][-1]))
        self.assertEqual(200300, int(results.outputs['alive']['GEN'][-1]), 'Adult population size at end of simulation period should be 200,300 without deaths or aging')
        self.proj.data['transfers']['migration_type_1'] = odict()
        return None
    
    def test_interinfectivity_model(self):
        """
        Assumptions:
            - as for SimpleModel
            - no births
            - no transfers
            - SAC: Everyone is in Susceptible
            - GEN: Significant amount of people are infected
        """
        origproj = dcp(self.proj)
        tomodify = ['lu_prog', 'tmt_a', 'rec_act']
        for par in tomodify:
            for popkey in self.proj.data['linkpars'][par]:
                self.proj.data['linkpars'][par][popkey]['y'][-1] = 0.05
            
        intermediateproj = dcp(self.proj)
        
        #Condition 1: No Interpopulation Infection:
        self.proj.makeParset(name='test_interinfectivity_condition1')
        results1 = self.proj.runSim(parset_name='test_interinfectivity_condition1')
        self.assertEqual(0, int(sum(results1.outputs['lt_inf']['SAC'])+sum(results1.outputs['ac_inf']['SAC'])), 'Children(0-14) population is getting infected with Infection rate = 0 and Interpopulation Infectivity turned off') 
        self.assertGreater(results1.outputs['lt_inf']['GEN'][1], results1.outputs['lt_inf']['GEN'][0], 'Adult population are not getting additional infections with Infection rate > 0 and Interpopulation Infectivity turned off')
        self.assertGreater(results1.outputs['ac_inf']['GEN'][1], results1.outputs['ac_inf']['GEN'][0], 'Adult population are not getting additional infections with Infection rate > 0 and Interpopulation Infectivity turned off')
        self.proj = dcp(intermediateproj)
        
        #Condition 2: Interpopulation Infection:
        self.proj.data['contacts']['from']['GEN']['SAC'] = 1.0
        self.proj.data['contacts']['from']['SAC']['GEN'] = 1.0
        self.proj.data['contacts']['into']['SAC']['GEN'] = 1.0
        self.proj.data['contacts']['into']['GEN']['SAC'] = 1.0
        self.proj.makeParset(name='test_interinfectivity_condition2')
        results2 = self.proj.runSim(parset_name='test_interinfectivity_condition2')
        self.assertGreater(int(results2.outputs['lt_inf']['SAC'][1]), 0, 'Children(0-14) population is not getting infected with Interpopulation Infectivity with equal weighting')
        self.assertGreater(int(results2.outputs['ac_inf']['SAC'][-1]), 0, 'Children(0-14) population is not progressing to active TB in case Interpopulation Infectivity with equal weighting')
        self.proj = dcp(intermediateproj)
        #Condition 3: Interpopulation Infection high rate of interinfectivity in children by adult population
        self.proj.data['contacts']['from']['GEN']['SAC'] = 10.0
        self.proj.data['contacts']['from']['SAC']['GEN'] = 10.0
        self.proj.data['contacts']['into']['SAC']['GEN'] = 10.0
        self.proj.data['contacts']['into']['GEN']['SAC'] = 10.0
        self.proj.makeParset(name='test_interinfectivity_condition3')
        results3 = self.proj.runSim(parset_name='test_interinfectivity_condition3')
        self.assertGreater(int(results3.outputs['lt_inf']['SAC'][1]), int(results2.outputs['lt_inf']['SAC'][1]), 'Increased latent infections in Children(0-14) population is not witnessed when compared to condition 2 with Interpopulation Infectivity with impact children ten weighting')
        self.assertGreater(int(results3.outputs['ac_inf']['SAC'][-1]),int(results2.outputs['ac_inf']['SAC'][-1]), 'Increased active infections in Children(0-14) population is not witnessed when compared to condition 2 with Interpopulation Infectivity with impact children ten times weighting')
        #Reset before exit
        self.proj = dcp(origproj)
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
     
databook = '../tests/databooks/databook_model_simple.xlsx'
cascade =  '../tests/cascade_spreadsheet/cascade_model_simple.xlsx'

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