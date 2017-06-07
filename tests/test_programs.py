from optima_tb.project import Project
import numpy as np
import pylab
import unittest
import os

#Spreadsheets to use
databook = os.path.abspath('../tests/databooks/databook_model_full.xlsx')
cascade =  os.path.abspath('../tests/cascade_spreadsheet/cascade_model_full.xlsx')

#Setup Unit Test class
class CreateProgsetTest(unittest.TestCase):
    def setUp(self):
        self.proj= Project(name = 'unittest', cascade_path = cascade)
        self.proj.loadSpreadsheet(databook_path=databook)
        self.proj.setYear([2000,2015], False)
        self.proj.makeParset()

    def tearDown(self):
        self.proj = None
        pylab.close('all')
        
class ProgramsTest(unittest.TestCase):
    def setUp(self):
        self.proj= Project(name = 'unittest', cascade_path = cascade)
        self.proj.loadSpreadsheet(databook_path=databook)
        self.proj.setYear([2000,2015], False)
        self.proj.makeParset()
        self.proj.makeProgset()
        self.prog_short_names = ['ds-tb', 'mdr-tb', 'vac-tb', 'cure-tb', 'fixed']
        self.prog_full_names  = ['DS Treatment Program', 'MDR Treatment Program', 
                            'Susceptible Vaccination Program', 'Latent Cure Program',
                            'Fixed Cost Program']

    def tearDown(self):
        self.proj = None
        pylab.close('all')
    
#Define Tests to run
class TestProgset(CreateProgsetTest):
    def test_progset_creation(self):
        '''
            - Intended to test whether progsets are created with the correct name
        '''
        progset_name = 'test_creation'
        self.proj.makeProgset(name=progset_name)
        num_progsets = len(self.proj.progsets)
        self.assertEqual(1, num_progsets, 'Incorrect number of progsets exist! Expected Number: 1, Existing Number: %i' %num_progsets)
        self.assertEqual(progset_name, self.proj.progsets[0].name, 'Progset names incorrectly created! Expected name: %s, Actual Name: %s' % (progset_name, self.proj.progsets[0].name))
        return None
    
class TestPrograms(ProgramsTest):
    def test_program_names(self):
        '''
            - Check program names are created properly
            - 
        '''
        for index, prog_name in enumerate(self.prog_short_names):
            self.assertEqual(self.prog_short_names[index], self.proj.progsets[0].progs[index].label, 
                             'Program labels are mismatched! Expected name: %s, Actual name: %s' %(self.prog_short_names[index], self.proj.progsets[0].progs[index].label))
            self.assertEqual(self.prog_full_names[index], self.proj.progsets[0].progs[index].name, 
                             'Program descriptive names are mismatched! Expected name: %s, Actual name: %s' %(self.prog_full_names[index], self.proj.progsets[0].progs[index].name))
        return None
        
    def test_progset_populationgroups(self):
        '''
            - Programs are created for the correct population groups
        '''
        prog_def = {'ds-tb':  ['0-14', '15-49', '50+'], 
                    'mdr-tb': ['15-49', '50+'], 
                    'vac-tb': ['0-14'], 
                    'cure-tb':['Pris'],
                    'fixed': []}

        for index, prog_name in enumerate(self.prog_short_names):
            self.assertListEqual(prog_def[prog_name], self.proj.progsets[0].progs[index].target_pops, 
                                 'Population definitions for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].target_pops))
        return None
    
    def test_progset_budgets(self):
        '''
            - Test whether budgets have been passed into progset correctly:
                - Length of array
                - Content of each array/list
        '''
        prog_def = {'ds-tb':  [7625000., 15250000.], 
                    'mdr-tb': [5100000., np.nan], 
                    'vac-tb': [10500.], 
                    'cure-tb':[0., 4500.],
                    'fixed': [1000000.]}
        
        for index, prog_name in enumerate(self.prog_short_names):
            self.assertEqual(len(prog_def[prog_name]), len(self.proj.progsets[0].progs[index].cost),
                             'Array length imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].cost))
            for indices, budget_value in enumerate(prog_def[prog_name]):
                if np.isnan(budget_value):
                    self.assertEqual(np.isnan(budget_value), np.isnan(self.proj.progsets[0].progs[index].cost[indices]),
                                     'Budget values imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].cost))
                else:
                    self.assertEqual(budget_value, self.proj.progsets[0].progs[index].cost[indices],
                                 'Budget values imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].cost))
        return None
    
    def test_progset_attributes(self):
        '''
            - Parameters are not cross-referenced (i.e. parameters for Program X are not existing in Program Y)
        '''
        pass
    
    def test_progset_attribvalues(self):
        '''
            - Attirbute values are stored within the progset correctly
        '''
        pass
    

    
    def test_costcoverage_curves(self):
        '''
            - Test saturation of curves (saturation not implemented in cost curves yet)
            - Budget to Coverage mapping
            - Coverage to Budget mapping
        '''
        pass
    
    def test_parset_progset(self):
        '''
            - Continuity between parset and progset
        '''
    
if __name__ == '__main__':
    unittest.main()


