from optima_tb.project import Project
import unittest
import os

#Spreadsheets to use
databook = os.path.abspath('../tests/databooks/databook_model_full (with programs).xlsx')
cascade =  os.path.abspath('../tests/cascade_spreadsheet/cascade_model_full (with programs).xlsx')

#Setup Unit Test class
class CreateProgsetTest(unittest.TestCase):
    def setUp(self):
        self.proj= Project(name = 'unittest', cascade_path = self.cascade)
        self.proj.makeParset()

    def tearDown(self):
        self.proj = None
        
class ProgramsTest(unittest.TestCase):
    def setUp(self):
        self.proj= Project(name = 'unittest', cascade_path = self.cascade)
        self.proj.makeParset()
        self.proj.makeProgset()

    def tearDown(self):
        self.proj = None
    
#Define Tests to run
class TestProgset(CreateProgsetTest):
    def test_progset_creation(self):
        '''
            - Intended to test whether progsets are created with the correct name
        '''
        self.proj.makeProgset(name='test_creation')
    
class TestPrograms(ProgramsTest):
    def test_progset_populationgroups(self):
        '''
            - Programs are created for the correct population groups
        '''
        pass
    
    def test_progset_attributes(self):
        '''
            - Parameters are not cross-referenced (i.e. parameters for Program X are not existing in Program Y)
        '''
        pass
    
    
    def test_costcoverage_curves(self):
        '''
            - Test saturation of curves
            - Budget to Coverage mapping
            - Coverage to Budget mapping
        '''
        pass
    
    def test_parset_progset(self):
        '''
            - Continuity between parset and progset
        '''
    
if __name__ == '__main__':
    CreateProgsetTest.cascade = cascade
    CreateProgsetTest.databook = databook
    ProgramsTest.cascade = cascade
    ProgramsTest.databook = databook
    unittest.main()


