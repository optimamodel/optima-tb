from optima_tb.project import Project
from optima_tb.parsing import FunctionParser
from optima_tb.utils import odict
import numpy as np
import pylab
import unittest
import os

#Spreadsheets to use
databook  = os.path.abspath('../tests/databooks/databook_model_full.xlsx')
cascade   =  os.path.abspath('../tests/cascade_spreadsheet/cascade_model_full.xlsx')
TOLERANCE = 1e-6
parser = FunctionParser(debug=False)

def calculateTreatment(effic=None, adher=None, dur=None):
    '''Convenience function to calculate impact value when verifying the getImpact functionality
    '''
    if adher is None:
        val = 1-(1-effic)**(12/dur)
    else:
        val = 1-adher
    return val
    
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
                - Budget created against correct year
        '''
        year_def = {'ds-tb':  [2010., 2015.], 
                    'mdr-tb': [2013., 2015.], 
                    'vac-tb': [2000.], 
                    'cure-tb':[2009., 2015.],
                    'fixed':  [2000.]}
        
        prog_def = {'ds-tb':  [7625000., 15250000.], 
                    'mdr-tb': [5100000., np.nan], 
                    'vac-tb': [10500.], 
                    'cure-tb':[0., 4500.],
                    'fixed': [1000000.]}
        
        for index, prog_name in enumerate(self.prog_short_names):
            self.assertEqual(len(prog_def[prog_name]), len(self.proj.progsets[0].progs[index].cost),
                             'Budget array length imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].cost))
            self.assertListEqual(year_def[prog_name], np.ndarray.tolist(self.proj.progsets[0].progs[index].t),
                                 'Year definitions for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, year_def[prog_name], self.proj.progsets[0].progs[index].t))
            for indices, budget_value in enumerate(prog_def[prog_name]):
                if np.isnan(budget_value):
                    self.assertEqual(np.isnan(budget_value), np.isnan(self.proj.progsets[0].progs[index].cost[indices]),
                                     'Budget values imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].cost))
                else:
                    self.assertEqual(budget_value, self.proj.progsets[0].progs[index].cost[indices],
                                 'Budget values imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].cost))
        return None
    
    def test_progset_coverage(self):
        '''
            - Test whether coverages have been passed into progset correctly:
                - Length of array
                - Coverage format
                - Content of each array/list
        '''
        prog_format = {'ds-tb':  'number', 
                       'mdr-tb': 'number', 
                       'vac-tb': 'fraction', 
                       'cure-tb':'number',
                       'fixed':  None}
        
        prog_def = {'ds-tb':  [76250., 152500.], 
                    'mdr-tb': [25500., np.nan], 
                    'vac-tb': [0.01], 
                    'cure-tb':[0., 0.9],
                    'fixed':  [np.nan]}

        for index, prog_name in enumerate(self.prog_short_names):
            self.assertEqual(len(prog_def[prog_name]), len(self.proj.progsets[0].progs[index].cov),
                             'Coverage array length imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].cov))
            self.assertEqual(prog_format[prog_name], self.proj.progsets[0].progs[index].cov_format,
                             'Coverage format imported for program: %s is incorrect! Expected Type: %s, Actual Type: %s' %(prog_name, prog_format[prog_name], self.proj.progsets[0].progs[index].cov_format))
            
            for indices, coverage_value in enumerate(prog_def[prog_name]):
                if np.isnan(coverage_value):
                    self.assertEqual(np.isnan(coverage_value), np.isnan(self.proj.progsets[0].progs[index].cov[indices]),
                                     'Coverage values imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].cov))
                else:
                    self.assertEqual(coverage_value, self.proj.progsets[0].progs[index].cov[indices],
                                 'Coverage values imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].cov))
        return None
    
    def test_progset_attributes(self):
        '''
            - Parameters are not cross-referenced (i.e. parameters for Program X are not existing in Program Y)
            - Check for attribute duplication
            NOTE: odict check is not rigourous enough, needs to be better incorporated
        '''
        prog_def = {'ds-tb':  ['effic', 'adher', 'dur'], 
                    'mdr-tb': ['effic', 'adher', 'dur'], 
                    'vac-tb': ['effic'], 
                    'cure-tb':['effic'],
                    'fixed':  odict()}
        for index, prog_name in enumerate(self.prog_short_names):
            for attrib_name in prog_def[prog_name]:
                self.assertEqual(True, attrib_name in self.proj.progsets[0].progs[index].attributes.keys(),
                                 'Attribute key: %s was not found in Program: %s. Expected Attributes: %s, Actual Attributes: %s' %(attrib_name, prog_name, prog_def[prog_name], self.proj.progsets[0].progs[index].attributes.keys()))
        return None
    
    def test_progset_attribvalues(self):
        '''
            - Attribute values are stored within the progset correctly
        '''
        prog_def = {'ds-tb':  {'effic': [0.75, 0.9], 'adher': [0.66666667, 0.9], 'dur': [12., 6.]}, 
                    'mdr-tb': {'effic': [np.nan, 2./3.], 'adher': [0.5, np.nan], 'dur': [18., 18.]}, 
                    'vac-tb': {'effic': [0.5]}, 
                    'cure-tb':{'effic': [0.2, 0.2]},
                    'fixed':  odict()}
        
        for index, prog_name in enumerate(self.prog_short_names):
            for attrib_name in prog_def[prog_name].keys():
                for indices, attrib_value in enumerate(prog_def[prog_name][attrib_name]):
                    if np.isnan(attrib_value):
                        self.assertEqual(np.isnan(attrib_value), np.isnan(self.proj.progsets[0].progs[index].attributes[attrib_name][indices]),
                                         'Attribute values loaded for "%s" and imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(attrib_name, prog_name, prog_def[prog_name][attrib_name], self.proj.progsets[0].progs[index].attributes[attrib_name]))
                    else:
                        difference = abs(attrib_value - self.proj.progsets[0].progs[index].attributes[attrib_name][indices])
                        self.assertLess(difference, TOLERANCE,
                                         'Attribute values loaded for "%s" and imported for program: %s is incorrect! Expected List: %s, Actual List: %s' %(attrib_name, prog_name, prog_def[prog_name][attrib_name], self.proj.progsets[0].progs[index].attributes[attrib_name]))
        return None

    def test_budgetToCoverage_mapping(self):
        '''
            - Budget to Coverage mapping
        '''
        budget_values = [0., 100000., 500000., 1000000.]
        unit_costs =    {'ds-tb'  : 100., 
                         'mdr-tb' : 200., 
                         'vac-tb' : 10.5, 
                         'cure-tb': 50.,
                         'fixed'  : np.nan}
        for index, prog_name in enumerate(self.prog_short_names):
            for budget in budget_values:
                if self.proj.progsets[0].progs[index].cov_format == 'number':
                    expected_value = budget / unit_costs[prog_name]
                else:
                    expected_value = 0.01 * budget / unit_costs[prog_name]
                    
                model_value = self.proj.progsets[0].progs[index].getCoverage(budget)
                
                if self.proj.progsets[0].progs[index].func_specs['type'] == 'cost_only':
                    self.assertEqual(np.isnan(expected_value), np.isnan(model_value),
                                'Coverage value for program %s was incorrectly calculated. Expected value: %s, Model value: %s' %(prog_name, expected_value, model_value))
                elif self.proj.progsets[0].progs[index].func_specs['type'] == 'linear':
                    self.assertEqual(expected_value, model_value,
                                'Coverage value for program %s was incorrectly calculated. Expected value: %s, Model value: %s' %(prog_name, expected_value, model_value))
                else: print "Incorrect entries detected!"
                    
        return None
    
    def test_coverageToBudget_mapping(self):
        '''
            - Budget to Coverage mapping
        '''
        coverage_values = [0., 100., 1667., 1000000.]
        unit_costs =    {'ds-tb'  : 100., 
                         'mdr-tb' : 200., 
                         'vac-tb' : 10.5, 
                         'cure-tb': 50.,
                         'fixed'  : np.nan}
        for index, prog_name in enumerate(self.prog_short_names):
            if prog_name == 'fixed': continue
            for coverage in coverage_values:
                if self.proj.progsets[0].progs[index].cov_format == 'number':
                    expected_value = coverage * unit_costs[prog_name]
                else:
                    expected_value = coverage * unit_costs[prog_name] / 0.01
                    
                model_value = self.proj.progsets[0].progs[index].getBudget(coverage)
                
                if self.proj.progsets[0].progs[index].func_specs['type'] == 'linear':
                    self.assertEqual(expected_value, model_value,
                                'Budget value for program %s was incorrectly calculated. Expected value: %s, Model value: %s' %(prog_name, expected_value, model_value))
                else: print "Incorrect entries detected!"
                    
        return None
    
    def test_getDefaultBudget(self):
        '''
            - Test interpolation of budget
            - Test extrapolation of budget
            - Test whether correct budget is returned for a given year
        '''
        test_years = [2000, 2005, 2010, 2015]
        expected_value = {'ds-tb'  : {'2000': 7625000., '2005': 7625000., '2010': 7625000., '2015': 15250000.},
                          'mdr-tb' : {'2000': 5100000., '2005': 5100000., '2010': 5100000., '2015':  5100000.},
                          'vac-tb' : {'2000':   10500., '2005':   10500., '2010':   10500., '2015':    10500.},
                          'cure-tb': {'2000':       0., '2005':       0., '2010':     750., '2015':     4500.},
                          'fixed'  : {'2000': 1000000., '2005': 1000000., '2010': 1000000., '2015':  1000000.}
                         }
        for year in test_years:
            for prog in self.proj.progsets[0].prog_ids:
                prog_index = self.proj.progsets[0].prog_ids[prog]
                year_budgetValue = self.proj.progsets[0].progs[prog_index].getDefaultBudget(year=year)
                difference = abs(year_budgetValue - expected_value[prog][str(year)]) #taking difference so that floating points can be compared
                self.assertLess(difference, TOLERANCE, 
                                'Budget value was incorrectly imported or there was some issue with interpolation!\nProgram: %s, Year: %i\n Model Value: %g, Expected Value: %g'%(prog,year, year_budgetValue, expected_value[prog][str(year)]))
        return None
    
    def test_getImpact(self):
        '''
            - Test impact of a program
            - Test whether impact saturates at some value
        '''
        test_years = [2000, 2005, 2012, 2015]
        expected_value = {'ds-tb'  : {'adtreat_prop'   : {'2000': 76250., '2005': 76250., '2012':106750., '2015': 152500.},
                                      'adtorec_prop'   : {'2000': 76250.*calculateTreatment(effic=0.75, dur=12.), '2005': 76250.*calculateTreatment(effic=0.75, dur=12.), '2012':106750.*calculateTreatment(effic=0.81, dur=9.6), '2015': 152500.*calculateTreatment(effic=0.9, dur=6.)},
                                      'adtomtreat_prop': {'2000': 76250.*calculateTreatment(adher=2./3.), '2005': 76250.*calculateTreatment(adher=2./3.), '2012':106750.*calculateTreatment(adher=0.76), '2015': 152500.*calculateTreatment(adher=0.9)},
                                      },
                          
                          'mdr-tb' : {'adtreat_prop'   : {'2000': 25500.*0.85, '2005': 25500.*0.85, '2012': 25500.*0.85, '2015':  25500.*0.85},
                                      'adtorec_prop'   : {'2000': 25500.*calculateTreatment(effic=2./3., dur=18.), '2005': 25500.*calculateTreatment(effic=2./3., dur=18.), '2012': 25500.*calculateTreatment(effic=2./3., dur=18.), '2015':  25500.*calculateTreatment(effic=2./3., dur=18.)},
                                      'adtomtreat_prop': {'2000': 25500.*calculateTreatment(adher=0.5), '2005': 25500.*calculateTreatment(adher=0.5), '2012': 25500.*calculateTreatment(adher=0.5), '2015':  25500.*calculateTreatment(adher=0.5)},
                                      },
                                      
                          'vac-tb' : {'v_rate'         : {'2000': 0.01*10500./10.5, '2005': 0.01*10500./10.5, '2012': 0.01*10500./10.5, '2015':  0.01*10500./10.5}},
                                      
                          'cure-tb': {'cure_rate'      : {'2000':     0., '2005':     0., '2012':    45., '2015':     90.}},
                                      
                          'fixed'  : {'cure_rate'      : {'2000':     0., '2005':     0., '2012':     0., '2015':      0.}},
                         }
        for year in test_years:
            for prog in self.proj.progsets[0].prog_ids:
                prog_index = self.proj.progsets[0].prog_ids[prog]
                year_budgetValue = self.proj.progsets[0].progs[prog_index].getDefaultBudget(year=year)
                for impact_label in self.proj.progsets[0].progs[prog_index].target_pars:
                    year_impactValue = self.proj.progsets[0].progs[prog_index].getImpact(year_budgetValue, parser=parser, impact_label=impact_label, years=[year])
                    calculated_impactValue = expected_value[prog][impact_label][str(year)]
                    difference = abs(year_impactValue - calculated_impactValue) #taking difference so that floating points can be compared
                    self.assertLess(difference, TOLERANCE, 
                                    'Impact value was incorrectly calculated!\nProgram: %s, Year: %i, Impact Parameter: %s\nModel Value: %g, Expected Value: %g'%(prog,year, impact_label, year_impactValue, calculated_impactValue))
        return None
    
    def test_costcoverage_curves(self):
        '''
            - Test saturation of curves (saturation not implemented in cost curves yet)
        '''
        pass
    
if __name__ == '__main__':
    unittest.main()


