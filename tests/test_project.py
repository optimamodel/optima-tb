from optima_tb.project import Project
from optima_tb.dataio import importObj
from optima_tb.utils import odict
import pylab
import unittest
import numpy as np
import os

#Initialize parameters
databook = '../tests/databooks/databook_model_full.xlsx'
cascade =  '../tests/cascade_spreadsheet/cascade_model_full.xlsx'
data = ['pops', 'contacts', 'transfers', 'characs', 'linkpars']
attributes = ['data', 'parsets', 'results', 'scenarios']
pars = ['pop_names', 'pop_labels', 'pars', 'par_ids', 'transfers', 'contacts']

class ModelTest(unittest.TestCase):
    def setUp(self):
        self.proj= Project(name = 'unittest', cascade_path = self.cascade)

    def tearDown(self):
        self.proj = None
    
class TestProject(ModelTest):
    def test_project_creation(self):
        '''
            1. Project is of correct Class type
            2. Project has all attributes
            3. Project attributes are of type odict and empty upon creation
        '''
        self.assertIsInstance(self.proj, Project, 'Project is not of instance type optima_tb.project')
        
        #Check existence of all attributes and then check whether they are empty odicts
        for attribute_type in self.attributes:
            self.assertTrue(hasattr(self.proj, attribute_type), 'Attribute "%s" does not exist in Project' % attribute_type)
            temp = getattr(self.proj, attribute_type)
            self.assertIsInstance(temp, odict, 'Attribute "%s" is not of type odict' % attribute_type)
            self.assertEqual(temp.keys(), [], 'Attribute "%s" of project is a non-empty odict' % attribute_type)
        return None
    
    def test_loadspreadsheet(self):
        '''
            1. All data types exist within project.data
            2. No other entries are created  when loadspreadsheet is called (e.g. resultset, scenario)
        '''
        self.proj.loadSpreadsheet(databook_path=databook)
        
        #Check data exists and is loaded into the project from spreadsheet
        for data_type in self.data:
            self.assertIn(data_type, self.proj.data.keys(), 'Key "%s" does not exist in Project' % data_type)
            
        #Secondary checks to ensure that no other data was created
        for attribute_type in self.attributes:
            self.assertTrue(hasattr(self.proj, attribute_type), 'Attribute "%s" does not exist in Project' % attribute_type)
            if attribute_type == 'data': continue
            else: 
                temp = getattr(self.proj, attribute_type)
                self.assertIsInstance(temp, odict, 'Attribute "%s" is not of type odict' % attribute_type)
                self.assertEqual(temp.keys(), [], 'Attribute "%s" of project is a non-empty odict' % attribute_type)
        return None
    
    def test_makeresetParset(self):
        '''
            1. Checks whether multiple parametersets are created or not
            2. Parameter sets are of correct length
            3. Parameter mapping checks (i.e. par_ids map to pars['cascade'] and pars['linkpars'])
            4. No other entries are created  when loadspreadsheet is called (e.g. resultset, scenario)
            5. Test whether reset parset manages to ensure all parsets are set as empty odicts
        '''
        self.proj.loadSpreadsheet(databook_path=databook)
        self.assertEqual(0, len(self.proj.parsets), 'A parset already exists even though makeParset was not called')
        
        #Test makeParset
        self.proj.makeParset()
        self.assertEqual(1, len(self.proj.parsets),'Multiple parsets (or no parset) exists when only one should exist')
        self.assertEqual(len(self.proj.parsets[-1].pop_labels), len(self.proj.parsets[-1].pop_names), 'Number of Population groups does not match between "pop_labels" and "pop_names"')
        for par_type in self.pars:
            self.assertTrue(hasattr(self.proj.parsets[-1], par_type), 'Key "%s" does not exist in parset' % par_type)
        for key in self.proj.parsets[-1].par_ids:
            for pars in self.proj.parsets[-1].par_ids[key]:
                self.assertTrue(pars == self.proj.parsets[-1].pars[key][self.proj.parsets[-1].par_ids[key][pars]].label, 'Parameter "%s" under parameter id "%s" does not exist at index "%i" as mentioned in parset.par_ids' % (pars, key, self.proj.parsets[-1].par_ids[key][pars]))
        #Secondary checks to ensure that no other data was created
        for attribute_type in self.attributes:
            self.assertTrue(hasattr(self.proj, attribute_type), 'Attribute "%s" does not exist in Project' % attribute_type)
            if attribute_type == 'data' or attribute_type == 'parsets': continue
            else: 
                temp = getattr(self.proj, attribute_type)
                self.assertIsInstance(temp, odict, 'Attribute "%s" is not of type odict' % attribute_type)
                self.assertEqual(temp.keys(), [], 'Attribute "%s" of project is a non-empty odict' % attribute_type)
        
        #Test resetParset
        self.proj.resetParsets()
        self.assertEqual(0, len(self.proj.parsets), 'ResetParsets() was called, but parset(s) still exist')
        self.assertEqual(self.proj.parsets.keys(), [], 'Attribute "parsets" of project is a non-empty odict after calling resetparsets()')
        return None
    
    def test_basic_project_settings(self):
        #Test year settings
        start_year = 2000.0
        end_year = 2030.0
        observed_end = 2020.0
        self.proj.setYear([start_year, end_year], observed_data=False)
        self.proj.setYear([start_year, observed_end], observed_data=True)
        self.assertEqual(self.proj.settings.tvec_observed_end, observed_end, 'Observed end year was not setup correctly, indicated start year is "%f", current value is: "%f"' % (self.proj.settings.tvec_observed_end, observed_end))
        self.assertEqual(self.proj.settings.tvec_start, start_year, 'Start year was not setup correctly, indicated start year is "%f", current value is: "%f"' % (self.proj.settings.tvec_start, start_year))
        self.assertEqual(self.proj.settings.tvec_end, end_year, 'End year was not setup correctly, indicated start year is "%f", current value is: "%f"' % (self.proj.settings.tvec_end, end_year))
        self.assertIsNotNone(self.proj.settings.tvec_dt, 'dt value is: %s' % self.proj.settings.tvec_dt)
        self.assertGreater(self.proj.settings.tvec_dt, 0., 'Error zero or negative value encountered for dt - value detected: %f' % (self.proj.settings.tvec_dt))
        return None
    
    def test_runSim(self):
        self.proj.loadSpreadsheet(databook_path=databook)
        self.proj.makeParset()
        try:
            plot_over = [self.proj.settings.tvec_start, self.proj.settings.tvec_end]
            self.proj.settings.plot_settings['x_ticks'] = [np.arange(plot_over[0],plot_over[1]+self.proj.settings.tvec_dt,5,dtype=int),np.arange(plot_over[0],plot_over[1]+self.proj.settings.tvec_dt,5,dtype=int)]
            self.proj.settings.plot_settings['xlim'] = (plot_over[0]-0.5,plot_over[1]+0.5)
            self.proj.makeParset()
            results = self.proj.runSim(plot=True)
        except: results = None
        self.assertIsNotNone(results, 'runSim failed/did not complete, check log for more details')
        pylab.close('all')
        return None

    def test_import_export_parset(self) :
        self.proj.loadSpreadsheet(databook_path=databook)
        self.proj.makeParset()
        parset_name = self.proj.parsets[-1].name
        default = 'default'
        new_name = 'imported_parset'
        self.proj.exportParset(parset_name=parset_name)
        self.proj.importParset(parset_filename=parset_name+'.csv', new_name=new_name)
        
        #Run tests on both parsets
        self.assertEqual(2, len(self.proj.parsets), 'Number of parsets is greater than two')
        for par_type in self.pars:
            self.assertTrue(hasattr(self.proj.parsets[default], par_type), 'Key "%s" does not exist in parset: %s' % (par_type, default))
            self.assertTrue(hasattr(self.proj.parsets[new_name], par_type), 'Key "%s" does not exist in parset: %s' % (par_type, new_name))
            temp  = getattr(self.proj.parsets[default], par_type)
            temp2 = getattr(self.proj.parsets[new_name], par_type)
            if par_type == 'pop_names' or par_type == 'pop_labels':
                self.assertListEqual(temp, temp2, 'Attribute "%s" Population Groups do not match between imported and exported parsets' % par_type)
            elif par_type == 'par_ids':
                self.assertDictEqual(temp, temp2, 'Attribute "%s" are not matching up between the imported/exported parsets' % par_type)
            elif par_type == 'pars':
                for key in temp:
                    for par in range(len(temp[key])):
                        for subpar in temp[key][par].__dict__:
                            newtemp  = getattr(temp[key][par], subpar)
                            newtemp2 = getattr(temp2[key][par], subpar)
                            if subpar == 'label':
                                self.assertEqual(newtemp, newtemp2, 'Parameter labels are mismatched between imported and exported parsets for par_type "%s". key "%s", par "%s". subpar "%s"' % (par_type, key, par, subpar))
                            else:
                                for popkey in newtemp:
                                    if isinstance(newtemp[popkey], str):
                                        self.assertEqual(newtemp[popkey], newtemp2[popkey], 'Parameter are mismatched between imported and exported parsets for par_type "%s". key "%s", par "%s". subpar "%s", popkey "%s"' % (par_type, key, par, subpar, popkey))
                                    elif isinstance(newtemp[popkey], np.ndarray):
                                        x = newtemp[popkey].tolist()
                                        y = newtemp2[popkey].tolist()
                                        if np.isnan(x) and np.isnan(y): continue
                                        elif np.isinf(x) and np.isinf(y): continue
                                        self.assertListEqual(x, y, 'Parameter are mismatched between imported and exported parsets for par_type "%s". key "%s", par "%s". subpar "%s", popkey "%s"' % (par_type, key, par, subpar, popkey))
                                    elif isinstance(newtemp[popkey], float):
                                        self.assertEqual(float(newtemp[popkey]), float(newtemp2[popkey]), 'Parameter are mismatched between imported and exported parsets for par_type "%s". key "%s", par "%s". subpar "%s", popkey "%s"' % (par_type, key, par, subpar, popkey))
        os.remove(parset_name+'.csv')
        return None
    
    def test_calculateFit(self):
        self.proj.loadSpreadsheet(databook_path=databook)
        self.proj.makeParset()
        results = self.proj.runSim()
        score = self.proj.calculateFit(results)
        self.assertGreaterEqual(max(score), 0., 'Fit Score is "%s" for selected datasheet' % score)
        return None
    
    def test_exportProject(self):
        self.proj.loadSpreadsheet(databook_path=databook)
        self.proj.makeParset()
        filename = self.proj.exportProject()
        proj2 = importObj(filename=filename)
        for key in self.proj.settings.linkpar_specs:
            for key2 in self.proj.settings.linkpar_specs[key]:
                self.assertEqual(self.proj.settings.linkpar_specs[key][key2], proj2.settings.linkpar_specs[key][key2], '%s is not matching up for %s between the two projects'%(key, key2))
        for key in self.proj.settings.charac_specs:
            for key2 in self.proj.settings.charac_specs[key]:
                self.assertEqual(self.proj.settings.charac_specs[key][key2], proj2.settings.charac_specs[key][key2], 'Validation failed! %s is not matching up for %s between the two projects'%(key,key2))
        os.remove(filename)
        return None
    
    def test_makeManualCalibration(self):
        self.proj.loadSpreadsheet(databook_path=databook)
        self.proj.makeParset()
        results1 = self.proj.runSim()
        score1 = self.proj.calculateFit(results1)
        
        dict_change_params = {'v_rate': [0.9],
                              'b_rate' : [0]}
        dict_change_params2 = {'phi_early': [0.9]}

        rate_dict = {'0-14' : dict_change_params,
                     '15-49' : dict_change_params2}
        
        pname2 = 'test_Manual_calibration'
        self.proj.makeManualCalibration(parset_name=pname2, rate_dict=rate_dict)
        results2 = self.proj.runSim(parset_name=pname2)
        score2 = self.proj.calculateFit(results2)        
        self.assertNotEqual(max(score1), max(score2), 'Fit Scores should not be the same since parameter rates are different!')
        return None
    
    def test_runAutofitCalibration(self):
        self.proj.loadSpreadsheet(databook_path=databook)
        self.proj.makeParset()
        return None
    
    

if __name__ == '__main__':
    ModelTest.cascade = cascade
    ModelTest.databook = databook
    ModelTest.data = data
    ModelTest.attributes = attributes
    ModelTest.pars = pars
    
    unittest.main()


