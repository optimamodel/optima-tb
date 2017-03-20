#%% Imports
import logging
import logging.config

logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()


from optima_tb.utils import tic, toc, odict, OptimaException
from optima_tb.model import runModel
from optima_tb.settings import Settings 
from optima_tb.parameters import ParameterSet, export_paramset, load_paramset
from optima_tb.programs import ProgramSet
from optima_tb.plotting import plotProjectResults
from optima_tb.databook import makeSpreadsheetFunc, loadSpreadsheetFunc
from optima_tb.optimization import optimizeFunc
from optima_tb.calibration import makeManualCalibration, calculateFitFunc, performAutofit
from optima_tb.scenarios import ParameterScenario, BudgetScenario, CoverageScenario
from optima_tb.dataio import exportObj, importObj

from uuid import uuid4 as uuid
from numpy import max


#%% Project class (i.e. one self-contained geographical unit)

class Project(object):
    ''' The main Optima project class. Almost all Optima functionality is provided by this class. '''

    def __init__(self, name = 'default', cascade_path = '../data/cascade.xlsx', **args):
        ''' Initialize project. '''

        self.name = name
        self.uid = uuid()
        
        self.settings = Settings(cascade_path = cascade_path, **args)        
        self.data = odict()

        self.parsets = odict()
        self.progsets = odict()
        self.results = odict()
        
        self.scenarios = odict()
        
        logger.info("Created project: %s"%self.name)
        
    def resetParsets(self):
        ''' Convenience function called externally to delete all parsets. '''
        self.parsets = odict()
        
    def setYear(self, yearRange,observed_data=True):
        '''
        
        @param yearRange: tuple or list 
        @param observed_data: bool indicating whether to change observed date end or simulation date end 
        '''
        self.settings.tvec_start = yearRange[0]
        if observed_data:
            self.settings.tvec_observed_end = yearRange[1]
        else:
            self.settings.tvec_end = yearRange[1]
    
    
    def runSim(self, parset = None, parset_name = 'default', progset = None, progset_name = None, options = None, plot = False, debug = False):
        ''' Run model using a selected parset and store/return results. '''
        
        if parset is None:
            if len(self.parsets) < 1: 
                raise OptimaException('ERROR: Project "%s" appears to have no parameter sets. Cannot run model.' % self.name)
            else:
                try: parset = self.parsets[parset_name]
                except: raise OptimaException('ERROR: Project "%s" is lacking a parset named "%s". Cannot run model.' % (self.name, parset_name))
                
        if progset is None:
            try: progset = self.progsets[progset_name]
            except: logger.info('Initiating a standard run of project "%s" (i.e. without the influence of programs).' % self.name)
        if progset is not None:
            if options is None:
                logger.info('Program set "%s" will be ignored while running project "%s" due to no options specified.' % (progset.name, self.name))
                progset = None

        tm = tic()

        #results = runModel(settings = self.settings, parset = parset)
        results = runModel(settings = self.settings, parset = parset, progset = progset, options = options)

        toc(tm, label = 'running %s model' % self.name)
        
        if plot:
            tp = tic()
            self.plotResults(results = results, debug = debug)
            toc(tp, label = 'plotting %s' % self.name)
        
        return results
        
        
    def optimize(self, parset = None, parset_name = 'default', progset = None, progset_name = 'default', options = None):
        ''' Optimize model using a selected parset and store/return results. '''
        
        if parset is None:
            if len(self.parsets) < 1: 
                raise OptimaException('ERROR: Project "%s" appears to have no parameter sets. Cannot optimize model.' % self.name)
            else:
                try: parset = self.parsets[parset_name]
                except: raise OptimaException('ERROR: Project "%s" is lacking a parset named "%s". Cannot optimize model.' % (self.name, parset_name))        
        
        if progset is None:
            if len(self.progsets) < 1: 
                raise OptimaException('ERROR: Project "%s" appears to have no program sets. Cannot optimize model.' % self.name)
            else:
                try: progset = self.progsets[progset_name]
                except: raise OptimaException('ERROR: Project "%s" is lacking a progset named "%s". Cannot optimize model.' % (self.name, progset_name))
                
        results = optimizeFunc(settings = self.settings, parset = parset, progset = progset, options = options)
        
        return results
        

    def plotResults(self, results, colormappings=None, colorlabels=None, debug=False, pop_labels=None, plot_observed_data=True,savePlot=False,figName=None):
        ''' Plot all available results '''

        plotProjectResults(results, settings=self.settings, data=self.data, title = self.name.title(), colormappings=colormappings, colorlabels=colorlabels, pop_labels=pop_labels, debug = debug, plot_observed_data=plot_observed_data, save_fig=savePlot, fig_name=figName)
            
    
    
    def makeSpreadsheet(self, databook_path = None, num_pops = 5, num_migrations = 2, num_progs = 0):
        ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''
        
        if databook_path is None: databook_path = '../data/' + self.name + '-data.xlsx'
        logger.info("Attempting to create databook %s"%databook_path)
        
        makeSpreadsheetFunc(settings = self.settings, databook_path = databook_path, num_pops = num_pops, num_migrations = num_migrations, num_progs = num_progs)        
        
    
    def loadSpreadsheet(self, databook_path = None):
        ''' Load data spreadsheet into Project data dictionary. '''
        if databook_path is None: databook_path = '../data/' + self.name + '-data.xlsx'
        logger.info("Attempting to load databook %s"%databook_path)
        
        self.data = loadSpreadsheetFunc(settings = self.settings, databook_path = databook_path) 
        


    def makeParset(self, name = 'default'):
        ''' Transform project data into a set of parameters that can be used in model simulations. '''

        if not self.data: raise OptimaException('ERROR: No data exists for project "%s".' % self.name)
        self.parsets[name] = ParameterSet(name = name)
        self.parsets[name].makePars(self.data)
        
    def makeProgset(self, name = 'default'):
        ''' Transform project data into a set of programs that can be used in budget scenarios and optimisations. '''

        if not self.data: raise OptimaException('ERROR: No data exists for project "%s".' % self.name)
        self.progsets[name] = ProgramSet(name = name)
        self.progsets[name].makeProgs(data = self.data, settings = self.settings)
        
    def exportParset(self, parset_name):
        ''' Exports parset to .csv file '''
        if not parset_name in self.parsets.keys():
            raise OptimaException("ERROR: no parameter set '%s' found"%parset_name)
        export_paramset(self.parsets[parset_name])
        
    def importParset(self, parset_filename, new_name=None):
        ''' Imports parameter set from .csv file '''
        paramset = load_paramset(parset_filename)
        if new_name is None:
            new_name = paramset.name
        else:
            paramset.name = new_name
            
        if not new_name in self.parsets.keys():
            logger.info("Imported new parameter set: %s"%new_name)
        else:
            logger.info("Imported and overwriting parameter set: %s"%new_name)
        self.parsets[new_name] = paramset
        
        
    def makeManualCalibration(self, parset_name, rate_dict):
        ''' Take in dictionary of updated values for rates that can be used for manual calibration
            Update values dict: {pop_name : {parameter : value} '''
        if not parset_name in self.parsets.keys():
            self.makeParset(name=parset_name)
        paramset = self.parsets[parset_name]
        logger.info("Updating parameter values in parset=%s"%(parset_name))
        
        makeManualCalibration(paramset,rate_dict)
    
    
    def calculateFit(self, results, metric=None):
        '''
        Calculates the score for the fit during manual calibration and prints to output. 
        
        Params:
            results    resultSet object
            metric     type of metric used (defaults to default specified in settings).
        
        Future: can consider saving results into resultset and accessing results by name / parset name
        '''
        if metric is None:
            metric = self.settings.fit_metric
        
        if results is None:
            raise OptimaException('ERROR: no result is specified. Cannot calculate fit.')
        
        if self.data is None or len(self.data)==0:
            raise OptimaException('ERROR: no data is specified. Cannot calculate fit.')
        
        datapoints, _, _ = results.getCharacteristicDatapoints()
        score = calculateFitFunc(datapoints,results.t_observed_data,self.data['characs'],metric)
        logger.info("Calculated scores for fit using %s: largest value=%.2f"%(metric,max(score)))
      
      
    def runAutofitCalibration(self, new_parset_name = None, old_parset_name="default", target_characs=None):
        """
        Runs the autofitting calibration routine, as according to the parameter settings in the 
        settings.autofit_params configuration.
        
        Params:
            new_parset_name    name to save the resulting autofit to
            old_parset_name    name of the parset to use as a base. Default value="default"
        """
        
        if not old_parset_name in self.parsets.keys():
            self.makeParset(name=old_parset_name)
        paramset = self.parsets[old_parset_name]
        
        if new_parset_name is None:
            # TODO: check that autofit doesn't already exist; if so, add suffix
            new_parset_name = "autofit" 
        
        logger.info("About to run autofit on parameters using parameter set = %s"%old_parset_name)
        new_parset = performAutofit(self,paramset,new_parset_name=new_parset_name,target_characs=target_characs,**self.settings.autofit_params)
        logger.info("Created new parameter set '%s' using autofit"%new_parset_name)
        self.parsets[new_parset_name] = new_parset
        
        
    def createScenarios(self, scenario_dict):
        """
        Creates the scenarios to be run, and adds them to this project's store
        of available scenarios to run. Each scenario is described as a (key, value)
        pair in scenario_dict, with key = scenario name. 
        
        Each dictionary value describing a Scenario contains required fields in 
        addition to optional fields.
            Required fields: 
                "type"             "Parameter", "Budget" or "Coverage"
                "scenario_values"   odict of values for scenario, used in initialising the corresponding Scenario object
            Optional fields:
                "run_scenario"     bool : indicating whether this scenario should be run. Default value if unspecified: False
                "overwrite"        bool indicating whether this scenario's values should replace (overwrite = True) or
                                    be added to (overwrite = False) the 'business-as-usual' simulation.
        
        Params:
            scenario_dict =  { name: {
                                     "type" : "Parameter",
                                     "run_scenario" : False,
                                     "overwrite" : True,
                                     "scenario_values" : {...}} , 
                                name: {
                                     "type" : "Budget",
                                     "run_scenario" : True,
                                     "overwrite" : False,
                                     "scenario_values" : {...}} , 
                                     
                                 }
        
        Returns:
            none
            
        
        Parameter Scenario Example: 
            scvalues = odict()
            param = 'birth_transit'
            scvalues[param] = odict()
            scvalues[param]['Pop1'] = odict()
            scvalues[param]['Pop1']['y'] = [3e6, 1e4, 1e4, 2e6]
            scvalues[param]['Pop1']['t'] = [2003.,2004.,2014.,2015.]
            scvalues[param]['Pop1']['y_format'] = 'number'
            scvalues[param]['Pop1']['y_factor'] = DO_NOT_SCALE
            scen_values = { 'test_scenario': {'type': 'Parameter',
                                  'run_scenario' : True,
                                  'scenario_values': scvalues}
               }
            proj = Project(name = 'sampleProject', cascade_path = 'data/cascade-simple.xlsx')
            proj.createScenarios(scen_values)
            
        Budget Scenario Example: 
            budget_options = {'HT-DS': 4e6,'SAT-DS':0,'HT-MDR': 3e4}
            scen_values = { 'test_scenario': {'type': 'Budget',
                                      'overwrite' : True, # it will overwrite scenario to the parset
                                      'run_scenario' : True,
                                      'scenario_values': budget_options}
                   }
            proj = Project(name = 'sampleProject', cascade_path = 'data/cascade-simple.xlsx')
            proj.makeParset(name = 'default_parset')
            proj.makeProgset(name = 'default_progset')
            proj.createScenarios(scen_values)
            resultset = proj.runScenarios(original_parset_name = 'default_parset',
                                  original_progset_name='default_progset',
                                  original_budget_options=options,
                                  include_bau=False)

  
        """
        logger.info("About to create scenarios")
        
        pop_labels = self.data['pops']['label_names']

        
        for scenario_name in scenario_dict.keys():
            vals = scenario_dict[scenario_name]
            
            if scenario_name in self.scenarios.keys():
                logger.warn("Attempting to add scenario '%s' to project %s, that already contains scenario of the same name. Will ignore."%(scenario_name,self.name))
                # TODO decide what to do if scenario with same name already exists. Update or ignore? SJ: prefer to ignore.
            
            if vals['type'].lower() == 'parameter':
                self.scenarios[scenario_name] = ParameterScenario(name=scenario_name,settings=self.settings,pop_labels=pop_labels,**vals)

            elif vals['type'].lower() == 'budget':
                self.scenarios[scenario_name] = BudgetScenario(name=scenario_name,pop_labels=pop_labels,**vals)
            
            elif vals['type'].lower() == 'coverage':
                self.scenarios[scenario_name] = CoverageScenario(name=scenario_name,pop_labels=pop_labels,**vals)

            else:
                raise NotImplementedError("ERROR: no corresponding Scenario type for scenario=%s"%scenario_name)
        
        logger.info("Successfully created scenarios")
        
        

    def runScenarios(self,original_parset_name,original_progset_name=None,original_budget_options=None,scenario_set_name=None,include_bau=False,plot=False,save_results=False):
        """
        Runs scenarios that are contained in this project's collection of scenarios (i.e. self.scenarios). 
        For each scenario run, using original_parset_name, the results generated are saved and 
        returns as a dictionary of results. 
        
        Optional flag include_bau indicates whether to additionally run a "business-as-usual" scenario
        i.e. no external changes. 
        
        
        # TODO implement ability to specify list of scenario names that will run if valid and regardless of scenario.run_scenario
        
        
        Params:
            original_parset_name    name of parameterSet to be used
            include_bau             bool indicating whether to include BAU (business as usual)
            plot                    bool flag indicating whether to plot results as we go
            
        Returns:
            results    dictionary of results obtained for each scenario, with key = scenario_name
        """
        orig_parset = self.parsets[original_parset_name]
        
        if original_progset_name is not None:
            orig_progset = self.progsets[original_progset_name]
        else:
            orig_progset = None
            
        if original_budget_options is None:
            original_budget_options = {}
            
        results = odict()
        
        if include_bau:
            results['BAU'] = self.runSim(parset_name = original_parset_name,progset=orig_progset,options=original_budget_options,plot=plot)

        
        for scen in self.scenarios.keys():
            if self.scenarios[scen].run_scenario:
                scen_name = 'scenario_%s'%self.scenarios[scen].name

                progset, budget_options = self.scenarios[scen].getScenarioProgset(orig_progset,original_budget_options)
            
                results[scen_name] = self.runSim(parset = self.scenarios[scen].getScenarioParset(orig_parset), progset=progset, options=budget_options, parset_name = scen_name, plot=plot)
                
                if scenario_set_name is None:
                    results[scen_name].name = '%s'%(scen_name)
                else:
                    results[scen_name].name = '%s:%s'%(scenario_set_name,scen_name)
                    
                if save_results:
                    results[scen_name].export()
                    export_paramset(self.scenarios[scen].getScenarioParset(orig_parset))
        
        return results
    
    
    
    
    def exportProject(self, filename=None, format='json', compression='zlib'):
        """
        
        This currently saves everything within a project, including results.
        
        Params:
            filename      filename to save to. If none is supplied, value is set to "<project.name>.project"
            format        string for supported format types (json)
            compression   string for supported compression types (zlib)
        
        Usage
            project = Project(name="sample", cascade="cascade.xlsx")
            project.exportProject()
            # saves to "sample.project.Z"
            project.exportProject(filename="special")
            # saves to "special.Z"
        
        """
        if filename is None:
            filename = "%s.project"%self.name
        
        logger.info("Attempting to save file in format=%s"%format)
        filename = exportObj(self,filename=filename,format=format,compression=compression)
        logger.info("Saved to file: %s"%filename)
        return filename
        
    
        
