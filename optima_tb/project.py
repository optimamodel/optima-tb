# %% Imports
import logging
import os.path

# NOTE - To configure logging in individual scripts, use the commands below to reset the
# logger settings
import logging.config
logging_conf = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {
            'format': '%(asctime)-20s %(levelname)-8s %(message)s',
            'datefmt': '%d-%m-%y %H:%M:%S'
        },
    },
    'handlers': {
        'default': {
            'level': 'DEBUG',
            'formatter': 'standard',
            'class': 'logging.StreamHandler',
        },
    },
    'loggers': {
        '': {
            'handlers': ['default'],
            'level': 'DEBUG',
        },
    }
}
logging.config.dictConfig(logging_conf)


logger = logging.getLogger()

from optima_tb.utils import tic, toc, odict, OptimaException
from optima_tb.model import runModel
from optima_tb.settings import Settings
from optima_tb.parameters import ParameterSet, export_paramset, load_paramset
from optima_tb.programs import ProgramSet
from optima_tb.databook import makeSpreadsheetFunc, loadSpreadsheetFunc
from optima_tb.optimization import optimizeFunc, parallelOptimizeFunc
from optima_tb.calibration import makeManualCalibration, calculateFitFunc, performAutofit
from optima_tb.scenarios import ParameterScenario, BudgetScenario, CoverageScenario
from optima_tb.reconciliation import reconcile

from uuid import uuid4 as uuid
import numpy as np
from copy import deepcopy as dcp


# %% Project class (i.e. one self-contained geographical unit)

class Project(object):
    ''' The main Optima project class. Almost all Optima functionality is provided by this class. '''

    def __init__(self, name='default', cascade_path='../data/cascade.xlsx', **args):
        ''' Initialize project. '''

        self.name = name
        self.uid = uuid()

        self.settings = Settings(cascade_path=cascade_path, **args)
        self.data = odict()

        self.parsets = odict()
        self.progsets = odict()
        self.results = odict()

        self.scenarios = odict()

        logger.info("Created project: %s" % self.name)

    def resetParsets(self):
        ''' Convenience function called externally to delete all parsets. '''
        self.parsets = odict()

    def setYear(self, yearRange, observed_data=True):
        '''
        
        @param yearRange: tuple or list 
        @param observed_data: bool indicating whether to change observed date end or simulation date end 
        '''
        self.settings.tvec_start = yearRange[0]
        if observed_data:
            self.settings.tvec_observed_end = yearRange[1]
        else:
            self.settings.tvec_end = yearRange[1]

    def runSim(self, parset=None, parset_name='default', progset=None, progset_name=None, options=None, plot=False, debug=False, store_results=True, result_type=None, result_name=None):
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

        # results = runModel(settings = self.settings, parset = parset)
        results = runModel(settings=self.settings, parset=parset, progset=progset, options=options)

        toc(tm, label='running %s model' % self.name)

        if plot:
            tp = tic()
            self.plotResults(results=results, debug=debug)
            toc(tp, label='plotting %s' % self.name)

        if store_results:
            if result_name is None:
                result_name = 'parset_' + parset.name
                if not progset is None:
                    result_name = result_name + '_progset_' + progset.name
                if result_type is not None:
                    result_name = result_type + '_' + result_name
                k = 1
                while k > 0:
                    result_name_attempt = result_name + '_' + str(k)
                    k = k + 1
                    if result_name_attempt not in self.results.keys():
                        result_name = result_name_attempt
                        k = 0
            self.results[result_name] = results

        return results


    def optimize(self, parset=None, parset_name='default', progset=None, progset_name='default', options=None, max_iter=None, max_time=None):
        ''' Optimize model using a selected parset and store/return results. '''

        optimization_params = dcp(self.settings.optimization_params)

        if max_time is not None:    # Update asd settings with new values ...
            optimization_params['maxtime'] = max_time
        if max_iter is not None:
            optimization_params['maxiters'] = max_iter

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

        results = optimizeFunc(settings=self.settings, parset=parset, progset=progset, options=options, **optimization_params)

        return results


    def parallelOptimize(self, parset=None, parset_name='default', progset=None, progset_name='default', options=None,
                         num_threads=4, block_iter=10, max_blocks=10, max_iter=None, doplot=False, fullfval=False, randseed=None, max_time=None):
        ''' Like optimize, but parallel '''
        parallel_optimization_params = dcp(self.settings.parallel_optimization_params)

        if max_time is not None:    # Update asd settings with new values ...
            parallel_optimization_params['maxtime'] = max_time
        if max_iter is not None:
            parallel_optimization_params['maxiters'] = max_iter
        if num_threads is not None:
            parallel_optimization_params['num_threads'] = num_threads
        if block_iter is not None:
            parallel_optimization_params['block_iter'] = block_iter
        if max_blocks is not None:
            parallel_optimization_params['max_blocks'] = max_blocks
        if fullfval is not None:
            parallel_optimization_params['fullfval'] = fullfval

        if parallel_optimization_params['num_threads'] > 1:
            logging.warning("WARNING: Setting reltol to None during parallel optimization")
            parallel_optimization_params['reltol'] = None

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


        results = parallelOptimizeFunc(settings=self.settings, parset=parset, progset=progset, options=options,
                                       randseed=randseed, doplot=doplot, **parallel_optimization_params)

        return results


    def makeSpreadsheet(self, databook_path=None, num_pops=5, num_migrations=2, num_progs=0):
        ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''

        if databook_path is None: databook_path = '../data/' + self.name + '-data.xlsx'
        logger.info("Attempting to create databook %s" % databook_path)

        makeSpreadsheetFunc(settings=self.settings, databook_path=databook_path, num_pops=num_pops, num_migrations=num_migrations, num_progs=num_progs)


    def loadSpreadsheet(self, databook_path=None):
        ''' Load data spreadsheet into Project data dictionary. '''
        if databook_path is None: databook_path = '../data/' + self.name + '-data.xlsx'
        logger.info("Attempting to load databook %s" % databook_path)

        self.data = loadSpreadsheetFunc(settings=self.settings, databook_path=databook_path)



    def makeParset(self, name='default'):
        ''' Transform project data into a set of parameters that can be used in model simulations. '''

        if not self.data: raise OptimaException('ERROR: No data exists for project "%s".' % self.name)
        self.parsets[name] = ParameterSet(name=name)
        self.parsets[name].makePars(self.data)
        return self.parsets[name]

    def makeProgset(self, name='default'):
        ''' Transform project data into a set of programs that can be used in budget scenarios and optimisations. '''

        if not self.data: raise OptimaException('ERROR: No data exists for project "%s".' % self.name)
        self.progsets[name] = ProgramSet(name=name)
        self.progsets[name].makeProgs(data=self.data, settings=self.settings)
        return self.progsets[name]

    def reconcile(self, parset_name=None, progset=None, progset_name=None, reconcile_for_year=2017, sigma_dict=None, unitcost_sigma=0.05, attribute_sigma=0.20, budget_sigma=0.0, impact_pars=None, constrain_budget=True, overwrite=True, max_time=None, save_progset=True):
        '''Reconcile identified progset with identified parset such that impact parameters are as closely matched as possible
           Default behaviour is to overwrite existing progset
        '''
        # Make a copy of the original simulation end date
        orig_tvec_end = self.settings.tvec_end
        # Checks and settings for reconcile
        if parset_name is None:
            try:
                parset_name = self.parsets.keys()[0]
                logger.info('Parameter set was not identified for reconciliation, using parameter set: "%s"' % parset_name)
            except:
                raise OptimaException('No valid parameter sets exist within the project')

        # If a progset is provided, it is assumed temporary and a reference is implanted into the current project to be deleted after reconciliation.
        # TODO: Make Project.reconcile() and similar methods much more aligned in handling parsets/progsets by label or by reference.
        link_progset_temporarily = False
        temp_progset_name = ''
        if not progset is None:
            link_progset_temporarily = True
            k = -1
            while link_progset_temporarily:
                if k == -1: temp_progset_name = progset.name
                elif k == 0: temp_progset_name = 'temp'
                else: temp_progset_name = 'temp' + str(k)
                if not temp_progset_name in self.progsets:
                    self.progsets[temp_progset_name] = progset
                    progset_name = temp_progset_name
                    break
                k += 1

        if progset_name is None:
            try:
                progset_name = self.progsets.keys()[0]
                logger.info('Program set was not identified for reconciliation, using program set: "%s"' % progset_name)
            except:
                raise OptimaException('No valid program sets exist within the project')

        if not parset_name in self.parsets.keys(): raise OptimaException("ERROR: No parameter set '%s' found" % parset_name)
        if not progset_name in self.progsets.keys(): raise OptimaException("ERROR: No program set '%s' found" % progset_name)
        # If overwrite selected, reconcile will overwrite the progset, otherwise a new progset is created
        if not overwrite:
            progset_name += '_reconciled'
            self.makeProgset(name=progset_name)
        logger.info('Reconciling progset "%s" as overwrite is set as "%s"' % (progset_name, overwrite))

        # Run reconcile functionality
        reconciled_progset, reconciled_output = reconcile(proj=self, reconcile_for_year=reconcile_for_year,
                                                                parset_name=parset_name, progset_name=progset_name, sigma_dict=sigma_dict,
                                                                unitcost_sigma=unitcost_sigma, budget_sigma=budget_sigma, attribute_sigma=attribute_sigma,
                                                                impact_pars=impact_pars, constrain_budget=constrain_budget, max_time=max_time)

        if save_progset:
            self.progsets[progset_name] = reconciled_progset

        # Deletes temporary progset reference.
        if link_progset_temporarily:
            del self.progsets[temp_progset_name]
            try: del self.progsets[progset_name]    # In case overwrite is False delete the 'reconcile' version of the temp progset. TODO: Rethink all the logic.
            except: pass

        return reconciled_progset, reconciled_output


    def compareOutcomes(self, parset_name=None, progset_name=None, budget_allocation=None, year=2017):
        '''Display how parameters for a progset and parset match up
        '''
        compareOutcomesFunc(proj=self, parset_name=parset_name, progset_name=progset_name,
                            budget_allocation=budget_allocation, year=year, compareoutcome=True)

    def exportParset(self, parset_name):
        ''' Exports parset to .csv file '''
        if not parset_name in self.parsets.keys():
            raise OptimaException("ERROR: no parameter set '%s' found" % parset_name)
        export_paramset(self.parsets[parset_name])

    def importParset(self, parset_filename, new_name=None):
        ''' Imports parameter set from .csv file '''
        paramset = load_paramset(parset_filename)
        if new_name is None:
            new_name = paramset.name
        else:
            paramset.name = new_name

        if not new_name in self.parsets.keys():
            logger.info("Imported new parameter set: %s" % new_name)
        else:
            logger.info("Imported and overwriting parameter set: %s" % new_name)
        self.parsets[new_name] = paramset


    def makeManualCalibration(self, parset_name, rate_dict):
        ''' Take in dictionary of updated values for rates that can be used for manual calibration
            Update values dict: {pop_name : {parameter : value} '''
        if not parset_name in self.parsets.keys():
            self.makeParset(name=parset_name)
        paramset = self.parsets[parset_name]
        logger.info("Updating parameter values in parset=%s" % (parset_name))

        makeManualCalibration(paramset, rate_dict)


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
            raise OptimaException('ERROR: No result is specified. Cannot calculate fit.')

        if self.data is None or len(self.data) == 0:
            raise OptimaException('ERROR: No data is specified. Cannot calculate fit.')

        datapoints, _, _ = results.getCharacteristicDatapoints()
        score = calculateFitFunc(datapoints, results.t_observed_data, self.data['characs'], metric)
        logger.info("Calculated scores for fit using %s: largest value=%.2f" % (metric, np.max(score)))
        return score


    def runAutofitCalibration(self, parset=None, new_parset_name=None, old_parset_name="default", target_characs=None, max_time=None, save_parset=True):
        """
        Runs the autofitting calibration routine, as according to the parameter settings in the 
        settings.autofit_params configuration.
        
        Params:
            new_parset_name    name to save the resulting autofit to
            old_parset_name    name of the parset to use as a base. Default value="default"
        """

        if parset is None:
            if not old_parset_name in self.parsets.keys():
                self.makeParset(name=old_parset_name)
            parset = self.parsets[old_parset_name]

        if new_parset_name is None:
            # TODO: check that autofit doesn't already exist; if so, add suffix
            new_parset_name = "autofit"

        logger.info("About to run autofit on parameters using parameter set = %s" % old_parset_name)

        if max_time is not None:    # Update autocalibration settings with new time limit...
            prev_max_time = self.settings.autofit_params['maxtime']
            self.settings.autofit_params['maxtime'] = max_time
        try: new_parset = performAutofit(self, parset, new_parset_name=new_parset_name, target_characs=target_characs, **self.settings.autofit_params)
        except Exception as e:
            raise OptimaException("ERROR: Autocalibration failed.")
        if max_time is not None:    # ...and revert.
            self.settings.autofit_params['maxtime'] = prev_max_time

        logger.info("Created new parameter set '%s' using autofit" % new_parset_name)
        if save_parset: self.parsets[new_parset_name] = new_parset

        return new_parset


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
                logger.warn("Attempting to add scenario '%s' to project %s, that already contains scenario of the same name. Will ignore." % (scenario_name, self.name))
                # TODO decide what to do if scenario with same name already exists. Update or ignore? SJ: prefer to ignore.

            if vals['type'].lower() == 'parameter':
                self.scenarios[scenario_name] = ParameterScenario(name=scenario_name, settings=self.settings, pop_labels=pop_labels, **vals)

            elif vals['type'].lower() == 'budget':
                self.scenarios[scenario_name] = BudgetScenario(name=scenario_name, pop_labels=pop_labels, **vals)

            elif vals['type'].lower() == 'coverage':
                self.scenarios[scenario_name] = CoverageScenario(name=scenario_name, pop_labels=pop_labels, **vals)
            else:
                raise NotImplementedError("ERROR: No corresponding Scenario type for scenario=%s" % scenario_name)

        logger.info("Successfully created scenarios")


    def runScenarios(self, original_parset_name, original_progset_name=None, original_budget_options=None,
                     run_scenario_names=None, bau_label="BAU", include_bau=False,
                     scenario_set_name=None, plot=False, save_results=False, store_results=True):
        """
        Runs scenarios that are contained in this project's collection of scenarios (i.e. self.scenarios). 
        For each scenario run, using original_parset_name, the results generated are saved and 
        returns as a dictionary of results. 
        
        Optional flag include_bau indicates whether to additionally run a "business-as-usual" scenario
        i.e. no external changes. 
        
        
        # TODO implement ability to specify list of scenario names that will run if valid and regardless of scenario.run_scenario
        
        
        Params:
            original_parset_name    name of parameterSet to be used
            run_scenario_names      list of names to be run. If None, the scenario's run_scenario status is used instead (default)
            scenario_set_name       string for name of returned resultset
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


        logger.info("Running scenarios: ")

        if include_bau:
#             results[bau_label] = self.runSim(parset_name=original_parset_name, progset=orig_progset, options=original_budget_options, plot=plot)
            results[bau_label] = self.runSim(parset_name=original_parset_name, plot=plot)
            logger.info("Completed scenario case: %s" % bau_label)


#         logger.info(run_scenario_names)
        for scen in self.scenarios.keys():
            scen_name = '%s' % self.scenarios[scen].name
            if (run_scenario_names is not None and self.scenarios[scen].name in run_scenario_names) or (run_scenario_names is None and self.scenarios[scen].run_scenario):
                progset, budget_options = self.scenarios[scen].getScenarioProgset(orig_progset, original_budget_options)
                logger.info("Starting scenario case: %s" % scen_name)
                results[scen_name] = self.runSim(parset=self.scenarios[scen].getScenarioParset(orig_parset), progset=progset, options=budget_options, parset_name=scen_name, plot=plot, store_results=False) # Don't store results here since stored later
                logger.info("Completed scenario case: %s" % scen_name)
                if scenario_set_name is None:
                    results[scen_name].name = '%s' % (scen_name)
                else:
                    results[scen_name].name = '%s_%s' % (scenario_set_name, scen_name)

                if store_results:
                    result_name = results[scen_name].name

                    k = 1
                    while k > 0:
                        result_name_attempt = result_name + '_' + str(k)
                        k = k + 1
                        if result_name_attempt not in self.results.keys():
                            result_name = result_name_attempt
                            k = 0
                    self.results[result_name] = dcp(results[scen_name])

                    if save_results:
                        results[scen_name].export()
                        export_paramset(self.scenarios[scen].getScenarioParset(orig_parset))

        return results




#    def exportProject(self, filename=None, format='json', compression='zlib'):
#        """
#
#        This currently saves everything within a project, including results.
#
#        Params:
#            filename      filename to save to. If none is supplied, value is set to "<project.name>.project"
#            format        string for supported format types (json)
#            compression   string for supported compression types (zlib)
#
#        Usage
#            project = Project(name="sample", cascade="cascade.xlsx")
#            project.exportProject()
#            # saves to "sample.project.Z"
#            project.exportProject(filename="special")
#            # saves to "special.Z"
#
#        """
#        if filename is None:
#            filename = "%s.project"%self.name
#
#        logger.info("Attempting to save file in format=%s"%format)
#        filename = exportObj(self,filename=filename,format=format,compression=compression)
#        logger.info("Saved to file: %s"%filename)
#        return filename



