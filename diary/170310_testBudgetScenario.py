import pylab

from optima_tb.scenarios import ParameterScenario
from optima_tb.utils import odict
from optima_tb.settings import Settings, DO_NOT_SCALE
from optima_tb.project import Project
from optima_tb.plotting import plotScenarios


"""

Example usage for parameter scenario: 

This demo shows the underlying functionality of how to 
create and run a parameter scenario, and illustrates the different 
options

"""
options = dict()
options['progs_start'] = 2010.0
options['init_alloc'] = {'HT-DS':1e6,'SAT-DS':1e6,'HT-MDR':1e6}
options['saturate_with_default_budgets'] = False     # Do NOT set as True for scenario unless calibration and program impacts are well-reconciled.


budget_options = {'HT-DS': 333}
# budget_options = {'HT-DS': 4e6,'SAT-DS':0,'HT-MDR': 3e4}



def testProject():
    databook = 'data/databook-simple-cascade-autocalibration.xlsx'
    proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple-calibration.xlsx')
    proj.setYear([2000.,2030.],False) 
    proj.loadSpreadsheet(databook_path = databook)
    proj.makeParset(name = 'default_parset')
    proj.makeProgset(name = 'default_progset')
    
    scen_values = { 'test_scenario': {'type': 'Budget',
                                      'overwrite' : True, # it will overwrite scenario to the parset
                                      'run_scenario' : True,
                                      'scenario_values': budget_options}
                   }
    proj.createScenarios(scen_values)
    print proj.scenarios
    
    resultset = proj.runScenarios(original_parset_name = 'default_parset',
                                  original_progset_name='default_progset',
                                  original_budget_options=options,
                                  include_bau=False)

        
    return resultset, proj


testProject()
# results, proj = testProject()





# 
# scen_labels = ['BAU','Test1']
# characs = ['at_treat','vaccin','lt_inf','ac_inf','alive','test_act_prev']
# plotScenarios(results,scen_labels,proj.settings,proj.data,plot_charac=characs,plot_pops=['Pop1'])
# pylab.show()

