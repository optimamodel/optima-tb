import pylab
pylab.close('all')
from optima_tb.project import Project
from optima_tb.reconciliation import reconcileFunc
from optima_tb.utils import odict
import os

#Initialise conditions
parset_name='default'
progset_name = 'default'
reconcile_for_year=2017
unitcost_sigma=0.0#5
attribute_sigma=0.#20
budget_sigma = 1.
impact_pars=None

#Setup project
cascade  = os.path.abspath('..//..//tb-ucl-analyses//belarus//Cascade Spreadsheets//cascade-belarus.xlsx')
databook = os.path.abspath('..//..//tb-ucl-analyses//belarus//Databook Spreadsheets//databook-belarus.xlsx')
proj= Project(name = 'Belarus-newmodelfit', cascade_path = cascade, validation_level = 'avert', plotting_level = 'dev')
proj.loadSpreadsheet(databook_path=databook)
proj.makeParset(name=parset_name)
proj.makeProgset(name=progset_name)

#Save pre-reconciled budget values of programs
values = odict()
values['original'] = odict()
for prog in proj.progsets[progset_name].prog_ids:
    index = proj.progsets[progset_name].prog_ids[prog]
    if prog not in values['original'].keys(): values['original'][prog] = odict()
    values['original'][prog]['budget'] = proj.progsets[progset_name].progs[index].getDefaultBudget(year=reconcile_for_year)

#test normal reconciliation
progset, output      = reconcileFunc(proj=proj, reconcile_for_year=reconcile_for_year, 
                             parset_name=parset_name, progset_name= progset_name,
                             unitcost_sigma=unitcost_sigma, attribute_sigma=attribute_sigma, budget_sigma = budget_sigma,
                             impact_pars=None, constrain_budget = False, orig_tvec_end=proj.settings.tvec_end)

#Save post-reconciled budget values of programs
values['reconciled'] = odict()
for prog in proj.progsets[progset_name].prog_ids:
    index = proj.progsets[progset_name].prog_ids[prog]
    if prog not in values['reconciled'].keys(): values['reconciled'][prog] = odict()
    values['reconciled'][prog]['budget'] = proj.progsets[progset_name].progs[index].getDefaultBudget(year=reconcile_for_year)
    print('Program Name: %s\n\tOriginal Budget: %s, Reconciled Budget: %s' %(prog, values['original'][prog]['budget'], values['reconciled'][prog]['budget']))

##Setup new program budgets to test with compare outcomes and one incorrect prog_label to see if compare outcomes ignore that prog_label
#budget_allocation = {'HF DS-TB': 0., 'AMB DS-TB': 1e12, 'wrong_prog_label': 1e6}
##Test reconcile with a new proposed budget
#progset      = reconcileFunc(proj=proj, reconcile_for_year=reconcile_for_year,
#                             parset_name=parset_name, progset_name= progset_name, 
#                             unitcost_sigma=unitcost_sigma, attribute_sigma= attribute_sigma, budget_sigma = budget_sigma,
#                             impact_pars = None, budget_allocation = budget_allocation, orig_tvec_end = proj.settings.tvec_end)
