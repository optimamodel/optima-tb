import pylab
pylab.close('all')
from optima_tb.project import Project
from optima_tb.reconciliation import reconcileFunc
import os

#Initialise conditions
parset_name='default'
progset_name = 'default'
reconcile_for_year=2017
unitcost_sigma=0.05
attribute_sigma=0.20
impact_pars=None

#Setup project
cascade  = os.path.abspath('..//..//tb-ucl-analyses//belarus//Cascade Spreadsheets//cascade-belarus.xlsx')
databook = os.path.abspath('..//..//tb-ucl-analyses//belarus//Databook Spreadsheets//databook-belarus.xlsx')
proj= Project(name = 'Belarus-newmodelfit', cascade_path = cascade, validation_level = 'avert', plotting_level = 'dev')
proj.loadSpreadsheet(databook_path=databook)
proj.makeParset(name=parset_name)
proj.makeProgset(name=progset_name)

#test normal reconciliation
progset      = reconcileFunc(proj=proj, reconcile_for_year=reconcile_for_year, 
                             parset_name=parset_name, progset_name= progset_name,
                             unitcost_sigma=unitcost_sigma, attribute_sigma=attribute_sigma, 
                             impact_pars=None,orig_tvec_end=proj.settings.tvec_end)

#Setup new program budgets to test with compare outcomes and one incorrect prog_label to see if compare outcomes ignore that prog_label
budget_allocation = {'HF DS-TB': 0., 'AMB DS-TB': 1e12, 'wrong_prog_label': 1e6}
#Test reconcile with a new proposed budget
progset      = reconcileFunc(proj=proj, reconcile_for_year=reconcile_for_year,
                             parset_name=parset_name, progset_name= progset_name, 
                             unitcost_sigma=unitcost_sigma, attribute_sigma= attribute_sigma,
                             impact_pars = None, budget_allocation = budget_allocation, orig_tvec_end = proj.settings.tvec_end)
