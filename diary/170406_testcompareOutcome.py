import pylab
pylab.close('all')
from optima_tb.project import Project
import os

#Initialise conditions
parset_name='default'
progset_name = 'default'

#Setup project
cascade  = os.path.abspath('..//..//tb-ucl-analyses//belarus//Cascade Spreadsheets//cascade-belarus.xlsx')
databook = os.path.abspath('..//..//tb-ucl-analyses//belarus//Databook Spreadsheets//databook-belarus.xlsx')
proj= Project(name = 'Belarus-newmodelfit', cascade_path = cascade, validation_level = 'avert', plotting_level = 'dev')
proj.loadSpreadsheet(databook_path=databook)
proj.makeParset(name=parset_name)
proj.makeProgset(name=progset_name)

#Test compare Outcomes general functionality
proj.compareOutcomes(parset_name=parset_name, progset_name=progset_name, year=2017)

#Setup new program budgets to test with compare outcomes and one incorrect prog_label to see if compare outcomes ignore that prog_label
budget_allocation = {'HF DS-TB': 0., 'AMB DS-TB': 1e12, 'wrong_prog_label': 1e6}
#Test compare outcomes with a new proposed budget
proj.compareOutcomes(parset_name=parset_name, progset_name=progset_name, budget_allocation=budget_allocation, year=2017)