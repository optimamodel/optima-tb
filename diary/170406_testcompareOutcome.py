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

#Set years for Simulation runs
proj.compareOutcomes(parset_name=parset_name, progset_name=progset_name, year=2017)
    