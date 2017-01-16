"""
Run model simulation and plot output:

Loads and runs model data

"""

from optima_tb.project import Project
import pylab

proj = Project(name = 'Belarus', cascade_path = './cascade-belarus.xlsx', validation_level = 'avert')
proj.settings.tvec_dt = 1.0

#set the year range we simulate over as starting in 1995:
proj.setYear([2000,2030],False)


proj.loadSpreadsheet(databook_path = './databook-belarus.xlsx')
proj.makeParset(name = 'default')

# run and plot simulations

results = proj.runSim(plot = True)
# calculate a score for how good a fit the model is to the observed data 
proj.calculateFit(results)

pylab.show()