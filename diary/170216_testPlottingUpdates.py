
from optima_tb.project import Project
from optima_tb.utils import odict

from optima_tb.plotting import plotCharacteristic

import pylab

path = '../tb-ucl-analyses/belarus/'
dt = 1.0/4
plot_over = (2000,2030)
proj = Project(name = 'Belarus', cascade_path = '%s/cascade-belarus.xlsx'%path, validation_level = 'avert', plotting_level = 'presentation')
proj.settings.tvec_dt = dt
#setup_belarusPlotting(proj,dt,plot_over)

#set the year range we simulate over as starting in 1995:
proj.setYear([1999,2030],False)


proj.loadSpreadsheet(databook_path = '%s/databook-belarus.xlsx'%path)
proj.makeParset(name = 'default')
settings = proj.settings
plotdict = settings.plot_settings
plotdict['ylim'] = 0
# run and plot simulations

#proj.parsets[0].pars['cascade']

results = proj.runSim()
"""
# plot all
plotCharacteristic(results, settings.charac_specs, proj.data, title='Test', outputIDs=['alive'], 
                       pop_labels = None, plot_total = False,
                       plot_observed_data=True, colors=None, plotdict=settings.plot_settings)
"""
# plot subset
plotCharacteristic(results, settings.charac_specs, proj.data, title='Test', outputIDs=['alive'], 
                       pop_labels = ['0-4','5-14'], plot_total = False,
                       plot_observed_data=True, colors=None, plotdict=plotdict)


# plot total of two populations
plotCharacteristic(results, settings.charac_specs, proj.data, title='Test', outputIDs=['alive'], 
                       pop_labels = ['0-4','5-14'], plot_total = True,
                       plot_observed_data=False, colors=None, plotdict=plotdict)

# plot total all populations

plotCharacteristic(results, settings.charac_specs, proj.data, title='Test', outputIDs=['alive'], 
                       pop_labels = None, plot_total = True,
                       plot_observed_data=False, colors=None, plotdict=plotdict)

pylab.show()