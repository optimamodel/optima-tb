from optima_tb.project import Project
from optima_tb.utils import odict
from optima_tb.analysis import calculateCumulativeDerivatives

import pylab

path = '../tb-ucl-analyses/belarus/'
dt = 1.0/4
plot_over = (2000,2020)
proj = Project(name = 'Belarus', cascade_path = '%s/Cascade Spreadsheets/cascade-belarus.xlsx'%path, validation_level = 'avert', plotting_level = 'presentation')
proj.settings.tvec_dt = dt
#setup_belarusPlotting(proj,dt,plot_over)

#set the year range we simulate over as starting in 1995:
proj.setYear([1999,2020],False)


proj.loadSpreadsheet(databook_path = '%s/Databook Spreadsheets/databook-belarus.xlsx'%path)
proj.makeParset(name = 'default')
settings = proj.settings


results = proj.runSim()


print results.getValueAt('alive',2017)
print results.getValueAt('spdu',2017)
print results.getValueAt('nddiag',2017)

