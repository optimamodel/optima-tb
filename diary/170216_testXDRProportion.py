
from optima_tb.project import Project
import numpy as np
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


dict_change_params = {}
rate_dict = {'15-64' : dict_change_params}


pname2='testCalib'
proj.makeManualCalibration(pname2,rate_dict)
results = proj.runSim(parset_name=pname2,plot=plot)
proj.calculateFit(results)
pylab.show()
