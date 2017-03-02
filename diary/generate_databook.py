from optima_tb.project import Project
import numpy as np
import os
 
#Setup Project and paths
cascade  = os.path.abspath('..//tests//cascade_spreadsheet//cascade_model_simple.xlsx')
databook = os.path.abspath('..//tests//databooks//databook_model_simple.xlsx')
proj= Project(name = 'Belarus-activetbfit', cascade_path = cascade)
proj.setYear([2000, 2030])

dt = proj.settings.tvec_dt
plot_over = [2000,2030]
proj.settings.plot_settings['x_ticks'] = [np.arange(plot_over[0],plot_over[1]+dt,5,dtype=int),np.arange(plot_over[0],plot_over[1]+dt,5,dtype=int)]
proj.settings.plot_settings['xlim'] = (plot_over[0]-0.5,plot_over[1]+0.5)

proj.loadSpreadsheet(databook_path=databook)
proj.makeParset()
results = proj.runSim(plot=True)
