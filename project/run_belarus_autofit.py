# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 13:12:24 2016

@author: Lara
"""

import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from project import Project
import pylab

proj= Project(name = 'Belarus', cascade_path = './cascade-belarus.xlsx')
#proj.makeSpreadsheet(databook_path = './databook_belarus_template.xlsx', num_pops = 7, num_migrations = 3)



proj.loadSpreadsheet(databook_path = './databook-belarus.xlsx')
proj.makeParset()

# Run the autofit calibration process. ----------------------------
# The default amount of runs is 500 iterations. If you find it needs to run for more iterations, then 
# uncomment out the next line:
#proj.settings.autofit_params['MaxIter'] = 1000
# This line runs the autofit calibration process:
proj.runAutofitCalibration(new_parset_name='autofit')
# End of the autofit calibration process. -------------------------

# Run our simulation using the new process
results2 = proj.runSim(parset_name='autofit',plot=True)


pylab.show()