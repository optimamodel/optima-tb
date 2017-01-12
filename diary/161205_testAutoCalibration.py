import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from optima.tb.project import Project
import numpy as np
import pylab 

"""
Aim: test autocalibration

Note: the target for this autocalibration are the values contained in 

./data/databook-simple-cascade-autocalibration-answers.xlsx

which was used for initally choosing non-identical cascade values, running the model
and taking output datapoints (via visual inspection) as characteristic datapoints.

These characteristic datapoints were then entered into databook-simple-cascade-autocalibration.xlsx, 
setting the initial starting values for the compartments to the nearest significant figure and similarly
resetting the cascade values to their default values. Therefore, the characteristic values can be achieved
by a model, as described in *-answers.xlsx

"""

num_pop = 2

databook = './data/databook-simple-cascade-autocalibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = './data/cascade-simple-calibration.xlsx')
proj.setYear([2000.,2030.],False) 

"""
# 1. make spreadsheet after implementing number format
proj.makeSpreadsheet(databook_path=databook, num_pops = num_pop)

"""
# 2, load spreadsheet with number format. 
# Test by setting aging --> n
#print proj.settings.linkpar_specs
proj.loadSpreadsheet(databook_path = databook)
proj.makeParset()

# Examine the results of the uncalibrated model:
#results = proj.runSim(plot=True)

# --------------------------
# UNCOMMENT the following sections as necessary
# --------------------------


"""
# Autofit over *all* characteristics:
proj.runAutofitCalibration(new_parset_name='bob')
results1 = proj.runSim(parset_name='bob',plot=True)
"""


"""
# Autofit over a subset of characteristics:
proj.runAutofitCalibration(new_parset_name='charlie',target_characs=['vaccin','at_treat','lt_treat','lt_inf','ac_inf'])
results2 = proj.runSim(parset_name='charlie',plot=True)
"""


"""
# Autofit (bootstrapping) over a subset of characteristics:
proj.runAutofitCalibration(new_parset_name='eve',target_characs=['vaccin','at_treat']) 
proj.runAutofitCalibration(new_parset_name='frank',old_parset_name='eve',target_characs=['lt_treat','lt_inf','ac_inf'])

# Run using the second autofitted parameters
results3 = proj.runSim(parset_name='frank',plot=True)
"""


pylab.show()
