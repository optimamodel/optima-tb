from optima_tb.project import Project
from optima_tb.calibration import performSensitivityAnalysis
import numpy as np
import pylab 

"""
Aim: Test and validate Sensitivity analysis
"""

num_pop = 2

databook = '../data/databook-simple-cascade-autocalibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = '../data/cascade-simple-calibration.xlsx')
proj.setYear([2000.,2030.]) 
proj.loadSpreadsheet(databook_path = databook)
proj.makeParset()
results1 = proj.runSim(plot=False)

performSensitivityAnalysis(proj = proj)
#proj.calculateFit(results1)
#
#
#dict_change_params = {'v_rate': [0.02],
#                      'birth_transit' : [5000],
#                      'rec_act': [0.8]}
#dict_change_params2 = {'v_rate': [0.2],
#                      'birth_transit' : [500],
#                      'rec_act': [0.08]}
#
#rate_dict = {'Pop1' : dict_change_params,
#             'Pop2' : dict_change_params2}
#
#
#pname2='testCalib'
#proj.makeManualCalibration(pname2,rate_dict)
#print len(proj.parsets)
#results = proj.runSim(parset_name=pname2,plot=plot)
#proj.calculateFit(results)
#pylab.show()