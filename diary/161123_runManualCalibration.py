
from optima.tb.project import Project
import numpy as np
import pylab 

"""
Aim: test and validation that implementation of births

"""

num_pop = 2
plot = True
#plot=False


databook = './data/databook-simple-cascade-calibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple-calibration.xlsx')
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
#print proj.data
results1 = proj.runSim(plot=plot)



proj.calculateFit(results1)


dict_change_params = {'v_rate': [0.02],
                      'birth_transit' : [5000],
                      'rec_act': [0.8]}
dict_change_params2 = {'v_rate': [0.2],
                      'birth_transit' : [500],
                      'rec_act': [0.08]}

rate_dict = {'Pop1' : dict_change_params,
             'Pop2' : dict_change_params2}


pname2='testCalib'
proj.makeManualCalibration(pname2,rate_dict)
results = proj.runSim(parset_name=pname2,plot=plot)
proj.calculateFit(results)
pylab.show()
