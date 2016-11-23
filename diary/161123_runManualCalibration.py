
from project import Project
import numpy as np
import pylab 
from calibration import runRateCalibration

"""
Aim: test and validation that implementation of births

"""

num_pop = 2
plot = True


databook = './data/databook-simple-cascade-calibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple-calibration.xlsx')
proj.setYear([2000.,2030.],False) 

"""
# 1. make spreadsheet after implementing number format
proj.makeSpreadsheet(databook_path=databook, num_pops = num_pop)

"""
# 2, load spreadsheet with number format. 
# Test by setting aging --> n
proj.loadSpreadsheet(databook_path = databook)
 
 
#proj.makeParset()
#print proj.data
#r1, o1 = proj.runSim(plot=plot)

"""
pname = 'newparset'
proj.makeParset(name=pname)
ps = proj.parsets[pname]
# make changes
dict_change_params = {'v_rate': [0.2],
                      'rec_act': [0.8]}
pop_id = 0 # only change for first population
for k,v in dict_change_params.iteritems():
    par_index = ps.par_ids['cascade'][k]
    print "par_index = %g"%par_index
    ps.pars['cascade'][par_index].y[0] = np.array(v)
proj.runSim(parset_name=pname,plot=plot)
pylab.show()

"""
dict_change_params = {'v_rate': [0.02],
                      'birth_transit' : [5000],
                      'rec_act': [0.8]}
dict_change_params2 = {'v_rate': [0.2],
                      'birth_transit' : [500],
                      'rec_act': [0.08]}

rate_dict = {'Pop1' : dict_change_params,
             'Pop2' : dict_change_params2}


pname2='testCalib'
runRateCalibration(proj,pname2,rate_dict,plot=plot)

pylab.show()
