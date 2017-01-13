from optima_tb.project import Project
import pylab 
"""
This script runs a simulation for two parameter settings. 

Note: this script should be run from the home folder for this project i.e. ~/git/tb-ucl/


@date:   23-Nov-2016
@author: sjarvis

"""


databook = './training/Belarus-databook-cascade_161122.xlsx'
cascade = './training/cascade_161122.xlsx'


proj= Project(name = 'testParameters', cascade_path = cascade)
proj.setYear([2000.,2030.],False) 
proj.loadSpreadsheet(databook_path = databook)


# ----------------------------------------------------------------
# 1. Change the name of populations. These have to match the one that are in 
# Belarus-databook-cascade_161122.xlsx

pop1 = "0-4"
pop2 = "5-14"
pop3 = "15-64"
pop4 = "65+"
pop5 = "Prisoners"
pop6 = "15-64 HIV"
pop7 = "65+ HIV"



# ----------------------------------------------------------------
# 2. Change parameter values for parameters you're interested in. Note you have to 
# use the name that is in the cascade_161122.xlsx in the 'Parameter' sheet, in the 
#'Code Label' column

# For 1st population: 0-4
dict_change_params1 = {'v_rate': [0.02],
                      'xi_vac' : [0.02],
                      'betaMDR': [0.8]}
# For 2nd population: 5-14
dict_change_params2 = {'v_rate': [0.1],
                      'xi_vac' : [0.2],
                      'betaMDR': [0.85]}
# For 3rd population: 15-64
dict_change_params3 = {'v_rate': [0.02],
                      'xi_vac' : [0.2]}
# For 4th population: 65+
dict_change_params4 = {'v_rate': [0.25],
                      'xi_vac' : [0.2],
                      'betaMDR': [0.8]}
# For 5th population: Prisoners
dict_change_params5 = {'v_rate': [0.01],
                      'xi_vac' : [0.2],
                      'betaMDR': [0.8]}
# For 6th population: 15-64 HIV
dict_change_params6 = {'v_rate': [0.2],
                      'xi_vac' : [0.2],
                      'betaMDR': [0.8]}
# For 7th population: 65+ HIV
dict_change_params7 = {'v_rate': [0.2],
                      'xi_vac' : [0.2],
                      'betaMDR': [0.8]}


# ----------------------------------------------------------------
# 3. Set a name for the parameter we're testing. 
pname='newName' 


# ----------------------------------------------------------------
# 4. This is where the parameters get set and sent into the project. 
# The project is run and the results are plotted. 
rate_dict = {pop1 : dict_change_params1,
             pop2 : dict_change_params2,
             pop3 : dict_change_params3,
             pop4 : dict_change_params4,
             pop5 : dict_change_params5,
             pop6 : dict_change_params6,
             pop7 : dict_change_params7}


proj.makeManualCalibration(pname2,rate_dict)
proj.runSim(parset_name=pname2,plot=plot)
pylab.show()
