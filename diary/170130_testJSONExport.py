import json

from optima_tb.dataio import dumps, loads, importObj
from optima_tb.project import Project



# -------------------------------------------------
# 0. Setup
# -------------------------------------------------
num_pop = 2
plot = True
#plot=False
databook = 'data/databook-simple-cascade-autocalibration.xlsx'
proj= Project(name = 'test-dataio', cascade_path = 'data/cascade-simple-calibration.xlsx')
proj.loadSpreadsheet(databook_path = databook)
proj.makeParset() # make a default parset for population labels etc.



def dumpLoadCheck(obj):
    print("Trying to dump")
    x = dumps(obj)
    print("Trying to load")
    y = loads(x)
    #print obj 
    print y 
    return y
    
"""   
# -------------------------------------------------
# 1. Check export / import of each part of a project
# -------------------------------------------------
dumpLoadCheck(proj.parsets[0].pars[0][0].t)
print "A"
dumpLoadCheck(proj.parsets[0].pars[0][0])
print "B"
dumpLoadCheck(proj.parsets[0].pars[0])
print "C"
dumpLoadCheck(proj.parsets[0])
print "D"
dumpLoadCheck(proj.parsets)
print "E"

dumpLoadCheck(proj.data)
print "Data"


#dumpLoadCheck(proj.settings.parser)
print "Settings.parser"
dumpLoadCheck(proj.settings.autofit_params)
print "Settings.autofit params"
dumpLoadCheck(proj.settings.plot_settings)
print "Settings.plot settings"
dumpLoadCheck(proj.settings)
print "Settings"

dumpLoadCheck(proj.uid)

dumpLoadCheck(proj.results)

proj2 = dumpLoadCheck(proj)
print proj2.settings.autofit_params
print "F"
"""

# -------------------------------------------------
# 2. Check export functionality via project.py 
# -------------------------------------------------
filename = proj.exportProject()

# -------------------------------------------------
# 2a. Check import functionality and compare against original
# -------------------------------------------------

proj2 = importObj(filename)


print proj2.name, proj.name
#print proj2.parsets
#print proj.parsets



