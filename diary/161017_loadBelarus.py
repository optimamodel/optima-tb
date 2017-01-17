import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from optima_tb.project import Project
import pylab 
import pprint # just during debugging
pp = pprint.PrettyPrinter(indent=4)


proj= Project(name = 'test-Belarus', cascade_path = './cascade-161025.xlsx')
##when we have to make the initial spreadsheet
#proj.makeSpreadsheet(num_pops = 5)

# take test-Belarus-data.xlsx, set values and save as Belarus-cascade-data.xlsx (so that no danger of overwriting)
proj.loadSpreadsheet(databook_path = './Belarus-cascade-data.xlsx')

proj.makeParset()
print proj.data
r1 = proj.runSim(plot=True)

pylab.show()


## legacy: checking that we can load in the data from the country spreadsheet
#data = load_spreadsheet('../data/Belarus data entry sheet v4a.xlsx')
#print data.keys()
#pp.pprint(data)
#print proj.settings
