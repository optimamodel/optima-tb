
from project import Project
from spreadsheet import load_spreadsheet
import pprint # just during debugging
pp = pprint.PrettyPrinter(indent=4)



#proj.makeSpreadsheet(num_pops = 5)


proj= Project(name = 'test-Belarus', cascade_path = '../data/cascade-simple.xlsx')
proj.loadSpreadsheet(databook_path = '../diary/test-Belarus-data.xlsx')

proj.makeParset()
print proj.data
#r1, o1 = proj.runSim()


data = load_spreadsheet('../data/Belarus data entry sheet v4a.xlsx')
print data.keys()
pp.pprint(data)
print proj.settings