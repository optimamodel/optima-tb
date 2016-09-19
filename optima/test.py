#%% Imports

from utils import tic, toc
from project import Project



#%% Test project functionality
tt = tic()
ts1 = tic()
p1 = Project(name = 'simple-test', cascade_path = './cascade-simple.xlsx')
toc(ts1, label = 'creating %s project' % p1.name)
r1 = p1.runSim()
p1.settings.plotCascade()
tsp1 = tic()
p1.makeSpreadsheet(num_pops = 2)
toc(tsp1, label = 'creating %s spreadsheet' % p1.name)
ts2 = tic()
p2 = Project(name = 'standard-test')
toc(ts2, label = 'creating %s project' % p2.name)
p2.runSim()
p2.settings.plotCascade()
tsp2 = tic()
p2.makeSpreadsheet(num_pops = 2)
toc(tsp2, label = 'creating %s spreadsheet' % p2.name)
toc(tt, label = 'entire process')