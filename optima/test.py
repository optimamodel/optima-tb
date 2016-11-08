#%% Imports

from utils import tic, toc
from project import Project



#%% Test project functionality

tt = tic()

#ts1 = tic()
#p1 = Project(name = 'simple-test', cascade_path = './cascade-simple.xlsx')
#toc(ts1, label = 'creating %s project' % p1.name)
#p1.settings.plotCascade()
#tsp1 = tic()
#p1.makeSpreadsheet(num_pops = 2)
#toc(tsp1, label = 'creating %s spreadsheet' % p1.name)
#tspl1 = tic()
#p1.loadSpreadsheet()
#toc(tspl1, label = 'loading %s spreadsheet' % p1.name)
#tp1 = tic()
#p1.makeParset()
#toc(tp1, label = 'making parset for %s' % p1.name)
#r1, o1 = p1.runSim()

ts2 = tic()
p2 = Project(name = 'standard-test')
toc(ts2, label = 'creating %s project' % p2.name)
p2.settings.plotCascade()
tsp2 = tic()
p2.makeSpreadsheet(num_pops = 4)
toc(tsp2, label = 'creating %s spreadsheet' % p2.name)
tspl2 = tic()
p2.loadSpreadsheet()#databook_path = './standard-test-data-modded.xlsx')
toc(tspl2, label = 'loading %s spreadsheet' % p2.name)
tp2 = tic()
p2.makeParset()
toc(tp2, label = 'making parset for %s' % p2.name)
r2, o2 = p2.runSim(plot = True)

toc(tt, label = 'entire process')


#t_create = tic()
#p = Project(name = 'standard-test')
#toc(t_create, label = 'creating %s project' % p.name)
#t_load = tic()
#p.loadSpreadsheet(databook_path = './standard-test-data-transit-test.xlsx')
#toc(t_load, label = 'loading %s spreadsheet' % p.name)
#t_make = tic()
#p.makeParset()
#toc(t_make, label = 'making parset for %s' % p.name)
#r, o = p.runSim()
#
#invdt = int(1/p.settings.tvec_dt)
#dpop = 1 - o['alive'][0][invdt]/o['alive'][0][0]
#dlat = r[2].getComp('lte').popsize[invdt] / r[2].getComp('sus').popsize[0]
#dvac = r[2].getComp('vac').popsize[invdt] / r[2].getComp('sus').popsize[0]
#print('Intended aging fraction per year: %f' % p.data['transfers']['aging'][0][0]['y'])
#print('Fraction aged after a year: %f' % dpop)
#print('Infected fraction after a year: %f' % dlat)
#print('Vaccinated fraction after a year: %f' % dvac)