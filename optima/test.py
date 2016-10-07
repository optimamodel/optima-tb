#%% Imports

from utils import tic, toc
from project import Project



#%% Test project functionality

#tt = tic()
#
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
#
#ts2 = tic()
#p2 = Project(name = 'standard-test')
#toc(ts2, label = 'creating %s project' % p2.name)
#p2.settings.plotCascade()
#tsp2 = tic()
#p2.makeSpreadsheet(num_pops = 4)
#toc(tsp2, label = 'creating %s spreadsheet' % p2.name)
#tspl2 = tic()
#p2.loadSpreadsheet(databook_path = './standard-test-data-modded.xlsx')
#toc(tspl2, label = 'loading %s spreadsheet' % p2.name)
#tp2 = tic()
#p2.makeParset()
#toc(tp2, label = 'making parset for %s' % p2.name)
#r2, o2 = p2.runSim()
#
#toc(tt, label = 'entire process')


t_create = tic()
p = Project(name = 'standard-test')
toc(t_create, label = 'creating %s project' % p.name)
t_load = tic()
p.loadSpreadsheet(databook_path = './standard-test-data-transit-test.xlsx')
toc(t_load, label = 'loading %s spreadsheet' % p.name)
t_make = tic()
p.makeParset()
toc(t_make, label = 'making parset for %s' % p.name)
r, o = p.runSim()

invdt = int(1/p.settings.tvec_dt)
dlat = r[2].getComp('lte').popsize[invdt] / r[2].getComp('sus').popsize[0]
dvac = r[2].getComp('vac').popsize[invdt] / r[2].getComp('sus').popsize[0]
print dlat
print dvac

import scipy as sp
import scipy.optimize as spo

def func(f1,F1,f2,dt):
    return (1-sp.exp(-(f1+f2)*dt)*f1/(f1+f2))-F1

F1 = 0.5
f2 = 0.5
x = spo.fsolve(func, 1, args=(F1, f2, p.settings.tvec_dt), full_output = True)
print x

#import scipy as sp
#num_comps = 32
#num_links = 104
#num_pops = 3
#num_transfers = num_pops**2
#
#tvec_len = int((p.settings.tvec_end - p.settings.tvec_start)/p.settings.tvec_dt + 1)
#
#t_sparse = tic()
#
#list_source = []
#list_sink = []
#list_data = []
#for link in r[2].links:
#    if link.vals[0] != 0:
#        list_source.append(link.index_from)
#        list_sink.append(link.index_to)
#        list_data.append(link.vals[0])
#        list_source.append(link.index_from)
#        list_sink.append(link.index_from)
#        list_data.append(-link.vals[0])
#for comp in r[2].comps:
#    list_source.append(comp.index)
#    list_sink.append(comp.index)
#    list_data.append(1)
#M = sp.sparse.coo_matrix((list_data, (list_source, list_sink)), shape=(num_comps, num_comps))
#M_new_full = sp.linalg.fractional_matrix_power(M.toarray(),p.settings.tvec_dt)
#M_new = sp.sparse.coo_matrix(M_new_full)
#print M.row
#print M.col
#print M.data
#print M_new.row
#print M_new.col
#print M_new.data
#
#toc(t_sparse, label = 'testing sparse matrix powers for %s' % p.name)