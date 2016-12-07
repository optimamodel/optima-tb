
from project import Project
from spreadsheet import load_old_spreadsheet

import pylab 

"""
Aim: test and validation that implementation of number format for input transitions 

"""

num_pop = 2
plot = False
validateNegNumbers = False


def runSimpleCascade(validation_lvl):
    
    databook = './data/databook-simple-cascadeNumber.xlsx'
     
    proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple.xlsx', validation_level = validation_lvl)
    proj.setYear([2000.,2031.],False) # simulate to 2031
    
    # 1. make spreadsheet after implementing number format
    #proj.makeSpreadsheet(databook_path=databook, num_pops = num_pop)
    
    
    # 2, load spreadsheet with number format
    proj.loadSpreadsheet(databook_path = databook)
     
     
    proj.makeParset()
    #print proj.data
    r1 = proj.runSim(plot=plot)
    



def runFullCascade(validation_lvl):
    # Check with full cascade
    databook2 = './data/databook-full-cascade.xlsx'
     
     
    proj2= Project(name = 'test-Belarus-full', cascade_path = 'data/cascade.xlsx', validation_level = validation_lvl) 
    proj2.setYear([2000.,2016.])
    
    # 1. make spreadsheet after implementing number format
    #proj2.makeSpreadsheet(databook_path=databook2, num_pops = num_pop)
     
    # 2, load spreadsheet with number format
    proj2.loadSpreadsheet(databook_path = databook2)
     
      
    proj2.makeParset()
    # print proj2.data
    ss = proj2.settings
    ps = proj2.parsets[0]
    
    
    r1 = proj2.runSim(plot=plot) 


for vlvl in ['ignore','avert','warn','error']:
    """ should only halt when vlvl = error """ 
    #plot=True
    runSimpleCascade(vlvl)

    pylab.show()




