# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 13:12:24 2016

@author: Lara
"""

import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from project import Project
import pylab

proj= Project(name = 'Belarus', cascade_path = './cascade_belarus.xlsx')

proj.makeSpreadsheet(databook_path = './databook_belarus_template.xlsx', num_pops = 7, num_migrations = 3)

proj.loadSpreadsheet(databook_path = './databook_belarus.xlsx')
proj.makeParset()

r1,o1,s1,results = proj.runSim(plot = True)
#proj.calculateFit(results)
pylab.show()