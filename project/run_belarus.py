# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 13:12:24 2016

@author: Lara
"""

#new file
import sys
sys.path.append('../optima')

from project import Project

proj= Project(name = 'Belarus', cascade_path = './cascade_belarus.xlsx')

proj.makeSpreadsheet(databook_path = './databook_belarus_template.xlsx', num_pops = 7, num_migrations = 3)

proj.loadSpreadsheet(databook_path = './databook_belarus.xlsx')
proj.makeParset()

proj.runSim(plot = True)