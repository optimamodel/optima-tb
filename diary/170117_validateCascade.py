"""
Cascade validation check

Required: 
    The following file is required for this test:
        1. 'tb-ucl\data\validate_cascade.xlsx'
          
@date:   17-Jan-2017
@author: hussainazfar

"""
from optima_tb.project import Project
import os

cascade =  os.path.abspath('../data/validate_cascade.xlsx')

#cascade validation is automatically run  when a project is created
proj= Project(name = 'validate-cascade', cascade_path = cascade)
#verify cascade through plot
proj.settings.plotCascade()