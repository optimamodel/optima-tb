"""
Databook validation check

Required: 
    The following files are required and placed at 'Homefolder'\data
        1. databook-simple-cascadeBirths.xlsx
        2. cascade-simple-births.xlsx
          
Generates: sample-output.xlsx, which will be the same format as databook-simple-cascade.xlsx


@date:   05-Jan-2017
@author: hussainazfar

"""
from project import Project
import os

num_pop = 2

databook = os.path.abspath('../tests/databooks/databook_model_simple.xlsx')
cascade =  os.path.abspath('../tests/cascade_spreadsheet/cascade_model_simple.xlsx')


proj= Project(name = 'validate-databook', cascade_path = cascade)
#validate databook is a function automatically run within the loadSpreadsheet
### Look at databook.databookValidation
proj.loadSpreadsheet(databook_path = databook)