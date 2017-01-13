"""
Set up of the Belarus project:

This file generates the databook spreadsheet from the cascade-belarus.xlsx, so that it can be
filled with values from the country datasheet. 

"""
from optima_tb.project import Project

proj = Project(name = 'Belarus', cascade_path = './cascade-belarus.xlsx', validation_level = 'error')
proj.makeSpreadsheet(databook_path = './databook-belarus-template.xlsx', num_pops = 6, num_migrations = 2)
