import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

import settings
from spreadsheet import export_spreadsheet

settings = settings.Settings('../data/cascade.xlsx')

export_spreadsheet(settings,num_pops=2,pop_names=['Pop1','Pep2'])

"""
# TODO 
debug and remove weird spaces within a property
formatting for populations
sheetnames = proper sheet names

then testing and treatment / economic

"""