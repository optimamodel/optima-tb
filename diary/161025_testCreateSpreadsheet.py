import settings
from spreadsheet import export_spreadsheet

settings = settings.Settings('../data/cascade.xlsx')

export_spreadsheet(settings,num_pops=2)