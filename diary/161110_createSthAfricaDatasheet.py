import settings
from spreadsheet import export_spreadsheet

"""

Aim: create the Sth Africa datasheet for NFH (WB)

16.11.10 Needs to be updated following discussion today with NFH and GAJ
 

"""

settings = settings.Settings('./diary/cascade-161025.xlsx')
settings.tvec_observed_end = 2016.0

cb = settings.countrybook
# Malnutrition, diabetes, HIV, alcohol use
cb['sheet_values']['comorbidity']['Silicosis prevalence'] = ['populations','smears','strains']

export_spreadsheet(settings,filename='./data-entry-SthAfrica.xlsx',num_pops=6,
                   pop_names=['0-4 years','5-14 years','15-64 years','65+ years','Miners','Prisoners'])
