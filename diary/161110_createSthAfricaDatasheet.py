import settings
from spreadsheet import export_spreadsheet

"""

Aim: create the Sth Africa datasheet for NFH (WB)

16.11.10 Needs to be updated following discussion today with NFH and GAJ
16.11.28: required updates have been made. 

"""

settings = settings.Settings('./diary/cascade-161025.xlsx')
settings.tvec_observed_end = 2016.0

cb = settings.countrybook
# create notified cases
cb['sheet_names']['total_cases'] = 'Notified cases'
cb['sheet_names'].pop('incident_cases')

# Malnutrition, diabetes, HIV, alcohol use
cb['sheet_values']['comorbidity'].pop('HIV prevalence')
cb['sheet_values']['comorbidity']['Silicosis prevalence'] = ['populations']
cb['sheet_values']['comorbidity']['Diabetes prevalence'] = ['populations']
cb['sheet_values']['comorbidity']['Alcohol prevalence'] = ['populations']
cb['sheet_values']['comorbidity']['Malnutrition prevalence'] = ['populations']

cb['constants']['num_default_programs'] = 35


pop_names = ['0-4 years','5-14 years','15-64 years','65+ years',\
                              'Miners','Prisoners','Health Care Workers',\
                              '15-64 years HIV','65+ years HIV', 'Miners HIV',\
                              'Prisoners HIV','Health Care Workers HIV']

pop_names = ['0-4 years','5-14 years','15-64 years HIV-','15-64 years HIV+', \
             '65+ years HIV-','65+ years HIV+', \
             'Miners HIV-', 'Miners HIV+',\
             'Prisoners HIV-','Prisoners HIV+',\
             'Health Care Workers HIV-','Health Care Workers HIV+']

export_spreadsheet(settings,filename='./data-entry-SouthAfrica2.xlsx',num_pops=len(pop_names),
                   pop_names=pop_names)
