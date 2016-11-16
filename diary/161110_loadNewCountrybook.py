import settings
from spreadsheet import load_spreadsheet

settings = settings.Settings('./diary/cascade-161025.xlsx')

## legacy: checking that we can load in the data from the country spreadsheet
data = load_spreadsheet(settings,'./diary/country-data.xlsx')