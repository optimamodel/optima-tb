from project import Project
import pylab 

"""
This script generates the country databook. It should be run after the default parameter
values have been filled in. 

Note: this script should be run from the home folder for this project i.e. ~/git/tb-ucl/


@date:   14-Nov-2016
@author: sjarvis

"""

"""
Populations are:
    - 0-4
    - 5-14
    - 15-64
    - 65+
    - Prisoners
    - 15-64 HIV
    - 65+ HIV
"""
num_pop = 7 

databook = './training/Belarus-databook-cascade_161122.xlsx'
cascade = './training/cascade_161122.xlsx'


proj= Project(name = 'test-full-cascade', cascade_path = cascade)
# Set the year range we wish to enter data points for: from 2000 to 2016 inclusive
proj.setYear([2000.,2016.])

# ----------------------------------------------------------------
# 1. make spreadsheet after implementing number format
proj.makeSpreadsheet(databook_path=databook, num_pops = num_pop)


# ----------------------------------------------------------------
# 2. open Belarus-databook-cascade_161122.xlsx and fill in population names 
# in 'Population Definitions' tab, with minimum and maximum ages. 
# Also, go to 'Transfer Definitions' and select 'y' for aging for relevant populations
# i.e. 0-4 --> 5-14 ; 5-14 --> 15-64, 15-64 --> 65+, and for 15-64 HIV --> 65+ HIV


