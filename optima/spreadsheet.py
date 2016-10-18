from utils import odict, printv, OptimaException

from numpy import nan, isnan, array, logical_or, nonzero, shape 
from xlrd import open_workbook, XLRDError 


"""

Functions for creation and loading of data entry sheets for Optima-TB

Note that as the data entry spreadsheets had already been created by an earlier
incarnation of the project, there is some legacy choices in spreadsheet design
which require some backwards compatibility magic to be performed.  

The code is based on Optima-HIV's loadspreadsheet.py and savespreadsheet.py combined,
and some functions are lifted directly from the code.

TODO:
    - checks and validation of loaded data
    - confirm data required out of spreadsheets
    - confirm data format for 'Testing and Treatment' tab
    - confirm whether we need Constants tab still 
    - export spreadsheet

@author:    sjarvis
            Sept 2016

"""


DATA_COL_TBSTRAINS      = 2 # Column index for location of TB strains
ASSUMPTION_COL_OFFSET   = 2 # after the number of years, where is the assumption column
BHL_INDICES             = {'best':0, 'low':1, 'high':2} # Define best-low-high indices
# Sheets have one of two structures: either listing only one value in the entire sheet, or many values
UNIVALUE_SHEETS         = ['Population size','TB prevalence','TB incidence','Total cases','Incident cases']
MULTIVALUE_SHEETS       = ['Cost and coverage','Other epidemiology']#,'Comorbidity','Cost and coverage','Testing and treatment','Constants',]


def forcebool(entry, location=''):
    """ Convert an entry to be Boolean """
    if entry in [1, 'TRUE', 'true', 'True', 't', 'T']:
        return 1
    elif entry in [0, 'FALSE', 'false', 'False', 'f', 'F']:
        return 0
    else:
        errormsg = 'Boolean data "%s" not understood in spreadsheet location "%s"' % (entry, location)
        raise OptimaException(errormsg)

def validate_data():
    #TODO implement
    pass

def blank2nan(thesedata):
    ''' Convert a blank entry to a nan '''
    return list(map(lambda val: nan if val=='' else val, thesedata))

def check_headings(sheetname,subparlist,parcount=0):
    # TODO make this function relevant again, or remove it
    try:
        thispar = subparlist[parcount] # Get the name of this parameter, e.g. 'popsize'
    except:
        errormsg = 'Incorrect number of headings found for sheet "%s"\n' % sheetname
        errormsg += 'Check that there is no extra text in the first two columns'
        raise OptimaException(errormsg)
    
    
    

def load_spreadsheet(filename="simple.xlsx",verbose=2):
    """
    Load data entry spreadsheet, check the data and return as a data object.
    
    Note: data object structure is expected to change in near-future.
    
    """
    printv('Loading data from %s...' % filename, 1, verbose)
    
    
    # Initial check: workbook can be loaded
    try: 
        workbook = open_workbook(filename) # Open workbook
    except: 
        errormsg = 'Failed to load spreadsheet: file "%s" not found or other problem' % filename
        raise OptimaException(errormsg)
    
       
    ###########################################################################################################
    ## Define the workbook and parameter names -- should match :
    ## -- export_spreadsheet in this module 
    ## /-- partable in parameters.py/ :TODO check and implement; potentially, move this section to parameters.py
    ###########################################################################################################

    sheets = odict()
    #sheets['Population'] = []
    sheets['Population size'] = ['popsize']
    sheets['Total cases'] = ['prevalence_case']
    sheets['TB prevalence'] = ['prevalence_rate']
    sheets['TB incidence'] = ['incidence_rate']
    sheets['Incident cases'] = ['incidence_case']
    sheets['Other epidemiology'] = ['all_deaths','tb_deaths','birth_rate']
    sheets['Cost and coverage'] = ['ccov']
    #sheets['Comorbidity'] = []
    #sheets['Programs'] = []
    #sheets['Testing and treatment'] = []
    #sheets['Unit cost'] = []
    #sheets['Constants'] = []


    ###########################################################################
    ## Setup data structures
    ###########################################################################
    

    ## Basic setup
    data = odict() # To hold data
    data['meta'] = odict()
    #data['meta']['date'] = today()
    data['meta']['sheets'] = sheets # Store parameter names
    
    ## Initialize populations
    # TODO: this isn't currently supplied in existing spreadsheets, so we'll have to enter it via aux methods for the moment
    data['pops'] = odict() # Initialize to empty list
    data['pops']['short'] = [] # Store short population/program names, e.g. "FSW"
    data['pops']['long'] = [] # Store long population/program names, e.g. "Female sex workers"
    # We'll leave in this for the moment, as it may make it easier when/if we set up age-based transitions
    data['pops']['age'] = [] # Store the age range for this population


    ## Initialize constants
    data['const'] = odict() # Initialize to empty list
    

    ###########################################################################
    ## Determine year range and populations
    ###########################################################################
    
    try:
        sheetdata = workbook.sheet_by_name('Populations')
        printv('Loading "Populations"...', 2, verbose)
        # TODO: implement for that happy day when we have a summary spreadsheet
    except XLRDError:
        ## -- Populations is not an explicit sheet for TB workbooks [... no comment]
        ## -- so we'll currently implement a workaround that we can infer the populations from 
        ## -- existing sheet 'Population size'
        sheetdata = workbook.sheet_by_name('Population size')
        thesedata = sheetdata.col_values(0) # data is in first column
        # take unique values and map to string
        thesedata = [str(p) for p in set(thesedata)]
        thesedata.remove('') # remove the empty value that has no doubt crept in ...
        data['pops']['short'] = thesedata
        data['pops']['long'] = thesedata
            
        data['years'] = [] # Initialize epidemiology data years
        for col in range(2,sheetdata.ncols): # from the 3rd col onwards
            thiscell = sheetdata.cell_value(0,col) # 0 is the 1st row which is where the year data should be
            if thiscell=='' and len(data['years'])>0: #  We've gotten to the end
                lastdatacol = col # Store this column number
                break # Quit
            elif thiscell != '': # Nope, more years, keep going
                data['years'].append(float(thiscell)) #TODO take years as a slice
        
        
    assumptioncol = lastdatacol + 1 # Figure out which column the assumptions are in; the "OR" space is in between
    n_years = len(data['years'])
        
    
    ###########################################################################
    ## Load data sheets
    ###########################################################################
        
    ## Loop over remaining sheets
    
    for sheetname in UNIVALUE_SHEETS:
        printv('Loading "%s"...' % sheetname, 2, verbose)
        # this is for all sheets where only one value is being read in
        sheetdata = workbook.sheet_by_name(sheetname)
        thispar = sheets[sheetname][0]
        data[thispar] = read_univalue_sheet(sheetdata, n_years)
     
    for sheetname in MULTIVALUE_SHEETS: 
        printv('Loading "%s"...' % sheetname, 2, verbose)
        sheetdata = workbook.sheet_by_name(sheetname)
        parameters = sheets[sheetname]
        if len(parameters)==1:
            # The unicorns: note that here we read in the programs data, and return it into 
            # the parameter of our data. Used in Cost Coverage, etc.
            readdata = read_multivalue_sheet(sheetdata, sheets[sheetname],n_years,isParameterized=False)
            data[parameters[0]] = readdata
        else:
            readdata = read_multivalue_sheet(sheetdata, sheets[sheetname],n_years)
            data.update(readdata)
        

    return data
        
        
           
def read_univalue_sheet(sheetdata, n_years):
    """
    Read in data for sheets that have only one parameter value for entire sheet, such as Population size, TB prevalence or Total cases
    
    Note that for some sheets, there will be disaggregation by TB strain. This is handled automatically by this function
    by determining the number of columns before the data itself begins (which is a somewhat fragile assumption).
    
    
    Returns:
        # for TBstrain data
        data = {'prevalence_rate' : {'0-4yrs' : {   'DS-TB' : [0.00,0.00, ...,0.09],
                                                    'MDR-TB': [0.00,0.00, ...,0.09],
                                                    'Total' : [0.00,0.00, ...,0.09],
                                                }
                                     '5-14yrs': {   'DS-TB' : [0.00,0.00, ...,0.09],
                                                    'MDR-TB': [0.00,0.00, ...,0.09],
                                                    'Total' : [0.00,0.00, ...,0.09],
                                                }
                }
        OR
        # for non-TBstrain data
        data = {'pop_size' : { '0-4yrs' : [45,76,800,...,1e6],
                                5-14yrs': [45,76,800,...,1e6],
                }
        
    
    """
    data = {}
    
    # Determine whether data starts in 2nd or 3rd column i.e. whether we have TB strain info or not
    if sheetdata.cell_value(0,2) != '': # the year data starts in 0,2 i.e. we DON'T have info about TB-strains
        dataColIndex = DATA_COL_TBSTRAINS
        dataTBStrain = None
    else:
        dataColIndex = DATA_COL_TBSTRAINS+1
        dataTBStrain = DATA_COL_TBSTRAINS
        
    lastColIndex = dataColIndex+n_years-1
    assumptioncol = lastColIndex + ASSUMPTION_COL_OFFSET 
    
    # Loop over rows in the workbook
    row = 0
    while row < sheetdata.nrows-1: 
        row += 1
        if sheetdata.cell_value(row,dataColIndex-1) == '': # check whether label for row or whether there's a space between 
            continue # as this is a filler row

        assumptiondata = sheetdata.cell_value(row, assumptioncol)
        if assumptiondata != '': 
            row_data = [assumptiondata] # Use if a non-blank assumption has been entered
        else:
            row_data = blank2nan(sheetdata.row_values(row, start_colx=dataColIndex, end_colx=dataColIndex+n_years)) # Data starts in 3rd column
        
        blh = str(sheetdata.cell_value(row, dataColIndex-1)).lower() # Read in whether indicator is best, low, or high
        pop = sheetdata.cell_value(row,0)
        if pop not in data.keys():
            data[pop] = {}

#         TODO: implement checks
#         if thispar=='tbprev': 
#             validatedata(thesedata, sheetname, thispar, row, checkblank=(blh=='best'), checkupper=True)  # Make sure at least the best estimate isn't blank
#         else: # popsize                
#             validatedata(thesedata, sheetname, thispar, row, checkblank=(blh=='best'))
        
        if dataTBStrain:
            tbstrain = str(sheetdata.cell_value(row,dataColIndex-2)).lower()
            if tbstrain not in data[pop].keys():
                data[pop][tbstrain] = {}
            data[pop][tbstrain][blh] = row_data
        else:
            data[pop][blh] = row_data 
            

    return data


def read_multivalue_sheet(sheetdata,parameters,n_years,isParameterized=True):
    """
    
    Params:
        sheetdata
        parameters        
        n_years            
        isParameterized   flag indicating whether individual parameter names will be specified or not
    
    
    Returns data in format:
        # if is parameterized, using the parameter names supplied
        data = { 'tb_deaths :  {'label': 'Percentage of TB-related deaths per year',
                                '0-4years' : [0.00,0.00, ...,0.09],
                                '5-14years': [0.00,0.00, ...,0.19],
                                '15+years' : [0.00,0.00, ...,0.14]
                                }
                 'all_deaths' :{'label': 'Percentage of people who die each year',
                                '0-4years' : [0.00,0.00, ...,0.09],
                                '5-14years': [0.00,0.00, ...,0.19],
                                '15+years' : [1.84,1.84, ...,2.14]
                                }
                }
        # OR if not parameterized, take the first column as the parameter name
        data = {'Active case finding among obligatory groups': {'Cost': [1.0],  
                                                                'Coverage': [200.],
                                                                'label': 'Active case finding among obligatory groups'},
                'Example treatment'                          : {'Cost': [10.0], 
                                                                'Coverage': [260.],
                                                                'label':'Example treatment'},
                }

    """
    data = {}
    # being overly-explicit here in case we have to generalize later ... 
    dataColIndex = DATA_COL_TBSTRAINS
    lastColIndex = dataColIndex+n_years-1
    assumptioncol = lastColIndex + ASSUMPTION_COL_OFFSET 

    row = 0
    parindex = 0 # what parameter we're up to ... 
    flag = False # indicates whether we're reading a parameter or not. Wouldn't need this if could read from multiple lines simultaneously :(
    while row < sheetdata.nrows-1: 
        row += 1

        if sheetdata.cell_value(row,0) != '':
            flag = True
            if isParameterized:
                par = parameters[parindex]
            else:
                par = str(sheetdata.cell_value(row,0))
            printv("creating label for %s"%par,2, 2)
            # TODO convert to odict
            data[par] = {'label': str(sheetdata.cell_value(row,0))} # for debugging purposes
            parindex += 1
            
        elif sheetdata.cell_value(row,dataColIndex-1) == '': # check whether label for row or whether there's a space between 
            flag = False
            continue # as this is a filler row

        # no longer an empty row? then we must have a new collection of data
        if flag: # can't have this as part of the same logic as above, as we might have a data continuing on same or next row
            key = str(sheetdata.cell_value(row,dataColIndex-1))
            if key is '':
                key = 'Total'
                
            assumptiondata = sheetdata.cell_value(row, assumptioncol)
            if assumptiondata != '': 
                row_data = [assumptiondata] # Use if a non-blank assumption has been entered
            else:
                row_data = blank2nan(sheetdata.row_values(row, start_colx=dataColIndex, end_colx=dataColIndex+n_years)) # Data starts in 3rd column
            
            # TODO: implement checks
            #validatedata(thesedata, sheetname, ....
            
            data[par][key] = row_data
                  
    return data
    
    
def export_spreadsheet(filename="simple.xlsx",verbose=2):
    """
    
    
    """
    # TODO implement
    pass


# Debugging and testing:
#data = load_spreadsheet('Belarus data entry sheet v4a.xlsx')
#print(data)