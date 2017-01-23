import logging
logger = logging.getLogger(__name__)

from optima_tb.utils import odict, printv, OptimaException

import xlrd
from xlrd import open_workbook, XLRDError 
import xlsxwriter as xw
import numpy as np
from xlsxwriter.utility import xl_rowcol_to_cell as rc
from numpy import nan, isnan, array, logical_or, nonzero, shape 
import itertools



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

#%% Constants used across spreadsheet functions

DEFAULT_PATH = './country-data.xlsx'    # Default path if not called via a project method.

DATA_COL_TBSTRAINS      = 2 # Column index for location of TB strains
ASSUMPTION_COL_OFFSET   = 2 # after the number of years, where is the assumption column
BHL_INDICES             = {'best':0, 'low':1, 'high':2} # Define best-low-high indices
# Sheets have one of two structures: either listing only one value in the entire sheet, or many values
UNIVALUE_SHEETS         = ['Population size','TB prevalence','TB incidence','Total cases','Incident cases']
MULTIVALUE_SHEETS       = ['Cost and coverage','Other epidemiology']#,'Comorbidity','Cost and coverage','Testing and treatment','Constants',]

WB_COLORS = {'blue'         : '#00CCFF',
             'light_blue'   : '#B7EBFF',
             'grey'         : '#BFBFBF',
             'light_grey'   : '#D9D9D9',
             'green'        : '#69AD45',
             'white'        : '#ffffff'}

WB_FORMATS = {  'header_format': {  'bold': True,
                                    'align': 'center',
                                    'text_wrap': True,
                                    'valign': 'vcenter'},
                'right_format' : {'align': 'right',
                                  'bold': True,
                                  'valign': 'vcenter'},
                'center_format': {'align': 'center',
                                                   'bold': True,
                                                   'bg_color': WB_COLORS['light_blue'],
                                                   'border': 1,
                                                   'border_color': WB_COLORS['white'],
                                                   'valign': 'vcenter'},
                'lgrey_center_format': {'align': 'center',
                                                   'bold': True,
                                                   'bg_color': WB_COLORS['light_grey'],
                                                   'border': 1,
                                                   'border_color': WB_COLORS['white'],
                                                   'valign': 'vcenter'},
                'lblue_bground': {'align': 'center',
                                        'valign': 'vcenter',
                                        'bg_color': WB_COLORS['light_blue'],
                                        'border': 1,
                                        'border_color': WB_COLORS['white']
                                        },
                'lblue_bground_right': {'align': 'right',
                                        'valign': 'vcenter',
                                        'bg_color': WB_COLORS['light_blue'],
                                        'border': 1,
                                        'border_color': WB_COLORS['white']
                                        },
                'lgrey_bground_right': {'align': 'right',
                                        'valign': 'vcenter',
                                        'bg_color': WB_COLORS['light_grey'],
                                        'border': 1,
                                        'border_color': WB_COLORS['white']
                                        },
                'blue_bground': {'align': 'center',
                                        'valign': 'vcenter',
                                        'bg_color': WB_COLORS['blue'],
                                        'border': 1,
                                        'border_color': WB_COLORS['white']
                                        },
                'grey_bground': {'align': 'center',
                                        'valign': 'vcenter',
                                        'bg_color': WB_COLORS['grey'],
                                        'border': 1,
                                        'border_color': WB_COLORS['white']
                                        },
                'lgrey_bground': {'align': 'center',
                                        'valign': 'vcenter',
                                        'bg_color': WB_COLORS['light_grey'],
                                        'border': 1,
                                        'border_color': WB_COLORS['white']
                                        },
                'main_label': {'align': 'left',
                                                 'bg_color': WB_COLORS['grey'],
                                                 'font_size': 14,
                                                 'bold': True},
                'lgrey_label': {'align': 'left',
                                                 'bg_color': WB_COLORS['light_grey'],
                                                 'font_size': 12,
                                                 'bold': True},
                'year_header': {'align': 'center',
                                'text_wrap': True,
                                                 'bg_color': WB_COLORS['green'],
                                                 'bold': True}
              }
    



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
    
    
    

def load_old_spreadsheet(filename="simple.xlsx",verbose=2):
    """
    Load data entry spreadsheet, check the data and return as a data object.
    
    Note: data object structure is expected to change in near-future.
    
    """
    import os 
    
    filename = os.path.abspath(filename)
    logging.info('Loading data from %s...' % filename)
    
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
        logging.info('Loading "Populations"...')
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
        logging.info('Loading "%s"...' % sheetname)
        # this is for all sheets where only one value is being read in
        sheetdata = workbook.sheet_by_name(sheetname)
        thispar = sheets[sheetname][0]
        data[thispar] = read_univalue_sheet(sheetdata, n_years)
     
    for sheetname in MULTIVALUE_SHEETS: 
        logging.info('Loading "%s"...' % sheetname, 2, verbose)
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
            logging.debug("creating label for %s"%par)
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
    


def _writeYearHeader(ws,row,col_offset,start_year,end_year,formats):
    for i,year in enumerate(range(int(start_year),int(end_year)+1)):
        ws.write(row,col_offset+i,year,formats['year_header'])

def _writeSinglePopulationBlock(ws,row_index_start,col_index_start,col_index_end,pop_labels,formats,col_header_start=None,disagg=False,pop_label_index=0):
    """
    
    Params:
        ws                 worksheet
        row_index_start    row index of where block of data will be
        col_index_start    col index of where block of data will be
        col_index_end      ending col index of where block of data will be
        pop_labels         labels for row headers (
        formats
        col_header_start
        disagg
        pop_label_index    index of populations
        
    
    """
    
    index_total = []
    if col_header_start is None: col_header_start = col_index_start - 1
    final_row = 0
    row_index = 0

    for i,pop_label in enumerate(pop_labels):
        if not disagg:
            ws.write(row_index_start+row_index,col_header_start,"=%s"%pop_label,formats['right_format'])
        else:
            if '' in pop_label[:-1]:
                continue
            
            if '' == pop_label[-1]:
                row_index += 1
                continue
            
            logging.debug("ADDING %s"%pop_label)
            
            
            for (j,v) in enumerate(pop_label):
                cell_value = '%s'%v
                if j == pop_label_index:
                    cell_value = '='+cell_value
                    
                    
                if v=='Total': #@TODO make this use constants
                    index_total.append(row_index)
                    ws.write(row_index_start+row_index,j+col_header_start,cell_value,formats['lgrey_bground_right'])
                else:
                    ws.write(row_index_start+row_index,j+col_header_start,cell_value,formats['lblue_bground_right'])
                
        ws.write(row_index_start+row_index,col_index_end,'OR',formats['center_format'])
        __formatBlock(ws, row_index_start+row_index, row_index_start+row_index+1, col_index_start, col_index_end, formats['blue_bground'])
        __formatBlock(ws, row_index_start+row_index, row_index_start+row_index+1, col_index_end+1, col_index_end+2, formats['blue_bground'])
        final_row = row_index
        row_index += 1
        
    
    for j in index_total:
        __formatBlock(ws, row_index_start+j, row_index_start+j+1, col_index_start, col_index_end, formats['grey_bground']) # main block
        __formatBlock(ws, row_index_start+j, row_index_start+j+1, col_index_end+1, col_index_end+2, formats['grey_bground']) # assumption
        ws.write(row_index_start+j,col_index_end,'OR',formats['lgrey_center_format'])
    
    return final_row


def _createUnivalueSheet(ws,ws_name,cb_settings,formats,start_year,end_year,pop_labels):
    """
    
    """
    disaggs = cb_settings['sheet_classes']['univalue'][ws_name]
    row_index_start = cb_settings['constants']['row_index_start']   
    col_index_start = cb_settings['constants']['col_index_start']+len(disaggs) -1
    col_index_end = col_index_start + int(end_year-start_year) + 1
    
    #determine disaggregations, and make it so that output is consistent format for later use
    dis_list = __genDisaggLabels(disaggs,cb_settings)
    #create row headers
    _writeYearHeader(ws, row=1, col_offset=col_index_start, start_year=start_year, end_year=end_year,formats=formats)
    ws.write(1,col_index_end+1,'Assumption',formats['year_header'])
    #create col headers, and block
    rows_added = _writeSinglePopulationBlock(ws,row_index_start,col_index_start,col_index_end,dis_list,formats,col_header_start=col_index_start-len(disaggs),disagg=True,pop_label_index=disaggs.index('populations'))
    # additional prettification
    ws.set_column(0, 0, 15)
    ws.set_column(col_index_start,col_index_end, 10)
    ws.set_column(col_index_end+1, col_index_end+1, 15)
    ws.set_column(col_index_end, col_index_end, 5)
    
def __genDisaggLabels(disaggs,cb_settings):
    disagg_list = []
    for dag_type in disaggs:
        if cb_settings['constants'].has_key('total_%s'%dag_type):
            disagg_list.append([cb_settings['constants']['total_%s'%dag_type]]+cb_settings['disaggregations'][dag_type]+[""])
        else:
            disagg_list.append(cb_settings['disaggregations'][dag_type]+[""])
    
    if len(disagg_list)>0:
        tmp_list = [list(a) for a in itertools.product(*disagg_list)] 
        disagg_list = tmp_list 
    return disagg_list

def _create_multivalue_block(ws,ws_name,row_index,col_index_start,col_index_end,multivalue_labels,cb_settings,formats,main_label=True):
    for label,item in multivalue_labels.iteritems():
        
        # create label
        ws.write(row_index,0,label)
        if main_label:
            ws.set_row(row_index,None,formats['lgrey_label'])
        row_index += cb_settings['constants']['spacing_multivalue_label']
        
        if isinstance(item,odict):
            row_index = _create_multivalue_block(ws, ws_name,row_index, col_index_start,col_index_end, item, cb_settings, formats,main_label=False)
        elif isinstance(item,list):
            
            # fix up row labels based on whether we're disaggregating values or not
            dis_list = __genDisaggLabels(item,cb_settings)
            # write col headers (poplabels and disaggs)
            
            try: 
                pop_label_index = item.index('populations') 
            except:
                pop_label_index =-1
            rows_added = _writeSinglePopulationBlock(ws,row_index,col_index_start,col_index_end,dis_list,formats,col_header_start=col_index_start-len(item),disagg=True,pop_label_index=pop_label_index)
            # interproperty spacing 
            row_index += rows_added + cb_settings['constants']['spacing_interproperty']
            
    return row_index
            

def _create_multivalue_sheet(ws,ws_name,cb_settings,formats,start_year,end_year,pop_labels):
    """
    Create sheets with multiblock values
    
    Format is specified in cb_settings, where each subsection follows the format:
    
        cb_settings['sheet_values'][sheetname][subsection_label] = value
        
        where value can be:
            - <a list of disaggregation properties>                                  #indicating a main list
            or
            - an odict of format "subsection": <list of disaggregation properties>   #indicating a sublist
    
    Examples:
    
        cb_settings['sheet_values']['example1'] = odict()
        cb_settings['sheet_values']['example1']['HIV prevalence'] = ['populations','smears','strains']
        cb_settings['sheet_values']['example1']['Diabetes prevalence'] = ['populations','strains']
    
        cb_settings['sheet_values']['example2'] = odict()
        cb_settings['sheet_values']['example2']['Background testing rates'] = ['populations']
        cb_settings['sheet_values']['example2']['Testing for latent TB'] = odict()
        cb_settings['sheet_values']['example2']['Testing for latent TB']['Percentage of population tested for latent TB per year'] = ['populations']
        cb_settings['sheet_values']['example2']['Testing for latent TB']['Number of people initiating treatment for latent TB each year'] = ['populations']
        cb_settings['sheet_values']['example2']['Testing for active TB'] = odict()
        cb_settings['sheet_values']['example2']['Testing for active TB']['Number of people lost to follow up for active TB each year'] = ['regimens','populations']
        cb_settings['sheet_values']['example2']['Testing for active TB']['Number of people successfully completing treatment for active TB each year'] = ['regimens','populations']
             
    
    """
    row_index_start = cb_settings['constants']['row_index_start']   
    # col_index_start: + number of max disaggregations for subcategories
    col_index_start =  len(max(cb_settings['sheet_values'][ws_name].values(),key=len))
    col_index_end = col_index_start + int(end_year-start_year) + 1
    # headers
    _writeYearHeader(ws, row=1, col_offset=col_index_start, start_year=start_year, end_year=end_year,formats=formats)
    ws.write(1,col_index_end+1,'Assumption',formats['year_header'])
    # write the blocks for each of the values as specified in cb_settings['sheet_values'][ws_name]
    _create_multivalue_block(ws,ws_name,row_index_start,col_index_start,col_index_end,cb_settings['sheet_values'][ws_name],cb_settings,formats)

def _create_populations(ws,ws_name,num_populations,formats, headers):

    #create headers
    for i,header in enumerate(headers):
        ws.write(1,i+1,header,formats['header_format'])
    #create columns
    for n in range(num_populations):
        ws.write('A%s'%(n+3),'%s'%(n+1),formats['right_format'])
    # format and prettify
    ws.set_column(0, 0, 5)
    ws.set_column(1, 3, 15)
    __formatBlock(ws, 2, 2+num_populations, 1, 4, formats['lblue_bground'])
    
    poplabels = ['%s!$B$%g'%(ws_name,i+3) for i in range(num_populations)]
    return poplabels
        
def _create_other_epidemiology(ws,ws_name,cb_settings,formats,start_year,end_year,pop_labels):
    _create_multivalue_sheet(ws,ws_name,cb_settings,formats,start_year,end_year,pop_labels)

def _create_comorbidity(ws,ws_name,cb_settings,formats,start_year,end_year,pop_labels):
    _create_multivalue_sheet(ws,ws_name,cb_settings,formats,start_year,end_year,pop_labels)

def _create_testing_treatment(ws,ws_name,cb_settings,formats,start_year,end_year,pop_labels):
    _create_multivalue_sheet(ws,ws_name,cb_settings,formats,start_year,end_year,pop_labels)

def _create_programs(ws, ws_name, num_programs, formats, headers):
    #create headers
    
    for i,header in enumerate(headers):
        ws.write(1,i+1,header,formats['year_header'])
    #create columns
    for n in range(num_programs):
        ws.write('A%s'%(n+3),'%s'%(n+1),formats['right_format'])
    # format and prettify
    ws.set_column(0, 0, 5)
    ws.set_column(1,4, 40)
    ws.set_column(2,2,15)
    ws.set_column(3,3,25)
    ws.set_column(5,6,25)
    __formatBlock(ws, 2, 2+num_programs, 1, 1+len(headers), formats['lgrey_bground'])
    __formatBlock(ws, 2, 2+num_programs, 1, 2, formats['blue_bground'])

    proglabels = ['%s!$B$%g'%(ws_name,i+3) for i in range(num_programs)]
    return proglabels

def _create_cost_coverage(ws,ws_name,cb_settings,formats,start_year,end_year,pop_labels):
    pass

def _create_unitcost(ws,cb_settings):
    pass

def _create_poptransitions(ws,cb_settings):
    pass

def __formatBlock(ws,from_row,to_row,from_col,to_col,format):
    for i in range(from_row,to_row):
        for j in range(from_col,to_col):
            ws.write(i,j,None,format)
            
            
def createFormats(workbook):
    formats = odict()
    for fkey,fmt in WB_FORMATS.iteritems():
        formats[fkey] = workbook.add_format(fmt)
    
    return workbook, formats
        
def export_spreadsheet(settings,filename=DEFAULT_PATH,num_pops=4,verbose=2,pop_names=None):
    """
    
    
    """
    workbook = xw.Workbook(filename)
    workbook,formats = createFormats(workbook)
    # functions
    availfns = globals().copy()
    availfns.update(locals())
    #year range:
    start_year = settings.tvec_start
    end_year= settings.tvec_observed_end

    for local_name in settings.countrybook['sheet_names']:
        name = settings.countrybook['sheet_names'][local_name]
        logging.debug("%s %s"%(local_name, name))
        logging.info("Creating sheet for %s"%name)
        ws = workbook.add_worksheet(name)
        #label for each data sheet
        ws.write('A1', settings.countrybook['labels'][local_name])
        ws.set_row(0, None, formats['main_label'])
        #populate
        if local_name == 'populations':
            poplabels = _create_populations(ws, name, num_pops, formats, settings.countrybook['headers']['populations'])
            settings.countrybook['disaggregations']['populations'] = poplabels
            ## TODO: remove (as is just a debugging measure)
            if pop_names is not None:
                for i,pop in enumerate(pop_names):
                    ws.write(i+2,1,pop)
        
        elif local_name == 'programs':
            proglabels = _create_programs(ws, local_name, settings.countrybook['constants']['num_default_programs'], formats, settings.countrybook['headers']['programs'])
            settings.countrybook['disaggregations']['programs'] = proglabels
                    
        elif local_name in settings.countrybook['sheet_classes']['univalue']:
            _createUnivalueSheet(ws, local_name, settings.countrybook, formats, start_year,end_year,poplabels)
        else:
            createSheet = availfns.get('_create_%s'%local_name)
            if not createSheet:
                raise NotImplementedError("No method associated in creating sheet '%s'"%name)
            createSheet(ws,local_name,settings.countrybook,formats, start_year,end_year,poplabels)
        
def _load_populations(ws,name,cb,num_cols):
    return _load_univalue_sheet(ws, name, cb, num_cols,col_index=1)

def _load_programs(ws,name,cb,num_cols):   
    return _load_univalue_sheet(ws, name, cb, num_cols,col_index=1)


def _load_univalue_sheet(ws,ws_name,cb_settings,num_cols,col_index):
    disaggs = cb_settings['sheet_classes']['univalue'][ws_name]
    row_index_start = cb_settings['constants']['row_index_start']   
    col_index_start = cb_settings['constants']['col_index_start']+len(disaggs) -1
    col_index_end = col_index_start + int(end_year-start_year) + 1
    
    row = 0
    while row < sheetdata.nrows-1: 
        row += 1
        if sheetdata.cell_value(row,col_index_start-1) == '': # check whether label for row or whether there's a space between 
            continue # as this is a filler row

        assumptiondata = sheetdata.cell_value(row, col_index_end)
        if assumptiondata != '': 
            row_data = [assumptiondata] # Use if a non-blank assumption has been entered
        else:
            row_data = blank2nan(sheetdata.row_values(row, start_colx=dataColIndex, end_colx=dataColIndex+n_years)) # Data starts in 3rd column
        print row_data
    return {'tmp':'tmp'}
    
def load_spreadsheet(settings,filename=DEFAULT_PATH):
    import os 
    
    filename = os.path.abspath(filename)
    try: 
        ws = open_workbook(filename) # Open workbook
    except: 
        errormsg = 'Failed to load spreadsheet: file "%s" not found or other problem' % filename
        raise OptimaException(errormsg)
    
    # functions
    availfns = globals().copy()
    availfns.update(locals())
    #year range:
    start_year = settings.tvec_start
    end_year= settings.tvec_observed_end
    
    cb = settings.countrybook
    
    for local_name in settings.countrybook['sheet_names']:
        name = settings.countrybook['sheet_names'][local_name]
        logging.info("Loading sheet %s"%name)
        if local_name == 'populations':
            popdata = _load_populations(ws, name, cb, len(settings.countrybook['headers']['populations']))
            poplabels = popdata.keys()
            settings.countrybook['disaggregations']['populations'] = poplabels
        
        elif local_name == 'programs':
            progdata = _load_programs(ws, local_name, cb, len(settings.countrybook['headers']['programs']))
            proglabels = progdata.keys()
            settings.countrybook['disaggregations']['programs'] = proglabels
                    
        elif local_name in settings.countrybook['sheet_classes']['univalue']:
            pass
            #_createUnivalueSheet(ws, local_name, settings.countrybook, formats, start_year,end_year,poplabels)
        else:
            pass

    print settings.countrybook



