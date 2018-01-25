# %% Imports

from optima_tb.utils import odict, OptimaException
import optima_tb.settings as project_settings

import logging
logger = logging.getLogger(__name__)

import xlrd
import xlsxwriter as xw
import numpy as np
from xlsxwriter.utility import xl_rowcol_to_cell as rc
from numbers import Number
from copy import deepcopy as dcp


WB_COLORS = {'blue'         : '#00CCFF',
             'light_blue'   : '#B7EBFF',
             'grey'         : '#BFBFBF',
             'light_grey'   : '#D9D9D9',
             'green'        : '#69AD45',
             'white'        : '#ffffff'}

WB_FORMATS = {
                'lblue_bground': {'align': 'right',
                                        'valign': 'vcenter',
                                        'bg_color': WB_COLORS['light_blue'],
                                        'border': 1,
                                        'border_color': WB_COLORS['white']
                                        },
              }



# %% Utility functions to generate sub-blocks of the project databook

def makeValueEntryArrayBlock(worksheet, at_row, at_col, num_arrays, tvec, cell_formats=None, assumption=0.0, assumption_overrides=None, data_formats=None, print_conditions=None, print_conditions_blank_preloading=True, no_header=False, only_assumption=False):
    '''
    Create a block where users choose data-entry format and enter values either as an assumption or time-dependent array.
    
    Args:
        worksheet                           -   The worksheet in which to produce the block.
        at_row                              -   The row at which the top-left corner of the block will be fixed.
        at_col                              -   The column at which the top-left corner of the block will be fixed.
        num_arrays                          -   The number of input arrays to generate.
        tvec                                -   The year range for which values are (optionally) to be entered.
        assumption                          -   The default value to read from the assumption column in case the file is not first opened in Excel before loading.
                                                Without opening file in Excel, custom equations are not calculated.
        assumption_overrides                -   A list of default values to use (or equations to calculate) in the assumption column, assuming the file is opened in Excel.
                                                List must be of num_arrays length.
        data_formats                        -   A list of formats that data can be entered as.
                                                This is not data validation, but will determine the form of calculations during model processing.
        print_conditions                    -   A list of Excel-string conditions that are used to test whether data-entry arrays should be shown.
                                                List must be of num_arrays length, but can include values of None to allow default printing behaviour for certain rows.
        print_conditions_blank_preloading   -   If True, pre-loaded Excel files have default format/assumption values set to ''. This can cause problems when auto-loading.
                                                If False, default format/assumption values are loadable, even if opening the Excel file then blanks them out.
        no_header                           -   A flag for whether the output excel block should contain assumption and year column headers.
        only_assumption                     -   Locks out value entry for anything but the assumption.
                                
    Note that if no data format is specified, data formats for 'Fraction' or 'Number' are available. The default choice of the 
    data format is assumed that if the assumption value is larger than 1., then it is likely to be a Number, otherwise it is assumed to be a Fraction.
    '''

    if data_formats is None:
        data_formats = ['Fraction', 'Number']  # , 'Probability', 'Number']
        data_format_assumption = data_formats[0]
        if assumption > 1.:
            data_format_assumption = data_formats[1]  # Number
    else:
        data_format_assumption = data_formats[0]

    # Handles issues where user wants to autoload assumptions rather print-condition blanks.
    if not print_conditions_blank_preloading:
        preloaded_data_format = data_format_assumption
        preloaded_assumption = assumption
    else:
        preloaded_data_format = ''
        preloaded_assumption = ''

    print preloaded_data_format

    offset = at_col + 3  # This is the column at which the input year range and corresponding values should be written.

    if not no_header:
        worksheet.write(at_row, at_col, 'Format')
        worksheet.write(at_row, at_col + 1, 'Assumption')
        for k in xrange(len(tvec)):
            worksheet.write(at_row, offset + k, tvec[k])
        header_offset = 1
    else:
        header_offset = 0

    # If no overrides are provided or they are of incorrect length, just revert to the original assumption.
    if assumption_overrides is None or len(assumption_overrides) != num_arrays:
        assumption_overrides = [assumption] * num_arrays

    value_default = ''
    if only_assumption: value_default = '...'

    for aid in xrange(num_arrays):
        row_id = at_row + aid + header_offset
        offset = at_col + 3

        worksheet.write(row_id, at_col, data_format_assumption)  # Default choice for data format.
        worksheet.data_validation('%s' % rc(row_id, at_col), {'validate': 'list', 'source': data_formats, 'ignore_blank': False})
        if not print_conditions is None and not print_conditions[aid] is None:
            worksheet.write(row_id, at_col, '=IF(%s,"%s","")' % (print_conditions[aid], data_format_assumption), cell_formats['lblue_bground'], preloaded_data_format)  # Default choice for data format.
            worksheet.write(row_id, at_col + 1, '=IF(%s,IF(SUMPRODUCT(--(%s:%s<>"%s"))=0,%s,"N.A."),"")' % (print_conditions[aid], rc(row_id, offset), rc(row_id, offset + len(tvec) - 1), value_default, assumption_overrides[aid]), None, preloaded_assumption)
            worksheet.write(row_id, at_col + 2, '=IF(%s,"OR","")' % print_conditions[aid], None, '')
            if only_assumption:
                for k in xrange(len(tvec)):
                    worksheet.write(row_id, at_col + 3 + k, '=IF(%s,"%s","")' % (print_conditions[aid], value_default), cell_formats['lblue_bground'], '')
                    worksheet.data_validation('%s' % rc(row_id, at_col + 3 + k), {'validate': 'list', 'source': ['%s' % value_default], 'ignore_blank': False})
        else:
            worksheet.write(row_id, at_col, data_format_assumption)  # Default choice for data format.
            worksheet.write(row_id, at_col + 1, '=IF(SUMPRODUCT(--(%s:%s<>"%s"))=0,%s,"N.A.")' % (rc(row_id, offset), rc(row_id, offset + len(tvec) - 1), value_default, assumption_overrides[aid]), cell_formats['lblue_bground'], assumption)
            worksheet.write(row_id, at_col + 2, 'OR')
            if only_assumption:
                for k in xrange(len(tvec)):
                    worksheet.write(row_id, at_col + 3 + k, '%s' % value_default, cell_formats['lblue_bground'])
                    worksheet.data_validation('%s' % rc(row_id, at_col + 3 + k), {'validate': 'list', 'source': ['%s' % value_default], 'ignore_blank': False})

        if not only_assumption: # for the cases that we do not have an assumption, but still need to format the cells for data entry
            for k in xrange(len(tvec)):
                worksheet.write(row_id, at_col + 3 + k, None, cell_formats['lblue_bground'])


def makeTagMatrix(worksheet, at_row, num_rows, at_col, labels, formula_labels=None, allowed_vals=None, no_validation=False):
    '''
    Create a matrix where users tag relations between one (left-column positioned) list and a potentially different (top-row positioned) list.
    The tags are forced to be from a finite-discrete allowed_vals set unless no_validation is turned on.
    Note that the left-column positioned list is not included in this matrix.
    
    Args:
        worksheet       -   The worksheet in which to produce the matrix.
        at_row          -   The row at which the top-left corner of the matrix will be fixed.
        num_rows        -   Excluding the matrix header row, the number of rows (i.e. taggable items) for this matrix.
        at_col          -   The column at which the top-left corner of the matrix will be fixed.
        labels          -   A list of labels for objects to be related to (e.g. populations for programs).
        formula_labels  -   An optional list corresponding to labels that contains Excel formulas.
                            The labels list is still required for default values.
                            In the case that the databook is not opened in Excel prior to loading, the defaults are read in.
        allowed_vals    -   A list of allowed values that matrix elements can be tagged by.
        no_validation   -   If true, does not enforce Excel validation for the allowed values.
    '''

    if allowed_vals is None: allowed_vals = ['n', 'y']
    if formula_labels is None: formula_labels = dcp(labels)
    if len(formula_labels) != len(labels): raise OptimaException('ERROR: The numbers of formula-based and default labels passed to a matrix-writing process do not match.')

    for k in xrange(len(labels)):
        worksheet.write(at_row, at_col + k, formula_labels[k], None, labels[k])

    for row_pre in xrange(num_rows):
        row_id = at_row + row_pre + 1
        for col_pre in xrange(len(labels)):
            col_id = at_col + col_pre

            worksheet.write(row_id, col_id, allowed_vals[0])  # Default choice for data format.
            if not no_validation:
                worksheet.data_validation('%s' % rc(row_id, col_id), {'validate': 'list', 'source': allowed_vals, 'ignore_blank': False})


def makeConnectionMatrix(worksheet, at_row, at_col, labels, formula_labels=None, allowed_vals=None, no_validation=False, symmetric=False, self_connections=''):
    '''
    Create a matrix where users tag connections from a (left-column positioned) list to the same (top-row positioned) list.
    The tags are forced to be from a finite-discrete allowed_vals set unless no_validation is turned on.
    Diagonals (self-connections) cannot be filled with values without giving self_connections a value.
    
    Args:
        worksheet       -   The worksheet in which to produce the matrix.
        at_row          -   The row at which the top-left corner of the matrix will be fixed.
        at_col          -   The column at which the top-left corner of the matrix will be fixed.
        labels          -   A list of labels for objects to be connected (e.g. population groups).
        formula_labels  -   An optional list corresponding to labels that contains Excel formulas.
                            The labels list is still required for default values.
                            In the case that the databook is not opened in Excel prior to loading, the defaults are read in.
        allowed_vals    -   A list of allowed values that matrix elements can be tagged by.
        no_validation   -   If true, does not enforce Excel validation for the allowed values.
        symmetric       -   If true, seeds half the matrix with Excel equations that equal the other half.
    '''

    if allowed_vals is None: allowed_vals = ['n', 'y']
    if formula_labels is None: formula_labels = dcp(labels)
    if len(formula_labels) != len(labels): raise OptimaException('ERROR: The numbers of formula-based and default labels passed to a matrix-writing process do not match.')

    for k in xrange(len(labels)):
        worksheet.write(at_row + k + 1, at_col, formula_labels[k], None, labels[k])
        worksheet.write(at_row, at_col + k + 1, formula_labels[k], None, labels[k])

    for row_pre in xrange(len(labels)):
        row_id = at_row + row_pre + 1
        for col_pre in xrange(len(labels)):
            col_id = at_col + col_pre + 1

            if row_pre != col_pre:
                worksheet.write(row_id, col_id, allowed_vals[0])  # Default choice for data format.
                if symmetric and row_pre > col_pre:
                    sym_rc = rc(at_row + col_pre + 1, at_col + row_pre + 1)
                    worksheet.write(row_id, col_id, '=IF(%s<>"",%s,"")' % (sym_rc, sym_rc), None, allowed_vals[0])
                if not no_validation:
                    worksheet.data_validation('%s' % rc(row_id, col_id), {'validate': 'list', 'source': allowed_vals, 'ignore_blank': False})
            else:
                worksheet.write(row_id, col_id, self_connections)
                if not no_validation:
                    worksheet.data_validation('%s' % rc(row_id, col_id), {'validate': 'list', 'source': [self_connections], 'ignore_blank': False})




# %% Function for generating project databook

# NOTE: Comment this better later, especially with the fact that all Excel formulae need a fall-back value.
def makeSpreadsheetFunc(settings, databook_path, num_pops=5, num_migrations=2, num_progs=0):
    ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''

    include_progs = False
    if num_progs > 0:
        if len(settings.progtype_specs.keys()) > 0:
            include_progs = True
        else:
            raise OptimaException('ERROR: Attempting to create a databook with program sections, despite no program types defined in loaded cascade.')

    workbook = xw.Workbook(databook_path)
    workbook, cell_formats = createFormats(workbook)

    ws_pops = workbook.add_worksheet(settings.databook['sheet_names']['pops'])
    ws_contact = workbook.add_worksheet(settings.databook['sheet_names']['contact'])
    ws_transmat = workbook.add_worksheet(settings.databook['sheet_names']['transmat'])
    ws_transval = workbook.add_worksheet(settings.databook['sheet_names']['transval'])
    if include_progs:
        ws_progmat = workbook.add_worksheet(settings.databook['sheet_names']['progmat'])
        ws_progval = workbook.add_worksheet(settings.databook['sheet_names']['progval'])

    # Regarding cascade parameters and characteristics, store sheets and corresponding row ids for further writing.
    ws_params = odict()
    for custom_label in settings.databook['custom_sheet_names']:
        ws_params[custom_label] = {'row_id':0, 'ws':workbook.add_worksheet(settings.databook['custom_sheet_names'][custom_label])}

    # Only produce standard sheets if there are any parameters left over that are not allocated to custom sheets.
    if settings.make_sheet_characs:
        ws_params['charac'] = {'row_id':0, 'ws':workbook.add_worksheet(settings.databook['sheet_names']['charac'])}
    if settings.make_sheet_linkpars:
        ws_params['linkpars'] = {'row_id':0, 'ws':workbook.add_worksheet(settings.databook['sheet_names']['linkpars'])}

    data_tvec = np.arange(settings.tvec_start, settings.tvec_observed_end + 1.0 / 2)
    ws_pops_width = 15
    ws_contact_width = 25
    ws_transmat_width = 15
    ws_transval_width = 15
    ws_progmat_width = 15
    ws_progval_width = 35
    name_width = 60
    assumption_width = 10


    # %% Population names sheet.
    age_min_col = 2
    age_max_col = 3
    ws_pops.write(0, 0, 'Name')
    ws_pops.write(0, 1, 'Abbreviation')
    ws_pops.write(0, age_min_col, 'Minimum Age')
    ws_pops.write(0, age_max_col, 'Maximum Age')
    ws_pops.write(0, 4, 'Disability Weight')
    ws_pops.write(0, 5, 'Average Life Expectancy')


    # While writing default population names and labels, they are stored for future reference as well.
    pop_names_default = []
    pop_labels_default = []
    for pid in xrange(num_pops):
        pop_name = 'Population ' + str(pid + 1)
        pop_label = pop_name[0:3] + str(pid + 1)
        pop_names_default.append(pop_name)
        pop_labels_default.append(pop_label)
        ws_pops.write(pid + 1, 0, pop_name)
        ws_pops.write(pid + 1, 1, '=LEFT(%s,3)&"%i"' % (rc(pid + 1, 0), pid + 1), None, pop_label)
    ws_pops.set_column(0, 3, ws_pops_width)

    # Excel formulae strings that point to population names and labels are likewise stored.
    pop_names_formula = []
    pop_labels_formula = []
    for pid in xrange(num_pops):
        pop_names_formula.append("='%s'!%s" % (settings.databook['sheet_names']['pops'], rc(pid + 1, 0, True, True)))
        pop_labels_formula.append("='%s'!%s" % (settings.databook['sheet_names']['pops'], rc(pid + 1, 1, True, True)))

    # %% Inter-population contact sheet to weight interactions between population groups

    ws_contact.write(0, 0, 'Interaction Impact Weights')
    makeConnectionMatrix(worksheet=ws_contact, at_row=0, at_col=0, labels=pop_labels_default, formula_labels=pop_labels_formula, allowed_vals=[''], no_validation=True, symmetric=True, self_connections=1)
    ws_contact.set_column(0, 0, ws_contact_width)

    # %% Inter-population transfers matrix sheet (from node to corresponding node).

    # Produce a matrix for each 'migration type' (e.g. prisoner-transfers, migration mixing, HIV infection, etc.).
    # Aging is the default first migration type and is always present.
    # Store type names and formulae strings that reference them. Also store the rows at which migration matrices start.
    row_offset = 0
    mig_types_default = []
    mig_types_formula = []
    mig_matrix_rows = []
    for mid in xrange(num_migrations + 1):
        if mid == 0:
            mig_type = 'Aging'
        else:
            mig_type = 'Migration Type ' + str(mid)
        mig_types_default.append(mig_type)
        mig_matrix_rows.append(row_offset)
        mig_types_formula.append("='%s'!%s" % (settings.databook['sheet_names']['transmat'], rc(row_offset, 0)))
        ws_transmat.write(row_offset, 0, mig_type)
        makeConnectionMatrix(worksheet=ws_transmat, at_row=row_offset, at_col=0, labels=pop_labels_default, formula_labels=pop_labels_formula)
        row_offset += num_pops + 2

    ws_transmat.set_column(0, 0, ws_transmat_width)

    # %% Inter-population transfer details sheet.
    sh_pops = settings.databook['sheet_names']['pops']  # Convenient abbreviation.

    row_id = 0
    for mid in xrange(num_migrations + 1):
        ws_transval.write(row_id, 0, mig_types_formula[mid], None, mig_types_default[mid])
        print_conditions = []
        assumption_overrides = []

        k = 0  # A counter for number of arrays, used to id print conditions.
        print_row = row_id
        for source_id in xrange(num_pops):
            for target_id in xrange(num_pops):
                if source_id != target_id:
                    row_id += 1
                    r = mig_matrix_rows[mid] + source_id + 1
                    c = target_id + 1
                    ws_transval.write(row_id, 0, "=IF('%s'!%s=%s,%s,%s)" % (settings.databook['sheet_names']['transmat'], rc(r, c), '"y"', pop_names_formula[source_id][1:], '"..."'), None, '...')
                    ws_transval.write(row_id, 1, "=IF('%s'!%s=%s,%s,%s)" % (settings.databook['sheet_names']['transmat'], rc(r, c), '"y"', '"--->"', '""'), None, '')
                    ws_transval.write(row_id, 2, "=IF('%s'!%s=%s,%s,%s)" % (settings.databook['sheet_names']['transmat'], rc(r, c), '"y"', pop_names_formula[target_id][1:], '""'), None, '')

                    print_conditions.append('%s<>"..."' % rc(print_row + k + 1, 0))
                    if mid == 0:  # Aging assumptions are equations that try to work out aging fraction based on Minimum and Maximum pop age.
                        assumption_equation = "IF(AND('%s'!%s<>%s,'%s'!%s<>%s),1/('%s'!%s-'%s'!%s+1),0)" % (sh_pops, rc(source_id + 1, age_max_col), '""', sh_pops, rc(source_id + 1, age_min_col), '""', sh_pops, rc(source_id + 1, age_max_col), sh_pops, rc(source_id + 1, age_min_col))
                        assumption_overrides.append(assumption_equation)
                    k += 1
        makeValueEntryArrayBlock(worksheet=ws_transval, at_row=print_row, at_col=3, num_arrays=num_pops * (num_pops - 1), tvec=data_tvec, assumption_overrides=assumption_overrides, print_conditions=print_conditions, cell_formats=cell_formats)

        row_id += 2

    ws_transval.set_column(0, 0, ws_transval_width)
    ws_transval.set_column(2, 2, ws_transval_width)
    ws_transval.set_column(3, 4, assumption_width)

    # %% Program definitions sheet.
    if include_progs:
        ws_progmat.write(0, 0, 'Name')
        ws_progmat.write(0, 1, 'Abbreviation')

        # While writing default program names and labels, they are stored for future reference as well.
        prog_names_default = []
        prog_labels_default = []
        for prid in xrange(num_progs):
            prog_name = 'Program ' + str(prid + 1)
            prog_label = prog_name[0:4] + str(prid + 1)
            prog_names_default.append(prog_name)
            prog_labels_default.append(prog_label)
            ws_progmat.write(prid + 1, 0, prog_name)
            ws_progmat.write(prid + 1, 1, '=LEFT(%s,4)&"%i"' % (rc(prid + 1, 0), prid + 1), None, prog_label)
        ws_progmat.set_column(0, 1, ws_progmat_width)

        # Excel formulae strings that point to program names and labels are likewise stored.
        prog_names_formula = []
        prog_labels_formula = []
        for prid in xrange(num_progs):
            prog_names_formula.append("='%s'!%s" % (settings.databook['sheet_names']['progmat'], rc(prid + 1, 0, True, True)))
            prog_labels_formula.append("='%s'!%s" % (settings.databook['sheet_names']['progmat'], rc(prid + 1, 1, True, True)))

        makeTagMatrix(worksheet=ws_progmat, at_row=0, num_rows=num_progs, at_col=2, labels=pop_labels_default, formula_labels=pop_labels_formula)

    # %% Program details sheet.

    if include_progs:
        rows_imp = settings.databook['format']['programs']['max_lines_impact']
        num_progtypes = len(settings.progtype_name_labels.keys())

        row_id = 0
        for prid in xrange(num_progs):

            print_conditions = []
            prog_cov_prefix = ''
            prog_cov_header = '"Program Coverage"'
            prog_unit_header = 'CONCATENATE("Unit Cost Estimate",IF(%s="Fraction"," (Per 1%%)",""))' % rc(row_id + 2, 2)
            for k in xrange(num_progtypes):
                progtype_name = settings.progtype_name_labels.keys()[k]
                progtype_label = settings.progtype_name_labels[k]
                if 'special' in settings.progtype_specs[progtype_label]:
                    progtype_special = settings.progtype_specs[progtype_label]['special']
                    if progtype_special == 'cost_only':
                        prog_cov_prefix = '='
                        prog_cov_header = 'IF(%s="%s","%s",%s)' % (rc(row_id, 1), progtype_name, '...', prog_cov_header)
                        prog_unit_header = 'IF(%s="%s","%s",%s)' % (rc(row_id, 1), progtype_name, '...', prog_unit_header)

            ws_progval.write(row_id, 0, prog_names_formula[prid], None, prog_names_default[prid])
            ws_progval.write(row_id, 1, settings.progtype_name_labels.keys()[0])
            ws_progval.data_validation('%s' % rc(row_id, 1), {'validate': 'list', 'source': settings.progtype_name_labels.keys(), 'ignore_blank': False})
            ws_progval.write(row_id + 1, 0, '...')
            ws_progval.write(row_id + 2, 0, 'Cost-Coverage Details')
            ws_progval.write(row_id + 3, 0, '...')
            ws_progval.write(row_id + 4, 0, '...')
            ws_progval.write(row_id + 5, 0, '...')
            ws_progval.write(row_id + 2, 1, prog_cov_prefix + prog_cov_header)
            ws_progval.write(row_id + 3, 1, 'Program Funding')
            ws_progval.write(row_id + 4, 1, '=' + prog_unit_header, None, 'Unit Cost Estimate')  # As long as this expression starts with 'Unit Cost Estimate', it doesn't matter what it's preloaded as.
            ws_progval.write(row_id + 5, 1, 'Program Saturation')
            makeValueEntryArrayBlock(worksheet=ws_progval, at_row=row_id + 1, at_col=2, num_arrays=1, tvec=data_tvec, data_formats=['Number', 'Fraction'], print_conditions=['%s<>"..."' % rc(row_id + 2, 1)], cell_formats=cell_formats)
            makeValueEntryArrayBlock(worksheet=ws_progval, at_row=row_id + 3, at_col=2, num_arrays=1, tvec=data_tvec, data_formats=['USD'], no_header=True, cell_formats=cell_formats)
            makeValueEntryArrayBlock(worksheet=ws_progval, at_row=row_id + 4, at_col=2, num_arrays=1, tvec=data_tvec, assumption='1.0E+300', data_formats=['USD'], print_conditions=['%s<>"..."' % rc(row_id + 4, 1)], no_header=True, only_assumption=True, cell_formats=cell_formats)
            makeValueEntryArrayBlock(worksheet=ws_progval, at_row=row_id + 5, at_col=2, num_arrays=1, tvec=data_tvec, assumption='0.95', data_formats=['Fraction'], print_conditions=['%s<>"..."' % rc(row_id + 5, 1)], no_header=True, only_assumption=True, cell_formats=cell_formats)

            print_conditions = []
            ws_progval.write(row_id + 6, 0, 'Impact Attributes')
            for row_imp in xrange(rows_imp):
                if row_imp == 0: ws_progval.write(row_id + 6 + row_imp, 0, 'Impact Attributes')
                else: ws_progval.write(row_id + 6 + row_imp, 0, '...')
                print_conditions.append('%s<>"..."' % rc(row_id + 6 + row_imp, 1))
                super_string = '"..."'
                assumption_attname = super_string
                for k in xrange(num_progtypes):
                    progtype_name = settings.progtype_name_labels.keys()[k]
                    progtype_label = settings.progtype_name_labels[k]
                    try: attrib_name = settings.progtype_specs[progtype_label]['attribute_name_labels'].keys()[row_imp]
                    except: attrib_name = '...'
                    if k == 0:
                        assumption_attname = attrib_name
                    super_string = 'IF(%s="%s","%s",%s)' % (rc(row_id, 1), progtype_name, attrib_name, super_string)
                ws_progval.write(row_id + 6 + row_imp, 1, '=' + super_string, None, assumption_attname)
            makeValueEntryArrayBlock(worksheet=ws_progval, at_row=row_id + 6, at_col=2, num_arrays=rows_imp, tvec=data_tvec, data_formats=['Unique'], print_conditions=print_conditions, print_conditions_blank_preloading=False, no_header=True, cell_formats=cell_formats)


            row_id += 6 + rows_imp

        ws_progval.set_column(0, 1, ws_progval_width)
        ws_progval.set_column(2, 3, assumption_width)

    # %% Combine characteristics and parameters into one dictionary, then print out all elements to appropriate datasheets.
    all_specs = dcp(settings.charac_specs)
    all_specs.update(settings.linkpar_specs)
    specs_ordered = sorted(dcp(all_specs.keys()), key=lambda x: all_specs[x]['databook_order'])
    for def_label in specs_ordered:
        opt_suffix = ''  # A tag at the end of characteristic/parameter label to indicate something about it in the databook.
        if not all_specs[def_label]['databook_order'] < 0:  # Do not print out characteristic/parameter if databook order is negative.
            ws = ws_params[all_specs[def_label]['sheet_label']]['ws']
            row_id = ws_params[all_specs[def_label]['sheet_label']]['row_id']

            def_name = all_specs[def_label]['name']
            default_val = 0.0
            if 'default' in all_specs[def_label]:
                default_val = all_specs[def_label]['default']

            # Make the data-entry blocks.
            # First handle 'count' and 'percentage' characteristics, as well as transitions for junctions.
            if def_label in settings.charac_specs.keys():
                if 'entry_point' in settings.charac_specs[def_label]:
                    opt_suffix = settings.databook['suffix']['seed']
                else:
                    opt_suffix = settings.databook['suffix']['output']
                if 'plot_percentage' in all_specs[def_label]:
                    data_formats = ['Fraction']
                else:
                    data_formats = ['Number']
            elif def_label in settings.linkpar_specs.keys():
                opt_suffix = settings.databook['suffix']['par']
                data_formats = None
                if 'tag' in settings.linkpar_specs[def_label]:
                    for pair in settings.links[settings.linkpar_specs[def_label]['tag']]:
                        if 'junction' in settings.node_specs[pair[0]].keys():
                            data_formats = ['Proportion']
            makeValueEntryArrayBlock(worksheet=ws, at_row=row_id, at_col=1, num_arrays=num_pops, tvec=data_tvec, assumption=default_val, data_formats=data_formats, cell_formats=cell_formats)

            # Make the population references.
            ws.write(row_id, 0, def_name + opt_suffix)
            for pid in xrange(num_pops):
                row_id += 1
                ws.write(row_id, 0, pop_names_formula[pid], None, pop_names_default[pid])

            row_id += 2
            ws_params[all_specs[def_label]['sheet_label']]['row_id'] = row_id

    # Adjust widths for all custom sheets.
    for ws_label in ws_params:
        ws = ws_params[ws_label]['ws']
        ws.set_column(0, 0, name_width)
        ws.set_column(1, 2, assumption_width)


    workbook.close()


# %% Function for loading project databook

# NOTE: This needs so much quality-assurance testing. Need to ensure that input data sheet aligns with cascade settings.
def loadSpreadsheetFunc(settings, databook_path):
    ''' Load data spreadsheet into Project data dictionary. '''
    import os

    databook_path = os.path.abspath(databook_path)
    try: workbook = xlrd.open_workbook(databook_path)
    except: raise OptimaException('ERROR: Project data workbook was unable to be loaded from... %s' % databook_path)
    ws_pops = workbook.sheet_by_name(settings.databook['sheet_names']['pops'])
    ws_transmat = workbook.sheet_by_name(settings.databook['sheet_names']['transmat'])
    ws_transval = workbook.sheet_by_name(settings.databook['sheet_names']['transval'])

    # Contact sheet can be optional.
    ws_contact_exists = True
    try: ws_contact = workbook.sheet_by_name(settings.databook['sheet_names']['contact'])
    except:
        ws_contact_exists = False
        logging.warning('There is no "%s" sheet in project data workbook.' % settings.databook['sheet_names']['contact'])

    # Program sheets can be optional. If either does not exist, the other is ignored.
    ws_prog_exists = True
    try:
        ws_progmat = workbook.sheet_by_name(settings.databook['sheet_names']['progmat'])
        ws_progval = workbook.sheet_by_name(settings.databook['sheet_names']['progval'])
    except:
        ws_prog_exists = False
        logging.warning('One or two program-based sheets, "%s" and "%s", are excluded in the project data workbook.' % (settings.databook['sheet_names']['progmat'], settings.databook['sheet_names']['progval']))

    # Extra validation. Even if there are program sheets in the databook, they cannot be used if no program type framework has been loaded from the cascade workbook.
    if len(settings.progtype_specs.keys()) == 0:
        ws_prog_exists = False
        logging.warning('Any program-based sheets in the databook will be ignored, due to no program type framework being loaded in from the cascade workbook.')

    # Regarding cascade parameters and characteristics, store sheets and corresponding row ids for further writing.
    ws_params = odict()
    for custom_label in settings.databook['custom_sheet_names'].keys():
        ws_params[custom_label] = workbook.sheet_by_name(settings.databook['custom_sheet_names'][custom_label])

    # Only produce standard sheets if there are any parameters left over that are not allocated to custom sheets.
    if settings.make_sheet_characs:
        ws_params['charac'] = workbook.sheet_by_name(settings.databook['sheet_names']['charac'])
    if settings.make_sheet_linkpars:
        ws_params['linkpars'] = workbook.sheet_by_name(settings.databook['sheet_names']['linkpars'])

    # %% Population names sheet.
    data = odict()
    data['meta'] = dict()  # NOTE: Consider making this much more important in the future. Maybe store pop definitions here.
    data['pops'] = odict()
    data['pops']['name_labels'] = odict()
    data['pops']['label_names'] = odict()  # A reverse of the previous dictionary.
    data['pops']['ages'] = odict()
    data['pops']['dw'] = odict()
    data['pops']['life_exp'] = odict()
    for row_id in xrange(ws_pops.nrows):
        if row_id > 0:
            if ws_pops.cell_value(row_id, 0) not in ['']:
                if ws_pops.cell_value(row_id, 1) not in ['']:
                    pop_label = str(ws_pops.cell_value(row_id, 1))
                else:
                    pop_label = 'pop' + str(row_id)
                data['pops']['name_labels'][str(ws_pops.cell_value(row_id, 0))] = pop_label
                data['pops']['label_names'][pop_label] = str(ws_pops.cell_value(row_id, 0))
                age_min = ws_pops.cell_value(row_id, 2)
                age_max = ws_pops.cell_value(row_id, 3)
                if isinstance(age_min, Number) and isinstance(age_max, Number):
                    data['pops']['ages'][pop_label] = {'min':float(age_min), 'max':float(age_max), 'range':1 + float(age_max) - float(age_min)}

                if 5 < ws_pops.ncols:
                    if ws_pops.cell_value(row_id, 4) != '':
                        data['pops']['dw'][pop_label] = float(ws_pops.cell_value(row_id, 4))
                    if ws_pops.cell_value(row_id, 5) != '':
                        data['pops']['life_exp'][pop_label] = float(ws_pops.cell_value(row_id, 5))


    # %% Population contacts sheet.
    data['contacts'] = dict()
    data['contacts']['into'] = dict()
    data['contacts']['from'] = dict()
    for pop in data['pops']['label_names'].keys():
        data['contacts']['into'][pop] = dict()
        data['contacts']['from'][pop] = dict()
    if ws_contact_exists:
        for row_id in xrange(ws_contact.nrows):
            for col_id in xrange(ws_contact.ncols):
                if row_id > 0 and col_id > 0:
                    source = str(ws_contact.cell_value(row_id, 0))
                    target = str(ws_contact.cell_value(0, col_id))
                    val = ws_contact.cell_value(row_id, col_id)
                    if val != '' and float(val) != 0:
                        data['contacts']['into'][target][source] = float(val)
                        data['contacts']['from'][source][target] = float(val)
    else:
        # Self-connections are the default if there is no contact worksheet. These can be turned off in an actual contacts sheet.
        for pop in data['pops']['label_names'].keys():
            data['contacts']['into'][pop][pop] = 1.0
            data['contacts']['from'][pop][pop] = 1.0
        logging.warning('No "%s" sheet means population groups only interact with themselves by default.' % settings.databook['sheet_names']['contact'])



    # %% Inter-population transfer definitions sheet.
    # Migration matrices must be divided from each other by an empty row.
    data['transfers'] = odict()
    mig_specified = False
    mig_type = None
    for row_id in xrange(ws_transmat.nrows):
        zero_col = ws_transmat.cell_value(row_id, 0)

        # If parser finds an empty row (technically empty first cell) migration-parsing resets, ready for the next custom type.
        if mig_specified and zero_col in ['']:
            mig_specified = False

        # Parse through migration matrix.
        if mig_specified:
            pop_source = str(zero_col)
            for col_id in xrange(ws_transmat.ncols):
                if col_id > 0:
                    val = ws_transmat.cell_value(row_id, col_id)
                    if val in ['y']:
                        pop_target = str(ws_transmat.cell_value(0, col_id))
                        if pop_source not in data['transfers'][mig_type].keys():
                            data['transfers'][mig_type][pop_source] = odict()
                        data['transfers'][mig_type][pop_source][pop_target] = odict()
                        if mig_type == 'aging':
                            if len(data['transfers'][mig_type][pop_source].keys()) > 1:
                                raise OptimaException('ERROR: There are too many outgoing "%s" transitions listed for population "%s".' % (mig_type, pop_source))

        # First row after a blank one must contain the new migration type as its first element. Parser re-activated.
        if not mig_specified and zero_col not in ['']:
            mig_specified = True
            mig_type = str(zero_col).lower().replace(' ', '_')
            data['transfers'][mig_type] = odict()

    # %% Inter-population transfer details sheet.
    mig_specified = False
    mig_type = None
    array_id = 0
    for row_id in xrange(ws_transval.nrows):
        zero_col = ws_transval.cell_value(row_id, 0)

        # If parser finds an empty row (technically empty first cell) migration-parsing resets, ready for the next custom type.
        if mig_specified and zero_col in ['']:
            mig_specified = False

        # Parse through migration-specific data rows.
        if mig_specified:
            pop_source = str(zero_col)
            if pop_source not in ['', '...']:
                pop_source_label = data['pops']['name_labels'][pop_source]
                pop_target = ws_transval.cell_value(row_id, 2)
                pop_target_label = data['pops']['name_labels'][pop_target]
                list_t = []
                list_y = []
                for col_id in xrange(ws_transval.ncols):
                    if col_id == 3:
                        data['transfers'][mig_type][pop_source_label][pop_target_label]['y_format'] = str(ws_transval.cell_value(row_id, col_id)).lower()
                    if col_id > 3 and isinstance(ws_transval.cell_value(row_id, col_id), Number):
                        val = ws_transval.cell_value(row_id, col_id)

                        list_y.append(float(val))
                        if not isinstance(ws_transval.cell_value(row_id - 1 - array_id, col_id), Number):
                            list_t.append(float(ws_transval.cell_value(row_id - 1 - array_id, col_id + 2)))
                            break
                        else:
                            list_t.append(float(ws_transval.cell_value(row_id - 1 - array_id, col_id)))

                data['transfers'][mig_type][pop_source_label][pop_target_label]['t'] = np.array(list_t)
                data['transfers'][mig_type][pop_source_label][pop_target_label]['y'] = np.array(list_y)
                data['transfers'][mig_type][pop_source_label][pop_target_label]['y_factor'] = project_settings.DEFAULT_YFACTOR  # NOTE: Quick hack to make migrations work with calibration branch.
            array_id += 1



        # First row after a blank one must contain the new migration type as its first element. Parser re-activated.
        if not mig_specified and zero_col not in ['']:
            mig_specified = True
            mig_type = str(zero_col).lower().replace(' ', '_')
            array_id = 0

    # %% Program definitions sheet.
    data['meta']['progs'] = dict()
    data['meta']['progs']['name_labels'] = odict()
    data['meta']['progs']['label_names'] = odict()  # A reverse of the previous dictionary.
    data['progs'] = odict()
    if ws_prog_exists:
        for row_id in xrange(ws_progmat.nrows):
            if row_id > 0:
                if ws_progmat.cell_value(row_id, 0) not in ['']:
                    if ws_progmat.cell_value(row_id, 1) not in ['']:
                        prog_label = str(ws_progmat.cell_value(row_id, 1))
                    else:
                        prog_label = 'prog' + str(row_id)
                    data['meta']['progs']['name_labels'][str(ws_progmat.cell_value(row_id, 0))] = prog_label
                    data['meta']['progs']['label_names'][prog_label] = str(ws_progmat.cell_value(row_id, 0))
                    data['progs'][prog_label] = dict()
                    data['progs'][prog_label]['name'] = str(ws_progmat.cell_value(row_id, 0))  # Label name linkage is also in metadata, but redundancy here is ok for now.
                    data['progs'][prog_label]['target_pops'] = []
                    for col_id in xrange(ws_progmat.ncols):
                        if col_id > 1:
                            if str(ws_progmat.cell_value(row_id, col_id)) == 'y':
                                data['progs'][prog_label]['target_pops'].append(str(ws_progmat.cell_value(0, col_id)))

    # %% Program details sheet.
    if ws_prog_exists:
        prog_specified = False
        prog_name = None
        prog_label = None
        progtype_name = None
        progtype_label = None
        prog_row_id = 0
        temp = odict()
        for row_id in xrange(ws_progval.nrows):
            zero_col = ws_progval.cell_value(row_id, 0)

            # If parser finds an empty row (technically empty first cell) program-parsing resets, ready for the next custom program.
            if prog_specified and zero_col in ['']:
                prog_specified = False

            # Parse through program-specific data rows.
            if prog_specified:
                get_data = False
                important_col = str(ws_progval.cell_value(row_id, 1))
                if important_col == 'Program Coverage':
                    tag = 'cov'
                    get_data = True
                if important_col == 'Program Funding':
                    tag = 'cost'
                    get_data = True
                if important_col == 'Program Saturation':
                    tag = 'sat'
                    get_data = True
                if important_col.startswith('Unit Cost Estimate'):
                    estim = str(ws_progval.cell_value(row_id, 3))
                    if not estim in ['']:
                        data['progs'][prog_label]['unit_cost'] = float(estim)
                if important_col in settings.progtype_specs[progtype_label]['attribute_name_labels'].keys():
                    tag = settings.progtype_specs[progtype_label]['attribute_name_labels'][important_col]
                    get_data = True
                if get_data:
                    for col_id in xrange(2, ws_progval.ncols):
#                        print tag
#                        print ws_progval.cell_value(row_id, col_id)
                        if col_id == 2:
                            if tag in ['cov', 'cost']:
                                data['progs'][prog_label]['%s_format' % tag] = str(ws_progval.cell_value(row_id, col_id)).lower()
                            else:  # Assume attributes have unique formats. No need to store at this time.
                                pass
#                                data['progs'][prog_label]['attributes']['%s_format' % tag] = str(ws_progval.cell_value(row_id, col_id)).lower()
                        elif col_id > 2 and isinstance(ws_progval.cell_value(row_id, col_id), Number):
                            val = ws_progval.cell_value(row_id, col_id)

                            if not isinstance(ws_progval.cell_value(prog_row_id + 1, col_id), Number):
                                tval = str(ws_progval.cell_value(prog_row_id + 1, col_id + 2))
                                temp[prog_label]['t_assumption'] = tval
                                if ('%s_assumption' % tag) not in temp[prog_label].keys():
                                    temp[prog_label]['%s_assumption' % tag] = dict()
                                temp[prog_label]['%s_assumption' % tag] = val
                                break
                            else:
                                tval = str(ws_progval.cell_value(prog_row_id + 1, col_id))
                                if tval not in temp[prog_label].keys(): temp[prog_label][tval] = dict()
                                temp[prog_label][tval][tag] = val

                        # special rule: assumption may also specify string referring to another program
                        # IMPORTANT NOTE: referenced programs must have been defined before!
                        elif col_id == 3 and ws_progval.cell_value(row_id, col_id) in temp.keys():
                            val = ws_progval.cell_value(row_id, col_id)
                            tval = str(ws_progval.cell_value(prog_row_id + 1, col_id + 2))
                            temp[prog_label]['t_assumption'] = tval
                            if ('%s_assumption' % tag) not in temp[prog_label].keys():
                                temp[prog_label]['%s_assumption' % tag] = dict()
                            temp[prog_label]['%s_assumption' % tag] = val
                            break
                        elif col_id == 3:
                            attr = ws_progval.cell_value(row_id, col_id - 2)
                            if attr not in settings.progtype_specs[progtype_label]['attribute_name_labels']:
                                continue
                            lab = settings.progtype_specs[progtype_label]['attribute_name_labels'][attr]
                            if lab.startswith('$ref_'):
                                logger.warn('If "{}" is to refer to the program '.format(prog_label) +
                                            '"{}", it must be defined '.format(ws_progval.cell_value(row_id, col_id)) +
                                            'in the databook before {}.'.format(prog_label) +
                                            ' Ignoring entry; this may result in unpredictable behaviour.')



            # First row after a blank one must contain the new program name as its first element. Parser re-activated.
            if not prog_specified and zero_col not in ['']:
                prog_specified = True
                prog_name = str(zero_col)
                prog_label = data['meta']['progs']['name_labels'][prog_name]
                prog_row_id = row_id
                progtype_name = str(ws_progval.cell_value(row_id, 1))
                progtype_label = settings.progtype_name_labels[progtype_name]
                data['progs'][prog_label]['prog_type'] = progtype_label
                data['progs'][prog_label]['attributes'] = odict()
                temp[prog_label] = odict()

        # Cost coverage data must be converted into time, cost and coverage lists.
        for prog_label in data['progs'].keys():
            list_t = []
            list_cost = []
            list_cov = []
            list_sat = []
            progtype_label = data['progs'][prog_label]['prog_type']
            special_tag = None
            if 'special' in settings.progtype_specs[progtype_label]:
                special_tag = settings.progtype_specs[progtype_label]['special']
                if special_tag == 'cost_only':
                    temp[prog_label]['cov_assumption'] = np.nan
                    data['progs'][prog_label]['cov_format'] = None
            num_attribs = len(settings.progtype_specs[progtype_label]['attribute_name_labels'])
            list_attribs = [[] for x in xrange(num_attribs)]

            # Run through all non-assumption temp-dict keys, i.e. years.
            for tval in sorted(temp[prog_label].keys()):
                try: t = float(tval)
                except: continue
                list_t.append(t)
                cost = np.nan
                cov = np.nan

                # If cost or coverage is missing for a year for which the other exists, make use of an assumption or, failing that, np.nan.
                if 'cost' in temp[prog_label][tval]:
                    cost = float(temp[prog_label][tval]['cost'])
                elif 'cost_assumption' in temp[prog_label]:
                    cost = float(temp[prog_label]['cost_assumption'])
                if 'cov' in temp[prog_label][tval]:
                    cov = float(temp[prog_label][tval]['cov'])
                elif 'cov_assumption' in temp[prog_label]:
                    cov = float(temp[prog_label]['cov_assumption'])
                if 'sat' in temp[prog_label][tval]:
                    sat = float(temp[prog_label][tval]['sat'])
                elif 'sat_assumption' in temp[prog_label]:
                    sat = float(temp[prog_label]['sat_assumption'])
                else:
                    sat = None

                list_cost.append(cost)
                list_cov.append(cov)
                list_sat.append(sat)

                for aid in xrange(num_attribs):
                    attrib_label = settings.progtype_specs[progtype_label]['attribute_name_labels'][aid]
                    attrib = np.nan
                    if attrib_label in temp[prog_label][tval]:
                        attrib = float(temp[prog_label][tval][attrib_label])
                    elif attrib_label + '_assumption' in temp[prog_label]:
                        if attrib_label.startswith('$ref_'): # special case of programs cross-referencing other programs
                            attrib = temp[prog_label][attrib_label + '_assumption']
                        else:
                            attrib = float(temp[prog_label][attrib_label + '_assumption'])
                    list_attribs[aid].append(attrib)


            if len(list_t) == 0:  # In the case that only assumptions are given...
                try:
                    list_t.append(float(temp[prog_label]['t_assumption']))
                    list_cost.append(float(temp[prog_label]['cost_assumption']))
                    list_cov.append(float(temp[prog_label]['cov_assumption']))
                    list_sat.append(float(temp[prog_label]['sat_assumption']))
#                    print num_attribs
#                    print list_attribs
#                    print temp[prog_label]
                    for aid in xrange(num_attribs):
                        attrib_label = settings.progtype_specs[progtype_label]['attribute_name_labels'][aid]
                        list_attribs[aid].append(float(temp[prog_label][attrib_label + '_assumption']))

                except:
                    raise OptimaException('ERROR: There is incomplete cost-coverage data provided in the databook for program "%s".' % prog_label)
            data['progs'][prog_label]['t'] = np.array(list_t)
            data['progs'][prog_label]['cost'] = np.array(list_cost)
            data['progs'][prog_label]['cov'] = np.array(list_cov)
            data['progs'][prog_label]['sat'] = np.array(list_sat)
            for aid in xrange(num_attribs):
                attrib_label = settings.progtype_specs[progtype_label]['attribute_name_labels'][aid]
                data['progs'][prog_label]['attributes'][attrib_label] = np.array(list_attribs[aid])


    # %% Gather data value inputs for epidemic characteristics and cascade parameters sheets, be they custom or standard.

    name_to_label = dcp(settings.charac_name_labels)
    name_to_label.update(settings.linkpar_name_labels)
    data['characs'] = odict()
    data['linkpars'] = odict()

    for ws_label in ws_params.keys():
        ws = ws_params[ws_label]
        current_def_name = None
        current_def_label = None
        data_label = None
        pop_id = 0
        for row_id in xrange(ws.nrows):
            val = str(ws.cell_value(row_id, 0))
            if val in ['']:
                current_def_name = None
                pop_id = 0
            elif current_def_name is None:
                for suffix in settings.databook['suffix'].values():  # Remove any optional suffixes. NOTE: Might need validation in cascade sheet in case somebody uses those suffixes.
                    if val.endswith(suffix):
                        val = val[:-len(suffix)]
                current_def_name = val
                current_def_label = name_to_label[val]
                # TODO: remove dependency in this clause: the following line would fail if current_def_label is in both charac_specs and linkpar_specs
                if current_def_label in settings.charac_specs.keys():
                    data_label = 'characs'
                elif current_def_label in settings.linkpar_specs.keys():
                    data_label = 'linkpars'
                else:
                    raise OptimaException('ERROR: Data entry in the "%s" sheet includes "%s", which has no defined characteristic/parameter specifications.' % (ws_label, current_def_name))
                data[data_label][current_def_label] = odict()
            else:
                current_pop_label = data['pops']['name_labels'][val]
                if current_pop_label != data['pops']['label_names'].keys()[pop_id]:
                    raise OptimaException('ERROR: Somewhere in the "%s" parameters sheet, populations are not ordered as in the population definitions sheet.' % data_label)
                data[data_label][current_def_label][current_pop_label] = odict()

                # Run through the rows beneath the year range, but only if there is not a number in the cell corresponding to assumption.
                # NOTE: Somewhat hard-coded. Improve.
                list_t = []
                list_y = []
#                list_comment = []

                for col_id in xrange(ws.ncols):
                    if col_id == 1:  # TODO: DK to remove hard-coded variables in code ...
                        data[data_label][current_def_label][current_pop_label]['y_format'] = str(ws.cell_value(row_id, col_id)).lower()
                    if col_id > 1 and isinstance(ws.cell_value(row_id, col_id), Number):
                        list_y.append(float(ws.cell_value(row_id, col_id)))
#                        # Store comments attached to databook parameter value cells within project data.
#                        # Requires an .xls file, which does not support the level of generic equation nesting the current design requires.
#                        # TODO: Consider sparse storage rather than bogging data down with a lot of None.
##                        if len(ws.cell_note_map) > 0: print ws.cell_note_map
#                        if (row_id,col_id) in ws.cell_note_map:
#                            list_comment.append(str(ws.cell_note_map[(row_id,col_id)].text))
#                        else:
#                            list_comment.append(None)
                        if not isinstance(ws.cell_value(row_id - 1 - pop_id, col_id), Number):
                            list_t.append(float(ws.cell_value(row_id - 1 - pop_id, col_id + 2)))
                            break
                        else:
                            list_t.append(float(ws.cell_value(row_id - 1 - pop_id, col_id)))
                data[data_label][current_def_label][current_pop_label]['t'] = np.array(list_t)
                data[data_label][current_def_label][current_pop_label]['y'] = np.array(list_y)
#                data[data_label][current_def_label][current_pop_label]['comment'] = np.array(list_comment)
                # if there is a corresponding data value for y_factor already in cascacade structure in settings, use that; else default to settings value
                try:
                    if data_label == 'linkpars':
                        data[data_label][current_def_label][current_pop_label]['y_factor'] = settings.linkpar_specs[current_def_label]['y_factor']
                    elif data_label == 'characs':
                        data[data_label][current_def_label][current_pop_label]['y_factor'] = settings.charac_specs[current_def_label]['y_factor']
                except:
                    logger.info("Couldn't read y_factor for parameter '%s'" % current_def_label)
                    data[data_label][current_def_label][current_pop_label]['y_factor'] = project_settings.DEFAULT_YFACTOR

                pop_id += 1

    # All parameters must be defined whether they are in the project databook or not.
    for label in settings.linkpar_specs.keys():
        if label not in data['linkpars'].keys():
            def_format = 'Fraction'.lower()  # NOTE: Hard-coded format assumption. Improve at some stage when allowing for default formats.
            def_val = np.nan
            if 'default' in settings.linkpar_specs[label]:
                def_val = settings.linkpar_specs[label]['default']
            if 'f_stack' in settings.linkpar_specs[label]:
                def_val = np.nan
            else:
                print('WARNING: Project data sheet does not contain required cascade parameter "%s".\n         Using default format "%s" and default value %f.' % (label, def_format, def_val))
            data['linkpars'][label] = odict()
            for pop in data['pops']['label_names'].keys():
                data['linkpars'][label][pop] = odict()
                data['linkpars'][label][pop]['y_format'] = def_format
                data['linkpars'][label][pop]['t'] = np.array([settings.tvec_start])
                data['linkpars'][label][pop]['y'] = np.array([float(def_val)])
                data['linkpars'][label][pop]['y_factor'] = settings.linkpar_specs[label]['y_factor']  # This shouldn't be overwritten as the DEFAULT_YFACTOR if it's already defined

    validation_level = settings.validation['databook_validation']
    validation = databookValidation(data=data)
    if validation:
        pass  # no inconsistencies detected
    elif validation_level == project_settings.VALIDATION_ERROR:
        raise OptimaException('ERROR: Databook entries incomplete or mismatched, please look at log for details')
    elif validation_level == project_settings.VALIDATION_WARN or validation_level == project_settings.VALIDATION_AVERT:
        logger.warn("Validating databook: possible inconsistencies observed (see log for details)")
    else:  # we ignore
        pass
    return data

def createFormats(workbook):
    formats = odict()
    for fkey, fmt in WB_FORMATS.iteritems():
        formats[fkey] = workbook.add_format(fmt)

    return workbook, formats

def __addCharacteristicData(data, charac_label, pop_label, ts, ys, y_format, y_factor=1.):
    """
    
    
    """
    if charac_label not in data['characs'].keys():
        data['characs'][charac_label] = odict()

    if pop_label not in data['characs'][charac_label].keys():
        data['characs'][charac_label][pop_label] = odict()


    data['characs'][charac_label][pop_label]['t'] = ts
    data['characs'][charac_label][pop_label]['y'] = ys
    data['characs'][charac_label][pop_label]['y_format'] = y_format
    data['characs'][charac_label][pop_label]['y_factor'] = y_factor

    return data

def getEmptyData():
    """
    Create empty data structure, following structure returned in loadSpreadsheetFunc()
    
    This skeleton structure is used elsewhere within the framework, as a placeholder.
    
    Clarification:
        pops        population definitions, with name, label and ages
        link_pars   inter-compartment (intra-population) transfer parameters 
        characs     characteristics
        transfers   inter-populations transfer parameters
    """
    data = odict()
    data['pops'] = odict()
    data['pops']['name_labels'] = odict()
    data['pops']['label_names'] = odict()

    data['contacts'] = dict()
    data['characs'] = odict()
    data['transfers'] = odict()
    data['linkpars'] = odict()
    return data

def databookValidation(data=None):
    '''
    Validate data fields within the databook:

    Current Operation:
        1. Checks to ensure data entered is in line:
             a) Format type: Fraction - Ensure data entered for parameters is not negative or greater than one
             b) Format type: Number - Ensure data entered for parameters is not negative or in range (0-1)
        2. Checks to ensure that age in in the correct range:
             a) Ensures Minimum age is lower than Maximum age
             b) Ensure that age is a non-negative number
    
    Output:
        Output is simply a complete log of errors found in relavant databook
        
    '''
    validation = True
    label = 'y'
    for key in data:
        for attribute in data[key]:
            if attribute in ['name_labels', 'label_names', 'life_exp']: pass
            else:
                for pop in data[key][attribute]:
                    if key == 'transfers':
                        for subpop in data[key][attribute][pop]:
                            for loop in range (len(data[key][attribute][pop][subpop][label])):
                                validation = validateFormatType(data[key][attribute][pop][subpop], label, loop, key, attribute, pop, validation)
                    elif key == 'pops':
                        if attribute == 'ages':
                            if data[key][attribute][pop]['max'] <= data[key][attribute][pop]['min'] or data[key][attribute][pop]['max'] <= 0 or data[key][attribute][pop]['min'] < 0:
                                logging.warning('Minimum and maximum age is defined incorrectly for Population: %s' % (pop))
                                validation = False
                        elif attribute == 'dw':
                            if not 0. <= data[key][attribute][pop] <= 1.:
                                logging.warning('Disability weight for Population %s must be in [0,1]!' % (pop))
                                validation = False
                        else:
                            logging.warning('Invalid attribute "%s" in population %s encountered' % (attribute, pop))
                            validation = False
                    elif key in ['meta', 'contacts', 'progs']:
                        # NOTE: Similar validation to be filled in for these keys at some point if required.
                        pass
                    else:
                        for loop in range (len(data[key][attribute][pop][label])):
                            validation = validateFormatType(data[key][attribute][pop], label, loop, key, attribute, pop, validation)
    return validation

def validateFormatType(data_to_validate, label, loop, key, attribute, pop, validation):
    '''
    Helper function which is called from databookValidation function. 
    It loops through the databook entries and makes sure they conform to the 
    format_type specified (fraction or number)
    '''
    if key == 'transfers': key = 'Transfer Details: '
    elif key == 'characs': key = 'Characteristic: '
    elif key == 'linkpars': key = 'Parameter: '

    if data_to_validate['y_format'] == 'fraction':
      if data_to_validate[label][loop] > 1. or data_to_validate[label][loop] < 0.:
          logging.warning('Please verify databook under %s%s and population %s as a number greater than 1 or negative number was entered for definition type "fraction" for Year: %i, value entered: %0.1f' % (key, attribute, pop, data_to_validate['t'][loop], data_to_validate['y'][loop]))
          validation = False
    elif data_to_validate['y_format'] == 'number':
      if data_to_validate[label][loop] < 0.:
          logging.warning('Please verify databook under %s%s and population %s as a negative number was entered for definition type "number" for Year: %i, value entered: %0.1f' % (key, attribute, pop, data_to_validate['t'][loop], data_to_validate['y'][loop]))
          validation = False
#       elif data_to_validate[label][loop] > 0. and data_to_validate[label][loop] < 1.:
#           logging.warning('Please verify databook under %s%s and population %s as a fraction or a negative number was entered for definition type "number" for Year: %i, value entered: %0.1f' % (key, attribute, pop, data_to_validate['t'][loop], data_to_validate['y'][loop]))
#           validation = False
    return validation



