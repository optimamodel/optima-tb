
import sys
sys.path.append('..')

from optima_tb.utils import odict

def new_countrybook():
    """
    Return a default Countrybook dictionary

    Typical workflow would be to modify these fields further, and then save 
    via write_countrybook()
    """

    countrybook = odict()
    countrybook['sheet_names'] = odict()
    countrybook['sheet_names']['populations'] = 'Populations'
    countrybook['sheet_names']['demographics'] = 'Demographics'
    countrybook['sheet_names']['prevalence'] = 'TB incidence & prevalence'
    countrybook['sheet_names']['notifications'] = 'Notifications'
    countrybook['sheet_names']['smear'] = 'Smear status'
    countrybook['sheet_names']['comorbidity'] = 'Comorbidities'
    countrybook['sheet_names']['morbidity'] = 'TB morbidity'
    countrybook['sheet_names']['testing_treatment_latent'] = 'Outcomes - Latent TB'
    countrybook['sheet_names']['testing_treatment_active'] = 'Outcomes - Active TB'
    # countrybook['sheet_names']['programs'] = 'Programs'
    # countrybook['sheet_names']['cost_coverage'] = 'Cost and coverage'
    # countrybook['sheet_names']['unitcost'] = 'Unit costs'
    # countrybook['sheet_names']['poptransitions'] = 'Population transitions'

    # headers for special sheets (i.e. headers aren't years)
    countrybook['headers'] = odict()
    countrybook['headers']['populations'] = ['Name', 'Minimum Age', 'Maximum Age']
    # countrybook['headers']['programs'] = ['Name', 'Short name', 'Intervention class', 'Coverage indicator (annual)', 'Duration of treatment (days per person on average)', 'Frequency of intervention (in years)']

    # labels for each sheet
    countrybook['labels'] = {'populations': '',
                                  'demographics':'',
                                  'prevalence'     : 'Estimated incidence and prevalence of TB',
                                  'notifications'  : 'TB notifications: please enter the number of notifications of active TB per population and drug-resistant strain per year',
                                  'smear'          : '',
                                  'comorbidity':'Comorbidities',
                                  'morbidity' : 'Number of TB-related deaths per year',
                                  'testing_treatment_latent': 'Testing and treatment outcomes for latent TB',
                                  'testing_treatment_active'        : 'Testing and treatment outcomes for latent TB'
                                  }

    # info
    countrybook['disaggregations'] = odict()
    # other values
    countrybook['disaggregations']['strains'] = ['DS-TB', 'MDR-TB', 'XDR-TB']
    countrybook['disaggregations']['smears'] = ['Smear-', 'Smear+']  # also potentially 'Extrapulmonary'
    # disagg for non-total smear and strain
    countrybook['disaggregations']['nt_strains'] = ['DS-TB', 'MDR-TB', 'XDR-TB']
    countrybook['disaggregations']['nt_smears'] = ['Smear-', 'Smear+']  # also potentially 'Extrapulmonary'
    countrybook['disaggregations']['populations'] = []  # determined dynamically at runtime
    countrybook['disaggregations']['regimens'] = ['DS-TB regimen', 'MDR-TB regimen', 'XDR-TB regimen']
    countrybook['disaggregations']['programs'] = []  # determined dynamically at runtime
    countrybook['disaggregations']['total_pop'] = ['Total population']
    countrybook['disaggregations']['total_birth'] = ['Total number of births']

    # for univalue sheets, includes information on how data should be disaggregated
    countrybook['sheet_classes'] = odict()
    countrybook['sheet_classes']['univalue'] = odict()
    # countrybook['sheet_classes']['univalue']['population_sizes'] = ['populations']
    # countrybook['sheet_classes']['univalue']['total_cases'] = ['populations', 'smears', 'strains']
    countrybook['sheet_classes']['univalue']['notifications'] = ['populations', 'smears', 'strains']
    countrybook['sheet_classes']['univalue']['morbidity'] = ['populations', 'strains']

    # sheet specific values
    countrybook['sheet_values'] = odict()

    countrybook['sheet_values']['demographics'] = odict()
    countrybook['sheet_values']['demographics']['Population size per year'] = ['populations']
    countrybook['sheet_values']['demographics']['Number of births per year'] = ['total_birth']
    countrybook['sheet_values']['demographics']['Percentage of people vaccinated per year'] = ['populations']
    countrybook['sheet_values']['demographics']['Percentage of people who die from non-TB related causes per year'] = ['populations']

    countrybook['sheet_values']['prevalence'] = odict()
    countrybook['sheet_values']['prevalence']['Estimated TB incidence (per 100,000)'] = ['populations']
    countrybook['sheet_values']['prevalence']['Estimated active TB prevalence'] = ['populations']
    countrybook['sheet_values']['prevalence']['Estimated MDR TB prevalence'] = ['populations']
    countrybook['sheet_values']['prevalence']['Estimated XDR TB prevalence'] = ['populations']
    countrybook['sheet_values']['prevalence']['Estimated latent TB prevalence'] = ['populations']

    countrybook['sheet_values']['smear'] = odict()
    countrybook['sheet_values']['smear']['Smear status by drug-resistant strain'] = ['nt_smears', 'nt_strains']
    countrybook['sheet_values']['smear']['Smear status by population'] = ['nt_smears', 'populations']


    countrybook['sheet_values']['comorbidity'] = odict()
    # countrybook['sheet_values']['comorbidity']['HIV prevalence'] = ['populations', 'smears', 'strains']
    # countrybook['sheet_values']['comorbidity']['Diabetes prevalence'] = ['populations', 'smears', 'strains']

    countrybook['sheet_values']['testing_treatment_latent'] = odict()
    countrybook['sheet_values']['testing_treatment_latent']['Percentage of population tested for latent TB per year'] = ['populations']
    countrybook['sheet_values']['testing_treatment_latent']['Number of people initiating treatment for latent TB per year'] = ['populations']
    countrybook['sheet_values']['testing_treatment_latent']['Number of people lost to follow up for latent TB per year'] = ['populations']
    countrybook['sheet_values']['testing_treatment_latent']['Number of people who succesfully completed treatment for latent TB'] = ['populations']

    countrybook['sheet_values']['testing_treatment_active'] = odict()
    countrybook['sheet_values']['testing_treatment_active']['Percentage of population tested for active TB per year'] = ['populations']
    countrybook['sheet_values']['testing_treatment_active']['Number of people initiating treatment for active TB per year'] = ['regimens', 'populations']
    countrybook['sheet_values']['testing_treatment_active']['Number of people lost to follow up for active TB per year'] = ['regimens', 'populations']
    countrybook['sheet_values']['testing_treatment_active']['Number of people who failed treatment for active TB'] = ['regimens', 'populations']
    countrybook['sheet_values']['testing_treatment_active']['Number of people who successfully completed treatment for active TB'] = ['regimens', 'populations']


    countrybook['constants'] = {'spacing_interpopulation':2,
                                     'spacing_intrapopulation':1,
                                     'spacing_interproperty'  :4,
                                     'spacing_multivalue_label':2,
                                     'total_strains': 'Total',  # All strains',
                                     'total_smears' : 'Total',
                                     'num_default_programs':28,
                                     'row_index_start':2,  # for when there are no disaggregations, etc.
                                     'col_index_start':1}  #

    return countrybook


if __name__ == '__main__':

    from optima_tb.spreadsheet import write_countrybook

    # Example usage
    cb = new_countrybook()
    pop_names = ['0-4 years', '5-14 years', '15+ years']
    write_countrybook(cb,2000,2016,len(pop_names),pop_names,filename='example.xlsx')
