from IPython import get_ipython
import matplotlib

ipython = get_ipython()
if ipython is not None:
    ipython.magic('load_ext autoreload')
    ipython.magic('autoreload 2')
matplotlib.use('Qt5Agg')

from optima_tb.project import Project
from optima_tb.defaults import defaultOptimOptions
from optima_tb.reconciliation import reconcile
from optima_tb.utils import odict
import optima_tb.plotting as oplt
import matplotlib.pyplot as plt
import os

#Initialise conditions
parset_name='default'
progset_name = 'default'
reconcile_for_year=2015
unitcost_sigma=0.5
attribute_sigma=0.0
budget_sigma = 0.0
impact_pars=None

#Setup project
cascade = '../tb-ucl-analyses/belarus/Cascade Spreadsheets/cascade-belarus.xlsx'
databook = '../tb-ucl-analyses/belarus/Databook Spreadsheets/databook-belarus.xlsx'
proj= Project(name = 'Belarus-newmodelfit', cascade_path = cascade, validation_level = 'avert', plotting_level = 'dev')
proj.loadSpreadsheet(databook_path=databook)
parset = proj.makeParset(name=parset_name)
progset = proj.makeProgset(name=progset_name)
#
#
# updated_vals = {
# ('HF XDR-TB', 'unit_cost' ) : 8566.66882845,
# ('HF XDR-TB', 'budget') : 3722607.0,
# ('HF XDR-TB', 'adher' ) : 0.399,
# ('HF XDR-TB', 'effic' ) : 0.25,
# ('INVO XDR-TB', 'unit_cost' ) : 75220.19999999445,
# ('INVO XDR-TB', 'budget') : 2849250.0,
# ('INVO XDR-TB', 'adher' ) : 1.0,
# ('INVO XDR-TB', 'effic' ) : 0.25,
# ('BCG (0-4)', 'unit_cost' ) : 1.9799999999999982,
# ('BCG (0-4)', 'budget') : 528000.0,
# ('BCG (0-4)', 'effect') : 0.12447998046875,
# ('Symp-Diag', 'unit_cost' ) : 53.55725,
# ('Symp-Diag', 'budget') : 770489.14,
# ('Symp-Diag', 'sens') : 0.11520886551428454,
# ('INVO MDR-TB', 'unit_cost' ) : 25073.400000000067,
# ('INVO MDR-TB', 'budget') : 3595525.56,
# ('INVO MDR-TB', 'adher' ) : 1.0,
# ('INVO MDR-TB', 'effic' ) : 0.75,
# ('Mass-Screening' , 'unit_cost' ) : 5072.745234375,
# ('Mass-Screening' , 'budget') : 9758596.0,
# ('Mass-Screening' , 'sens') : 0.00024070051365995783,
# ('IPT-PLHIV', 'unit_cost' ) : 19.004999999999374,
# ('IPT-PLHIV', 'budget') : 2477.0,
# ('IPT-PLHIV', 'adher' ) : 0.82,
# ('IPT-PLHIV', 'effic' ) : 0.33,
# ('IPT Gen Pop', 'unit_cost' ) : 18.9555078125,
# ('IPT Gen Pop', 'budget') : 10575.0,
# ('IPT Gen Pop', 'adher' ) : 0.5,
# ('IPT Gen Pop', 'effic' ) : 0.39,
# ('Pall-Care', 'unit_cost' ) : 2643.0600000000004,
# ('Pall-Care', 'budget') : 2894780.0,
# ('Pall-Care', 'adher' ) : 0.5,
# ('Pall-Care', 'effic' ) : 0.0,
# ('HF MDR-TB', 'unit_cost' ) : 9474.102,
# ('HF MDR-TB', 'budget') : 4251974.76,
# ('HF MDR-TB', 'adher' ) : 0.5760000000000001,
# ('HF MDR-TB', 'effic' ) : 0.75,
# ('HF DS-TB' , 'unit_cost' ) : 4306.0380000000005,
# ('HF DS-TB' , 'budget') : 5098606.1,
# ('HF DS-TB' , 'adher' ) : 0.897867187,
# ('HF DS-TB' , 'effic' ) : 0.85,
# ('C-Tracing', 'unit_cost' ) : 812.7678906250001,
# ('C-Tracing', 'budget') : 7046.0,
# ('C-Tracing', 'sens') : 0.003437188474607226}
#
# updated_vals = {
# ('HF XDR-TB', 'unit_cost' ) : 8566.66882845,
# ('INVO XDR-TB', 'unit_cost' ) : 75220.19999999445,
# ('BCG (0-4)', 'unit_cost' ) : 1.9799999999999982,
# ('Symp-Diag', 'unit_cost' ) : 53.55725,
# ('INVO MDR-TB', 'unit_cost' ) : 25073.400000000067,
# ('Mass-Screening' , 'unit_cost' ) : 5072.745234375,
# ('IPT-PLHIV', 'unit_cost' ) : 19.004999999999374,
# ('IPT Gen Pop', 'unit_cost' ) : 18.9555078125,
# ('Pall-Care', 'unit_cost' ) : 2643.0600000000004,
# ('HF MDR-TB', 'unit_cost' ) : 9474.102,
# ('HF DS-TB' , 'unit_cost' ) : 4306.0380000000005,
# ('C-Tracing', 'unit_cost' ) : 812.7678906250001,
# }

# for attribute, val in updated_vals.items():
#     prog = progset.getProg(attribute[0])
#     if attribute[1] == 'budget':
#         continue
#     elif attribute[1] == 'unit_cost':
#         prog.func_specs['pars']['unit_cost'] = val
#     else:
#         continue
#         prog.insertValuePair(t=2017, y=val, attribute=attribute[1], rescale_after_year=True)
#


parset_results = proj.runSim(parset_name=parset_name, plot=False)
options = defaultOptimOptions(settings=proj.settings, progset=proj.progsets[0])
default_results = proj.runSim(parset_name=parset_name, progset_name=progset_name, options=options, plot=False)
#

from optima_tb.reconciliation import reconcile
reconciled_progset,outcome = reconcile(proj,parset_name,progset_name,reconcile_for_year,impact_pars=None,unitcost_sigma=unitcost_sigma, attribute_sigma=attribute_sigma, budget_sigma = budget_sigma)


options = defaultOptimOptions(settings=proj.settings, progset=reconciled_progset)
options['init_alloc'] = {}
reconciled_results = proj.runSim(parset_name=parset_name, progset=reconciled_progset, options=options, plot=False)

#results = proj.runSim(parset_name=name, plot=False)

#%% Extrapolation test
# print proj.parsets[0].getPar('nu_snx').interpolate(tvec=proj.results[0].t_observed_data, pop_label='Miners')
# print list(proj.parsets[0].getPar('nu_snx').interpolate(tvec=proj.results[0].t_observed_data, pop_label='Miners', extrapolate_nan = True))
# oplt.plotResult(proj, reconciled_results, output_labels=['ltlno_rate'], pop_labels=['15-64'])
# oplt.plotResult(proj, default_results, output_labels=['ltlno_rate'], pop_labels=['15-64'])
oplt.plotCompareResults(proj, odict({'parset':parset_results,'default':default_results,'reconciled':reconciled_results}), output_labels=['spxyes_rate'], pop_labels=['15-64'])

# oplt.plotCompareResults(proj, odict({'parset':parset_results,'default':default_results}), output_labels=['spxyes_rate'], pop_labels=['15-64'])

plt.show()

# proj.progsets['default2'] = reconciled_progset
# reconciled_progset = reconcile(proj,parset_name,'default2',reconcile_for_year,impact_pars=None,unitcost_sigma=unitcost_sigma, attribute_sigma=attribute_sigma, budget_sigma = budget_sigma)
