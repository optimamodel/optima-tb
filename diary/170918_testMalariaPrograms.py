from optima_tb.project import Project
from optima_tb.utils import odict
from optima_tb.utils import OptimaException
import os
import pylab as pl
import numpy as np
import csv
from optima_tb.asd import asd



cascade_path = os.path.abspath('.//cascade-m-170918.xlsx')
databook_path = os.path.abspath('.//databook-m-170918.xlsx')

sim_step = 365.
sim_start = -3000.
sim_stop = 3825.

np.seterr(all='raise')



# Function for inspecting result outputs
def plotResults(res_dict, dt, t_start, t_end, labels, den_labels = None, pop_labels = None):
    
    figr, axr = pl.subplots()

    tvec = np.arange(t_start, t_end, dt)
    
    for result_label in res_dict.keys():
        result = res_dict[result_label]
        for label in labels:
            vals = result.getValuesAt(label, t_start, t_end, pop_labels=pop_labels)[0]
            axr.plot(tvec, vals, '-', label=result_label)
    
    axr.legend()
    
    return figr


# Forces the result of a function evaluation to be interpreted as number not as fraction
def interpret_value_as_number(proj, config='default'):
    # force the transitions from the birth compartments to be interpreted as 'number', not as 'fraction'
    for pop in proj.parsets[config].pop_labels:
        for label in ['num_hb', 'num_mb', 'num_htx']:
            idx = proj.parsets[config].par_ids['cascade'][label]
            proj.parsets[config].pars['cascade'][idx].y_format[pop] = 'number'
    return proj


def setup_project():
    proj = Project(name='malaria', cascade_path=cascade_path, validation_level='warn')

    proj.loadSpreadsheet(databook_path=databook_path)

    dt = proj.settings.tvec_dt = 1. / 1.
    t_start = proj.settings.tvec_start = sim_start
    proj.settings.tvec_observed_end = sim_stop
    t_end = proj.settings.tvec_end = sim_stop

    print('Simulating from ' + str(t_start) + ' to ' + str(t_end) + ' in ' + str(dt) + ' increments')

    proj.makeParset()
    proj.makeProgset()

    proj = interpret_value_as_number(proj)

    return proj


# wrapper function for plotResults
def visualize(proj, res_dict):
    outputs = ['t_hinf']

    start = 2000
    stop = sim_stop

    for pop_label in proj.data['pops']['label_names'].keys():
        for output in outputs:
            fig = plotResults(res_dict=res_dict, dt=proj.settings.tvec_dt, t_start=start, t_end=stop, labels=[output],
                              pop_labels=[pop_label])
            fig.suptitle('Impact of all programs on ' + pop_label)

    pl.show()


if __name__ == "__main__":
    proj = setup_project()
    proj = interpret_value_as_number(proj, 'default')


    # define various program combinations
    opt_all = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'TXg': 1.0, 'TXp': 1.0, 'TXc': 1.0,
                              'BCCg': 0.99, 'BCCp': 0.99, 'BCCc': 0.99,
                              'IPTp': 0.61, 'IRS': 0.036, 'LAV': 0.05,
                              'LLINg': 0.201, 'LLINp': 0.655, 'LLINc': 0.63,
                              'MDA': 0.064, 'SMC': 0.05}}

    opt_TX = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'TXg': 1.0, 'TXp': 1.0, 'TXc': 1.0}}
    opt_BCC = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'BCCg': 0.99, 'BCCp': 0.99, 'BCCc': 0.99}}
    opt_IPTp = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'IPTp': 0.61}}
    opt_IRS = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'IRS': 0.036}}
    opt_LAV = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'LAV': 0.05}}
    # LLINp is intentionally missing to validate that even though it is missing the influence of LLINg and LLINc is also applied to the preg population
    opt_LLIN = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'LLINg': 0.201, 'LLINc': 0.63}}
    opt_MDA = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'MDA': 0.064}}
    opt_SMC = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'SMC': 0.05}}

    opt_BCC_LLIN = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'BCCg': 0.99, 'BCCp': 0.99, 'BCCc': 0.99,
                              'LLINg': 0.201, 'LLINp': 0.655, 'LLINc': 0.63}}
    opt_BCC_IPTp = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'BCCg': 0.99, 'BCCp': 0.99, 'BCCc': 0.99, 'IPTp': 0.61}}
    opt_BCC_IPTp_LLIN = {'progs_start': 2000., 'progs_end': sim_stop, 'alloc_is_coverage': True,
               'init_alloc': {'BCCg': 0.99, 'BCCp': 0.99, 'BCCc': 0.99,
                              'IPTp': 0.61, 'LLINg': 0.201, 'LLINp': 0.655, 'LLINc': 0.63}}

    res_dict = odict()
    res_dict['no_progs'] = proj.runSim(parset_name='default')
    res_dict['all_progs'] = proj.runSim(parset_name='default', progset_name='default', options=opt_all)
    # res_dict['TX'] = proj.runSim(parset_name='default', progset_name='default', options=opt_TX)
    # res_dict['BCC'] = proj.runSim(parset_name='default', progset_name='default', options=opt_BCC) # will not have any effect
    # res_dict['IPTp'] = proj.runSim(parset_name='default', progset_name='default', options=opt_IPTp)
    # res_dict['IRS'] = proj.runSim(parset_name='default', progset_name='default', options=opt_IRS)
    # res_dict['LAV'] = proj.runSim(parset_name='default', progset_name='default', options=opt_LAV)
    # res_dict['LLIN'] = proj.runSim(parset_name='default', progset_name='default', options=opt_LLIN)
    # res_dict['MDA'] = proj.runSim(parset_name='default', progset_name='default', options=opt_MDA)
    # res_dict['SMC'] = proj.runSim(parset_name='default', progset_name='default', options=opt_SMC)
    # res_dict['BCC_LLIN'] = proj.runSim(parset_name='default', progset_name='default', options=opt_BCC_LLIN)
    # res_dict['BCC_IPTp'] = proj.runSim(parset_name='default', progset_name='default', options=opt_BCC_IPTp)
    # res_dict['BCC_IPTp_LLIN'] = proj.runSim(parset_name='default', progset_name='default', options=opt_BCC_IPTp_LLIN)

    visualize(proj, res_dict)
