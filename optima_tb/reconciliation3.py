from optima_tb.utils import odict, OptimaException
import optima_tb.settings as settings
from copy import deepcopy as dcp
from optima_tb.asd import asd
import numpy as np
from optima_tb.parsing import FunctionParser
import optima_tb.settings as project_settings
from optima_tb.defaults import defaultOptimOptions

import logging
logger = logging.getLogger(__name__)

def _extract_target_vals(parset_results, progset_results, impact_pars):
    # Store the target parset values in the same form as the progset computed parameter values
    # This will also mean that target_vals contains the par UIDs only for the impact parameters
    # requested for reconcilation (if None, then all will be used)
    if impact_pars is not None:
        prog_par_uids = [par.uid for par in progset_results.model.pset.pars if par.label in impact_pars]
    else:
        prog_par_uids = [par.uid for par in progset_results.model.pset.pars]

    # Then, get the pop_label/parameter_label combination
    mapping = []
    for pop in progset_results.model.pops:
        for par in pop.pars:
            if par.uid in prog_par_uids:
                mapping.append((pop.label, par.label, par.uid))

    # Finally, extract the values from the parset_results for the corresponding labels and store them in a dict
    # using the UIDs of the progset_results run
    target_vals = dict()  # This is a dict mapping par_uid:par_value for the reconciliation year in the Parset - the UIDs correspond to the progset run though
    initial_prog_vals = dict()
    for pop_label, par_label, par_uid in mapping:
        target_vals[par_uid] = parset_results.model.getPop(pop_label).getPar(par_label).vals
        initial_prog_vals[par_uid] = progset_results.model.getPop(pop_label).getPar(par_label).vals

    return target_vals,initial_prog_vals


def _update_progset(progset, attribute_list, attribute_dict, tval):
    # Update the ModelProgramSet/ProgramSet in place and update the cache ready for parameter computation

    # Load the attribute list into the attribute dict
    for idx, val in enumerate(attribute_list):
        attribute_dict[idx] = val

    # Update the programs and update the new allocation
    progset = dcp(progset)
    for attribute, val in attribute_dict.items():
        prog = progset.getProg(attribute[0])

        if attribute[1] == 'unit_cost':
            prog.func_specs['pars']['unit_cost'] *= val
        else:
            prog.insertValuePair(t=tval[0], y=val, attribute=attribute[1], rescale_after_year=True) # TODO - Multiply
    return progset

def _createAttributeDict(progset):
    '''Creates an attribute dictionary on a per program basis from the identified progset
       for all parameters/impact labels that can be reconciled

       Params:
            settings                Project settings to identify year range for interpolation of progset
            progset                 progset on which interpolation needs to be conducted

        Returns:
            attributes_dict         dict mapping (prog_label,attribute):value

    '''

    attributes_dict = odict()

    for prog in progset.progs:
        if 'unit_cost' in prog.func_specs['pars'].keys():
            attributes_dict[(prog.label, 'unit_cost')] = 1.0
            interpolated_attributes = prog.interpolate(tvec=tval)
            for key, val in interpolated_attributes.items():
                if key not in ['cov', 'dur', 'time', 'cost']:
                    attributes_dict[(prog.label, key)] = 1.0
    return attributes_dict

def _prepare_arrays(attribute_dict, sigma_dict, unitcost_sigma, attribute_sigma):
    xmin_d = odict()
    xmax_d = odict()
    attribute_dict = dcp(attribute_dict)

    for attribute, val in attribute_dict.items():
        # Attribute is e.g. ('HF DS-TB', 'budget')
        if sigma_dict is not None and attribute in sigma_dict:
            xmin_d[attribute] = (1 - sigma_dict[attribute])
            xmax_d[attribute] = (1 + sigma_dict[attribute])
        elif attribute[1] == 'unit_cost':
            xmin_d[attribute] = (1 - unitcost_sigma)
            xmax_d[attribute] = (1 + unitcost_sigma)
        else:
            xmin_d[attribute] = (1 - attribute_sigma)
            xmax_d[attribute] = (1 + attribute_sigma)

    for attribute in attribute_dict.keys():
        if xmin_d[attribute] == xmax_d[attribute]: # No allowed change means cannot optimize this quantity
            attribute_dict.pop(attribute)
            xmin_d.pop(attribute)
            xmax_d.pop(attribute)

    x0 = [x for _,x in attribute_dict.items()]
    xmin = np.array([x for _, x in xmin_d.items()])
    xmax = np.array([x for _, x in xmax_d.items()])

    return x0, xmin, xmax, attribute_dict

def _objective(attribute_list, pset, tval,attribute_dict,parset):
    new_pset = _update_progset(pset.progset, attribute_list, attribute_dict)
    pset.update_cache(alloc, tval, dt)

    args = {
        'pset': progset_results.model.pset,  # This is the ProgramSet to be optimized
        'tval': np.array([reconcile_for_year]),
        'attribute_dict': dcp(attribute_dict),
        'parset': parset,
    }


    proposed_vals, proposed_coverage = pset.compute_pars(0)  # As there is only one timepoint, that is the one we are using

    'tvec':parset_results.model.t


    obj = 0.0
    for par_uid in target_vals:
        obj += (target_vals[par_uid] - proposed_vals[par_uid]) ** 2  # Add squared difference in parameter value
        obj += sum([x for x in proposed_coverage[par_uid] if x > 1])  # Add extra penalty for excess coverage

    return obj

# ASD takes in a list of values. So we need to map all of the things we are optimizing onto 
def reconcile(proj, parset_name, progset_name, sigma_dict=None, unitcost_sigma=0.05, attribute_sigma=0.20, impact_pars=None, constrain_budget=False, max_time=15):
        """
        Reconciles progset to identified parset, the objective being to match the parameters as closely as possible with identified standar deviation sigma
        
        Params:
            proj                    Project object to run simulations for reconciliation process (type: Python object)
            reconcile_for_year      Year for which reconciliation needs to be done to ensure continuity of parset/progset (type: int)
            parset_name             Parameter set name to match/reconcile  (type: string)
            progset_name            Program set name to match/reconcile  (type: string)
            sigma_dict              Dictionary that links program labels with sigmas for 'unit_cost', 'budget' and 'attribute' (type: dict)
                                    If not provided, sigmas default across all programs to the ones provided by the other kwargs
                                    For example, sigma_dict = {('HF DS-TB','unit_cost'):0.4,('HF MDR-TB','budget'):1000}
            unitcost_sigma          Standard deviation allowable for Unit Cost (type: float)
            attribute_sigma         Standard deviation allowable for attributes identified in impact_pars (type: float)
            budget_sigma            Standard deviation allowable for budget for program (type: float)
            impact_pars             Impact pars to be reconciled (type: list or None)
            constrain_budget        Flag to inform algorithm whether to constrain total budget or not (type: bool)

        Returns:
            progset                 Updated progset with reconciled values
            outcome                 String denoting original and reconciled parset/progset impact comparison
            
        """

        # First, get baseline values
        options = defaultOptimOptions(settings=proj.settings, progset=proj.progsets[0])
        parset = proj.parsets[parset_name]
        parset_results = proj.runSim(parset=parset, store_results=False) # parset values we want to match
        progset_results = proj.runSim(parset=parset, progset_name=progset_name, store_results=False,options=options) # Default results - also, instantiates a ModelProgramSet for use in next steps

        # Extract the target values, make an attribute dict, etc.
        target_vals, initial_prog_vals = _extract_target_vals(parset_results,progset_results,impact_pars)
        attribute_dict = _createAttributeDict(progset_results.model.pset.progset,reconcile_for_year)
        x0,xmin,xmax,attribute_dict = _prepare_arrays(attribute_dict,sigma_dict,unitcost_sigma,attribute_sigma,budget_sigma)

        # Now, make the original attribute dict
        pset = progset_results.model.pset
        original_alloc = pset.get_alloc(progset_results.model.t,progset_results.model.dt,{'progs_start':reconcile_for_year,'init_alloc':{},'alloc_is_coverage':False,'saturate_with_default_budgets':True})

        # Copy over the popsizes from the original simulation
        for pop in progset_results.model.pops:
            for comp in pop.comps:
                comp.vals = parset_results.model.getPop(pop.label).getComp(comp.label).vals[-1:]

        args = {
            'pset': progset_results.model.pset, # This is the ProgramSet to be optimized
            'attribute_dict': dcp(attribute_dict),
            'parset':parset,
            'tvec':parset_results.model.t
            }

        optim_args = {
                     'stepsize': proj.settings.autofit_params['stepsize'],
                     'maxiters': proj.settings.autofit_params['maxiters'],
                     'maxtime': proj.settings.autofit_params['maxtime'],
                     'sinc': proj.settings.autofit_params['sinc'],
                     'sdec': proj.settings.autofit_params['sdec'],
                     'fulloutput': False,
                     'reltol': None,
                     'xmin': xmin,
                     'xmax': xmax,
                     }

        if not max_time is None:
            optim_args['maxtime'] = max_time

        best_attribute_list, _, _ = asd(_objective, x0, args, **optim_args)

        # Now, make the original attribute dict
        alloc = _update_progset(pset.progset, best_attribute_list, args['attribute_dict'], args['constrain_budget'], args['original_budget'], original_alloc, args['tval'], args['dt'])
        pset.update_cache(alloc, args['tval'], args['dt'])
        proposed_vals, proposed_coverage = pset.compute_pars(0)

        outcome = '%s %s %s %s %s\n' % ('Population'.ljust(15),'Parameter'.ljust(15),'Parset'.rjust(10), '    Initial'.ljust(21),'      Reconciled'.ljust(21))
        residual_old = 0.0
        residual_new = 0.0
        for pop in progset_results.model.pops:
            for par in pop.pars:
                if par.uid in target_vals:
                    outcome += '%s %s %10.4f %10.4f (%+10.4f) %10.4f (%+10.4f)\n' % (pop.label.ljust(15),par.label.ljust(15),target_vals[par.uid],initial_prog_vals[par.uid],initial_prog_vals[par.uid]-target_vals[par.uid],proposed_vals[par.uid],proposed_vals[par.uid]-target_vals[par.uid])
                    residual_old += (initial_prog_vals[par.uid]-target_vals[par.uid])**2
                    residual_new += (proposed_vals[par.uid]-target_vals[par.uid])**2
        outcome += 'Old residual (no coverage penalty) = %10.4f, New residual = %10.4f\n' % (residual_old,residual_new)

        return pset.progset, outcome


