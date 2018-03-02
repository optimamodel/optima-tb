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
    for pop_label, par_label, par_uid in mapping:
        target_vals[par_uid] = parset_results.model.getPop(pop_label).getPar(par_label).vals[-1]

    return target_vals

def updateProgset(pset, attribute_dict, year):
    '''Update the progset with the new values as obtained from reconciliation process, but only for the reconciled year
       Params:
            pset  - ModelProgramSet instance
            new_pars_dict          Dictionary form of the optimized list which is used to update progset
            year    - Year at which to apply the input values
    '''
    progset = pset.progset # Access the traditional values

    for attribute in attribute_dict:
        prog = progset.getProg(attribute[0])
        switch
    for prog_label in progset.prog_ids.keys():
        if prog_label not in new_pars_dict.keys(): continue
        else:
            index = progset.prog_ids[prog_label]
            for attribute in new_pars_dict[prog_label].keys():
                if attribute in progset.progs[index].attributes.keys():
                    progset.progs[index].insertValuePair(t=year, y=new_pars_dict[prog_label][attribute], attribute=attribute, rescale_after_year=True)
                else:
                    continue
            progset.progs[index].func_specs['pars']['unit_cost'] = new_pars_dict[prog_label]['unit_cost']
            progset.progs[index].insertValuePair(t=year, y=new_pars_dict[prog_label]['budget'], attribute='cost', rescale_after_year=True)
    return progset


def objective(attribute_list, pset, tval, dt, target_vals, attribute_dict, constrain_budget, original_budget):
    # tval - is the year we are reconciling at

    # Load the values into the attribute dict
    for idx,val in enumerate(attribute_list):
        attribute_dict[idx] = val

    # Constrain the budget if required
    if constrain_budget:
        normalization = original_budget / sum([val for attrib, val in attribute_dict if attrib[1] == 'budget'])
        for attrib in attribute_dict:
            if attrib[1] == 'budget':
                attribute_dict[attrib] *= normalization

    # Update the programs
    updateProgset(pset.progset, attribute_dict, tval)

    # Compute the new program values
    pset.update_cache(alloc, tval, dt)
    proposed_vals, proposed_coverage = pset.compute_pars(
        0)  # As there is only one timepoint, that is the one we are using

    # Compare them to the target values and compute objective
    obj = 0
    for par_uid in target_vals:
        obj += (target_vals[par_uid] - proposed_vals[par_uid]) ** 2  # Add squared difference in parameter value
        obj += proposed_coverage[par_uid] if proposed_coverage[
                                                 par_uid] > 1 else 0.0  # Add extra penalty for excess coverage
    return obj


# ASD takes in a list of values. So we need to map all of the things we are optimizing onto 
def reconcile(proj, parset_name, progset_name, reconcile_for_year, sigma_dict=None, unitcost_sigma=0.05, attribute_sigma=0.20, budget_sigma=0.0, impact_pars=None, constrain_budget=True, budget_allocation=None, orig_tvec_end=None, max_time=None):
        """
        Recodnciles progset to identified parset, the objective being to match the parameters as closely as possible with identified standar deviation sigma
        
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
            budget_allocation       Dictionary of programs with new budget allocations (type: dict)
            
        Returns:
            progset                 Updated progset with reconciled values
            outcome                 String denoting original and reconciled parset/progset impact comparison
            
        """

        # First, get baseline values
        proj.setYear([2000, reconcile_for_year], False) # This is the easiest way to evaluate the dependent parameter values
        options = defaultOptimOptions(settings=proj.settings, progset=proj.progsets[0])
        parset = proj.parsets[parset_name]
        progset = proj.progsets[progset_name]
        parset_results = proj.runSim(parset=parset, store_results=False) # parset values we want to match
        progset_results = proj.runSim(parset=parset, progset=progset, store_results=False,options=options) # Default results - also, instantiates a ModelProgramSet for use in next steps

        target_vals = _extract_target_vals(parset_results,progset_results,impact_pars)
        attribute_dict = _createAttributeDict(progset_results.model.pset.progset,reconcile_for_year)
        xmin,xmax = _compute_limits(attribute_dict,sigma_dict,unitcost_sigma,attribute_sigma,budget_sigma)

        # Now, make the original attribute dict
        args = {
            'pset': progset_results.model.pset.progset,
            'tval': np.array([reconcile_for_year]),
            'dt': progset_results.model.sim_settings['tvec_dt'],
            'target_vals': target_vals,
            'attribute_dict': dcp(attribute_dict),
            'constrain_budget': constrain_budget,
            'original_budget': sum([val for attrib,val in attribute_dict if attrib[1]=='budget']), # Original budget is sum all all program budget values
            }

        optim_args = {
                     'stepsize': proj.settings.autofit_params['stepsize'],
                     'maxiters': proj.settings.autofit_params['maxiters'],
                     'maxtime': proj.settings.autofit_params['maxtime'],
                     'sinc': proj.settings.autofit_params['sinc'],
                     'sdec': proj.settings.autofit_params['sdec'],
                     'fulloutput': False,
                     'reltol': None
                     }
        if not max_time is None:
            optim_args['maxtime'] = max_time

        best_attribute_list, _, _ = asd(objective, [x for _,x in attribute_dict.items()], args, xmin=[x for _, x in xmin.items()], xmax=[x for _, x in xmax.items()], **optim_args)

        pset = progset_results.model.pset
        original_alloc,tval,_ = pset.get_alloc({'progs_start':reconcile_for_year,'init_alloc':{},'tvec':np.array([reconcile_for_year]),'tvec_dt':progset_results.model.sim_settings['tvec_dt'],'alloc_is_coverage':False,'saturate_with_default_budgets':True})






        impact = {}
        logger.info('Reconciling for year: %i' % reconcile_for_year)
        if impact_pars is None:
            logger.info('No impact pars defined for reconciliation, using all impact parameters')
            impact_pars = progset.impacts.keys()
        else:
            new_pars = [z for z in impact_pars if z in progset.impacts.keys()]
            impact_pars = dcp(new_pars)

        results = proj.runSim(parset_name=parset_name, store_results=False)

        # Get original comparison between progset and parset
        impact['original'] = compareOutcomesFunc(proj=proj, parset_name=parset_name, progset_name=progset_name, year=reconcile_for_year, compareoutcome=True, display=False)


        # Run optimisation
        args = {'proj': proj, 'parset': parset.pars['cascade'], 'progset': progset, 'parset_name': parset_name,
                'impact_pars': impact_pars, 'results': results, 'attribute_dict': attribute_dict,
                'reconcile_for_year': reconcile_for_year, 'compareoutcome': False, 'prog_budget_alloc': budget_allocation, 'constrain_budget': constrain_budget}

        optim_args = {
                     'stepsize': proj.settings.autofit_params['stepsize'],
                     'maxiters': proj.settings.autofit_params['maxiters'],
                     'maxtime': proj.settings.autofit_params['maxtime'],
                     'sinc': proj.settings.autofit_params['sinc'],
                     'sdec': proj.settings.autofit_params['sdec'],
                     'fulloutput': False,
                     'reltol': None
                     }
        if not max_time is None:
            optim_args['maxtime'] = max_time

#         print reconciliationMetric
#         print "att list", attribute_list
#         print "args", args
#         print "xmin", xmin,
#         print "xmax", xmax
#         print "optim args", optim_args
        logging.info("About to use ASD during reconciliation")
        best_attribute_list, _, _ = asd(reconciliationMetric, attribute_list, args, xmin=xmin, xmax=xmax, **optim_args)
        logging.info("Regenerating attributes dictionary")
        best_attribute_dict = regenerateAttributesDict(attribute_list=best_attribute_list, orig_attribute_dict=attribute_dict)
        if constrain_budget:
            best_attribute_dict, _, _, _ = rescaleAllocation(best_attribute_dict, attribute_dict)
        print best_attribute_dict
        logging.info("Updating progset")
        progset = updateProgset(new_pars_dict=best_attribute_dict, progset=progset, year=reconcile_for_year)
        impact['reconciled'] = compareOutcomesFunc(proj=proj, parset_name=parset_name, progset_name=progset_name, year=reconcile_for_year, compareoutcome=True, display=False)
#        print impact['original']['snmno_rate']
#        print impact['reconciled'].keys()
        # Display comparison between old progset and new reconciled progset
        parset_value = 'Parset Impact'
        origprogset_value = 'Original Impact'
        reconcileprogset_value = 'Reconciled Impact'
        print('Comparing outcomes for year: %i' % reconcile_for_year)
        outcome = '\n\t\t\t%s\t\t%s\t\t%s\n' % (parset_value, origprogset_value, reconcileprogset_value)
        for par_label in impact['original'].keys():
            if par_label == 'net_difference': continue
            else:
                outcome += '%s\n' % (par_label)
                for popkey in impact['original'][par_label]:
                    # print('Pop key: %s, Par_label: %s\nType(original parset value): %s\nType(original progset value): %s\nType(reconciled parset value): %s\nType(reconciled progset value): %s\n' %(popkey, par_label, impact['original'][par_label][popkey]['parset_impact_value'], impact['original'][par_label][popkey]['progset_impact_value'], impact['reconciled'][par_label][popkey]['parset_impact_value'], impact['reconciled'][par_label][popkey]['progset_impact_value']))
#                    print impact['original'][par_label][popkey].keys()
#                    print impact['reconciled'][par_label][popkey].keys()
                    outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\n'.format(popkey, impact['original'][par_label][popkey]['parset_impact_value'], impact['original'][par_label][popkey]['progset_impact_capped'], impact['reconciled'][par_label][popkey]['progset_impact_capped'])
                outcome += '\n'
#        print outcome
        # Reset back to original runSim durations
        proj.setYear([2000, orig_tvec_end], False)
        return progset, outcome

def compareOutcomesFunc(proj, year, parset_name=None, progset_name=None, budget_allocation=None, compareoutcome=None, display=True, constrain_budget=True):
    """
    Compares impact parameters as informed by progset to identified parset, and display the comparison
    
    Params:
        proj                    Project object to run simulations for reconciliation process (type: Python object)
        year                    Year for which compaarison needs to be displayed/established (type: int)
        parset_name             Parameter set name to compare  (type: string)
        progset_name            Program set name to compare  (type: string)
        budget_allocation       Dictionary of programs with new budget allocations (type: dict)
        compareoutcome          Flag to pass into reconciliationMetric() to inform that this is not an optimization (type: bool)
        display                 Flag to indicate whether to print out comparison results (type: bool)
        
        
    Returns:
        impact                  dictionary of original and reconciled parset/progset impact comparison
        
    """
    # Make a copy of the original simulation end date
    orig_tvec_end = proj.settings.tvec_end
    # Checks and settings for reconcile
    if parset_name is None:
        try:
            parset_name = proj.parsets.keys()[0]
            logger.info('Parameter set was not identified for impact parameter comparison, using parameter set: "%s"' % parset_name)
        except: raise OptimaException('No valid parameter sets exist within the project')

    if progset_name is None:
        try:
            progset_name = proj.progsets.keys()[0]
            logger.info('Program set was not identified for impact parameter comparison, using program set: "%s"' % progset_name)
        except: raise OptimaException('No valid program sets exist within the project')

    if not parset_name in proj.parsets.keys(): raise OptimaException("ERROR: No parameter set '%s' found" % parset_name)
    if not progset_name in proj.progsets.keys(): raise OptimaException("ERROR: No program set '%s' found" % progset_name)

    # Set years for Simulation runs
    proj.setYear([2000, year], False)
    # Setup compareOutcomes kwargs
    parset = proj.parsets[parset_name].pars['cascade']
    progset = proj.progsets[progset_name]
    results = proj.runSim(parset_name=parset_name, store_results=False)
    # Declare impact parameters to use for display (use all)
    impact_pars = progset.impacts.keys()
    # use reconcilemetric to get desired result
    impact = reconciliationMetric(new_attributes=[], proj=proj, parset=parset, progset=progset,
                                  parset_name=parset_name, impact_pars=impact_pars,
                                  results=results, attribute_dict={}, reconcile_for_year=year,
                                  compareoutcome=compareoutcome, prog_budget_alloc=budget_allocation, constrain_budget=constrain_budget)
    # display output
    if display:
        print('Comparing outcomes for year: %i' % year)
        parset_value = 'parset_impact_value'
        progset_value_uncapped = 'progset_impact_uncapped'
        progset_value_capped = 'progset_impact_capped'
        progset_overflow = 'overflow_list'
        outcome = '\n\t\t\t%s\t%s\t%s\t%s\n' % (parset_value, progset_value_uncapped, progset_value_capped, progset_overflow)
        for par_label in impact.keys():
            if par_label == 'net_difference': continue
            else:
                outcome += '%s\n' % (par_label)
                for popkey in impact[par_label]:
                    try: outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\t\t{:<10}\n'.format(popkey, impact[par_label][popkey][parset_value], impact[par_label][popkey][progset_value_uncapped], impact[par_label][popkey][progset_value_capped], impact[par_label][popkey][progset_overflow])
                    except:
                        try: outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\t\t{:<10}\n'.format(popkey, impact[par_label][popkey][parset_value], impact[par_label][popkey][progset_value_uncapped][0], impact[par_label][popkey][progset_value_capped][0], impact[par_label][popkey][progset_overflow])
                        except:
                            try: outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\t\t{:<10}\n'.format(popkey, impact[par_label][popkey][parset_value], impact[par_label][popkey][progset_value_uncapped][0], impact[par_label][popkey][progset_value_capped], impact[par_label][popkey][progset_overflow])
                            except:
                                outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\t\t{:<10}\n'.format(popkey, impact[par_label][popkey][parset_value], impact[par_label][popkey][progset_value_uncapped], impact[par_label][popkey][progset_value_capped][0], impact[par_label][popkey][progset_overflow])
                outcome += '\n'
        print outcome
    # Reset back to original runSim durations
    proj.setYear([2000, orig_tvec_end], False)
    return impact

def _createAttributeDict(progset,reconcile_for_year):
    '''Creates an attribute dictionary on a per program basis from the identified progset 
       for all parameters/impact labels that can be reconciled
       
       Params:
            settings                Project settings to identify year range for interpolation of progset
            progset                 progset on which interpolation needs to be conducted
            
        Returns:
            attributes_dict         dict mapping (prog_label,attribute):value
            
    '''

    attributes_dict = odict()
    tval = np.array([reconcile_for_year])
    for prog in progset.progs:
        if 'unit_cost' in prog.func_specs['pars'].keys() and prog.getDefaultBudget(year=tval) > 0.:
            attributes_dict[(prog.label,'unit_cost')] = prog.func_specs['pars']['unit_cost']
            attributes_dict[(prog.label,'budget')] = prog.getDefaultBudget(year=tval)
            interpolated_attributes = prog.interpolate(tvec=tval)
            for key,val in interpolated_attributes.items():
                if key in ['cov' , 'dur' , 'time', 'cost']:
                    continue
                else:
                    attributes_dict[(prog.label,key)] = val[0]
    return attributes_dict

def _compute_limits(attribute_dict, sigma_dict, unitcost_sigma, attribute_sigma, budget_sigma):
    xmin = odict()
    xmax = odict()

    for attribute, val in attribute_dict.items():
        # Attribute is e.g. ('HF DS-TB', 'budget')
        if sigma_dict is not None and attribute in sigma_dict:
            xmin[attribute] = val*(1 - sigma_dict[attribute])
            xmax[attribute] = val*(1 + sigma_dict[attribute])
        elif attribute[1] == 'unit_cost':
            xmin[attribute] = val*(1 - unitcost_sigma)
            xmax[attribute] = val*(1 + unitcost_sigma)
        elif attribute[1] == 'budget':
            xmin[attribute] = val*(1 - budget_sigma)
            xmax[attribute] = val*(1 + budget_sigma)
        else:
            xmin[attribute] = val*(1 - attribute_sigma)
            xmax[attribute] = val*(1 + attribute_sigma)
    return xmin, xmax

