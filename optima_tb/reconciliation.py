from optima_tb.utils import odict, OptimaException
import optima_tb.settings as settings
from copy import deepcopy as dcp
from optima_tb.asd import asd
import numpy as np
from optima_tb.parsing import FunctionParser
from optima_tb.model import *

import logging
logger = logging.getLogger(__name__)




def reconcileFunc(proj, reconcile_for_year, parset_name, progset_name, sigma_dict = None, unitcost_sigma = 0.05, attribute_sigma = 0.20, budget_sigma = 0.0, impact_pars = None, constrain_budget=True, budget_allocation = None, orig_tvec_end = None, max_time = None):
        """
        Reconciles progset to identified parset, the objective being to match the parameters as closely as possible with identified standard deviation sigma
        
        Params:
            proj                    Project object to run simulations for reconciliation process (type: Python object)
            reconcile_for_year      Year for which reconciliation needs to be done to ensure continuity of parset/progset (type: int)
            parset_name             Parameter set name to match/reconcile  (type: string)
            progset_name            Program set name to match/reconcile  (type: string)
            sigma_dict              Dictionary that links program labels with sigmas for 'unit_cost', 'budget' and 'attribute' (type: dict)
                                    If not provided, sigmas default across all programs to the ones provided by the other kwargs
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

        #Set years for Simulation runs
        proj.setYear([2000, reconcile_for_year], False)
        
        #Setup parameters for reconciliation
        parset  = proj.parsets[parset_name].pars['cascade']
        progset = proj.progsets[progset_name]
        impact = {}
        logger.info('Reconciling for year: %i' %reconcile_for_year)
        if impact_pars is None: 
            logger.info('No impact pars defined for reconciliation, using all impact parameters')
            impact_pars = progset.impacts.keys()
        else:
            new_pars = [z for z in impact_pars if z in progset.impacts.keys()]
            impact_pars = dcp(new_pars)
        
        results = proj.runSim(parset_name=parset_name, store_results=False)
        
        #Get original comparison between progset and parset
        impact['original'] = compareOutcomesFunc(proj=proj, parset_name=parset_name, progset_name=progset_name, year=reconcile_for_year, compareoutcome=True, display=False)
        
        #Convert into an optimisable list
        attribute_dict = createAttributeDict(settings=proj.settings, progset=progset)
#        attribute_list, unitcost_index, budget_index = createAttributeList(attribute_dict=attribute_dict)
        attribute_list, index_dict = createAttributeList(attribute_dict=attribute_dict)
        #Setup min-max bounds for optimisation
        xmin, xmax = dcp(attribute_list), dcp(attribute_list)
        
        #Setup xmin and xmax conditions, unit_cost indexes are used in case unit_costs have a difference standard deviation value
        for index in range(len(attribute_list)):
            par_type = index_dict[index]['type']
            if sigma_dict is None:
                if par_type == 'unit_cost':
                    xmin[index] *= (1-unitcost_sigma)
                    xmax[index] *= (1+unitcost_sigma)
                elif par_type == 'budget':
                    xmin[index] *= (1-budget_sigma)
                    xmax[index] *= (1+budget_sigma)
                else:
                    xmin[index] *= (1-attribute_sigma)
                    xmax[index] *= (1+attribute_sigma)
                    if xmax[index] > 1.: xmax[index] = 1.   # TODO: Work out if this restriction is necessary here.
            else:
                prog_label = index_dict[index]['prog_label']
                if not par_type in ['unit_cost','budget']:
                    par_type = 'attribute'
                try: sigma = sigma_dict[prog_label][par_type]
                except Exception as E: raise OptimaException('ERROR: Sigma dictionary has been provided to program reconciliation process in incomplete form.\nDETAILS: %s' % repr(E))
                xmin[index] *= (1-sigma)
                xmax[index] *= (1+sigma)
                
            
            if xmin[index] <= 0.: xmin[index] = settings.TOLERANCE
        
        #Run optimisation
        args = {'proj': proj, 'parset': parset, 'progset': progset, 'parset_name': parset_name,
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

        # TODO change optimisation method
        best_attribute_list = asd(reconciliationMetric, attribute_list, args, xmin = xmin, xmax = xmax, **optim_args)
        best_attribute_dict = regenerateAttributesDict(attribute_list = best_attribute_list, orig_attribute_dict = attribute_dict)
        print best_attribute_dict
        progset = updateProgset(new_pars_dict=best_attribute_dict, progset=progset, year=reconcile_for_year)
        impact['reconciled'] = compareOutcomesFunc(proj=proj, parset_name=parset_name, progset_name=progset_name, year=reconcile_for_year, compareoutcome=True, display=False)
#        print impact['original']['snmno_rate']
#        print impact['reconciled'].keys()
        #Display comparison between old progset and new reconciled progset
        parset_value = 'Parset Impact'
        origprogset_value = 'Original Impact'
        reconcileprogset_value = 'Reconciled Impact'
        print('Comparing outcomes for year: %i' %reconcile_for_year)
        outcome = '\n\t\t\t%s\t\t%s\t\t%s\n' %(parset_value, origprogset_value, reconcileprogset_value)

        for par_label in impact['original'].keys():
            if par_label == 'net_difference': continue
            else:
                outcome += '%s\n' %(par_label)
                for popkey in impact['original'][par_label]:
                    outcome += '\t%s' % popkey
                    outcome += '\t%10.2f\t' % impact['original'][par_label][popkey]['parset_impact_value']
                    outcome += '\t%10.2f\t' % impact['original'][par_label][popkey]['progset_impact_capped']
                    outcome += '\t%10.2f\n' % impact['reconciled'][par_label][popkey]['progset_impact_capped']
                outcome += '\n'
#        print outcome
        #Reset back to original runSim durations
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
    #Make a copy of the original simulation end date
    orig_tvec_end = proj.settings.tvec_end
    #Checks and settings for reconcile
    if parset_name is None: 
        try: 
            parset_name = proj.parsets.keys()[0]
            logger.info('Parameter set was not identified for impact parameter comparison, using parameter set: "%s"' %parset_name)
        except: raise OptimaException('No valid parameter sets exist within the project')
        
    if progset_name is None: 
        try:
            progset_name = proj.progsets.keys()[0]
            logger.info('Program set was not identified for impact parameter comparison, using program set: "%s"' %progset_name)
        except: raise OptimaException('No valid program sets exist within the project')
    
    if not parset_name in proj.parsets.keys(): raise OptimaException("ERROR: No parameter set '%s' found"%parset_name)
    if not progset_name in proj.progsets.keys(): raise OptimaException("ERROR: No program set '%s' found"%progset_name)
    
    #Set years for Simulation runs
    proj.setYear([2000, year], False)
    #Setup compareOutcomes kwargs
    parset  = proj.parsets[parset_name].pars['cascade']
    progset = proj.progsets[progset_name]
    results = proj.runSim(parset_name=parset_name, store_results=False)
    #Declare impact parameters to use for display (use all)
    impact_pars = progset.impacts.keys()
    #use reconcilemetric to get desired result
    impact = reconciliationMetric(new_attributes=[], proj=proj, parset=parset, progset=progset, 
                                  parset_name=parset_name, impact_pars=impact_pars, 
                                  results=results, attribute_dict={}, reconcile_for_year=year, 
                                  compareoutcome=compareoutcome, prog_budget_alloc=budget_allocation, constrain_budget=constrain_budget)
    #display output
    if display:
        print('Comparing outcomes for year: %i' %year)
        parset_value  = 'parset_impact_value'
        progset_value_uncapped = 'progset_impact_uncapped'
        progset_value_capped = 'progset_impact_capped'
        progset_overflow = 'overflow_list'
        outcome = '\n\t\t\t%s\t%s\t%s\t%s\n' %(parset_value, progset_value_uncapped, progset_value_capped, progset_overflow)
        for par_label in impact.keys():
            if par_label == 'net_difference': continue
            else:
                outcome += '%s\n' %(par_label)
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
    #Reset back to original runSim durations
    proj.setYear([2000, orig_tvec_end], False)
    return impact

def createAttributeDict(settings, progset):
    '''Creates an attribute dictionary on a per program basis from the identified progset 
       for all parameters/impact labels that can be reconciled
       
       Params:
            settings                Project settings to identify year range for interpolation of progset
            progset                 progset on which interpolation needs to be conducted
            
        Returns:
            attributes_dict         A dictionary of all program_labels and respective attributes that are to be reconciled
    '''
    attributes_dict = odict()
    for prog_label in progset.prog_ids.keys():
        #print('Program Name: %s\n' %prog_label)
        mark_for_delete = False
        if prog_label not in attributes_dict.keys(): attributes_dict[prog_label] = odict()
        index = progset.prog_ids[prog_label]
        try:    
            attributes_dict[prog_label]['unit_cost'] = progset.progs[index].func_specs['pars']['unit_cost']
            #print('Unit Cost: %g\n' %attributes_dict[prog_label]['unit_cost'])
            attributes_dict[prog_label]['budget'] = progset.progs[index].getDefaultBudget(year=settings.tvec_end)
            #print('Budget: %g\n' %attributes_dict[prog_label]['budget'])
            interpolated_attributes = progset.progs[index].interpolate(tvec=np.arange(settings.tvec_start, settings.tvec_end + settings.tvec_dt/2, settings.tvec_dt))
            #TODO : generalise key
            for key in interpolated_attributes:
                if key == 'cov' or key == 'dur' or key == 'time': continue
                elif interpolated_attributes[key][-1] <= 0: 
                    mark_for_delete = True
                    break
                elif key == 'cost': continue
                else: attributes_dict[prog_label][key] = interpolated_attributes[key][-1]
        except: mark_for_delete = True
        
        if mark_for_delete:
            del attributes_dict[prog_label]
    return attributes_dict

def createAttributeList(attribute_dict):
    '''Converts the attribute dictionary into a list so that it can be passed into the asd function for reconciliation/optimization
       
       Params:
            attribute_dict          A dictionary of all program_labels and respective attributes that are to be reconciled
            
        Returns:
            attribute_list          A list form of the attributes dictionary
            index_dict              Maps indices of the attribute_list to program labels by key 'prog_label' and attribute type by key 'type'
    '''
    attribute_list = []
    index_dict = {}
#    unitcost_index = []
#    budget_index = []
    index = 0
    for prog_label in attribute_dict:
        for par in attribute_dict[prog_label]:
            attribute_list.append(attribute_dict[prog_label][par])
            index_dict[index] = {'prog_label':prog_label,'type':par}
#            if par == 'unit_cost':
#                unitcost_index.append(index)
#                index += 1
#            elif par == 'budget':
#                budget_index.append(index)            
#                index += 1
#            else: index += 1
            index += 1
    return attribute_list, index_dict #unitcost_index, budget_index


# special behaviour of program with the tag 'supp'
def processSuppTagReconcile(par_label, prog, progset, budget_alloc, ti):
    impacts = []
    parser = FunctionParser(debug=False)
    # extract all attributes which refer to another program
    refs = filter(lambda x: x.startswith('$ref_'), prog.attributes)
    # list of all attributes other than programs
    var = list(set(prog.attributes.keys()).difference(set(refs)))

    # loop over all referenced programs and increase their impact on impact parameters
    for p in refs:
        # label of referenced program
        ref_prog = prog.attributes[p][0]
        # suffix of the program, everything after the last '_' in program label (incl. '_')
        suff = p[p.rfind('_'):]

        if ref_prog in progset.impacts[par_label] and ref_prog in budget_alloc.keys():
            # all parameters specified in the spreadsheet are multiplied
            coeff = 1.
            for f in filter(lambda x: x.endswith(suff), var):
                try: coeff *= prog.attributes[f][ti]
                except: coeff *= prog.attributes[f]

            # obtain coverage of referenced program ..
            cov = progset.getProg(ref_prog).getCoverage(budget_alloc[ref_prog])
            scov = progset.getProg(prog.label).getCoverage(budget_alloc[prog.label])

            # .. and determine the minimum coverage for the impact
            net_cov = min(cov, scov)

            # apply impact
            impacts.append(progset.getProg(ref_prog).getImpact(budget_alloc[ref_prog], impact_label=par_label, parser=parser, years=[ti]) * net_cov * coeff)

    return np.sum(impacts)


# special behaviour for programs with the tag 'scale_props'
def processScalePropsTagReconcile(par_label, pops, pop_ids, progset, budget_alloc, parser, ti):
    # determine total population
    total_pop = sum([p.getDep('h_alive').vals[ti] for p in pops])
    impacts = []  # impacts are treated as fractions!
    # TODO allow impact numbers, too

    # find all programs of the same type as current one and weigh its impact proportional
    # to the population size before summing up
    rel_progs = progset.impacts[par_label] # all programs impacting considered parameter
    rel_progs = filter(lambda x: x in progset.GPs + progset.SOPs, rel_progs) # remove programs which are not GP or SOP

    for p in rel_progs:
        imp = 0.0
        target_pop_size = np.sum([pops[pop_ids[pp]].getDep('h_alive').vals[ti] for pp in progset.getProg(p).target_pops])
        if p in progset.GPs:
            progset.getProg(p).getImpact(budget_alloc[p], impact_label=par_label, parser=parser, years=[ti])
        elif p in progset.SOPs:
            imp = processSuppTagReconcile(par_label, progset.getProg(p), progset, budget_alloc, ti)

        impacts.append(imp * target_pop_size / total_pop)

    return np.sum(impacts)


def regenerateAttributesDict(attribute_list, orig_attribute_dict):
    '''Reverse process, where the attributes list is converted back into the attributes dictionary after optimization/reconciliation
       
       Params:
            attribute_list              A list of reconciled/optimized parameters that are to be used to update progset with new values
            orig_attribute_dict         Reference dictionary to map the list back onto the dictionary
            
       Returns:
            attribute_dict              Dictionary form of the optimized list which is used to update progset
    '''
    attribute_dict = dcp(orig_attribute_dict)
    index = 0
    for prog_label in attribute_dict.keys():
        for par in attribute_dict[prog_label]:
            attribute_dict[prog_label][par] = attribute_list[index]
            index += 1
    return attribute_dict

def updateProgset(new_pars_dict, progset, year):
    '''Update the progset with the new values as obtained from reconciliation process, but only for the reconciled year
    
       Params:
            new_pars_dict          Dictionary form of the optimized list which is used to update progset
            progset                Program set that is to be updated with the new values
            
       Returns:
            attribute_dict         Dictionary form of the optimized list which is used to update progset
    '''
    for prog_label in progset.prog_ids.keys():
        if prog_label not in new_pars_dict.keys(): continue
        else:
            # print prog_label
            index = progset.prog_ids[prog_label]
            for attribute in new_pars_dict[prog_label].keys():
                if attribute in progset.progs[index].attributes.keys():
#                    progset.progs[index].attributes[attribute][-1] = new_pars_dict[prog_label][attribute]
#                    print year
#                    print new_pars_dict[prog_label]['budget']
#                    print attribute
                    progset.progs[index].insertValuePair(t=year, y=new_pars_dict[prog_label][attribute], attribute=attribute, rescale_after_year=True)
                else:
                    continue
            progset.progs[index].func_specs['pars']['unit_cost'] = new_pars_dict[prog_label]['unit_cost']
#            progset.progs[index].cost[-1] = new_pars_dict[prog_label]['budget']
#            print year
#            print new_pars_dict[prog_label]['budget']
#            print 'cost'
            progset.progs[index].insertValuePair(t=year, y=new_pars_dict[prog_label]['budget'], attribute='cost', rescale_after_year=True)
    return progset

def rescaleAllocation(proposed_dict, original_dict):
    '''This function normalises the proposed allocated budget to make sure that the total allowable budget is maintained
    '''
    proposed_total_budget = 0.
    orig_total_budget = 0.
    constrained_total_budget = 0.
    
    for prog_name in proposed_dict.keys():
        orig_total_budget += original_dict[prog_name]['budget']
        proposed_total_budget += proposed_dict[prog_name]['budget']
    
    for prog_name in proposed_dict.keys():
#        scale_factor = proposed_dict[prog_name]['budget'] / proposed_total_budget 
        proposed_dict[prog_name]['budget'] = proposed_dict[prog_name]['budget'] * orig_total_budget / proposed_total_budget
        constrained_total_budget += proposed_dict[prog_name]['budget']
        
    return proposed_dict, orig_total_budget, proposed_total_budget, constrained_total_budget

def reconciliationMetric(new_attributes, proj, parset, progset, parset_name, impact_pars, results, attribute_dict, reconcile_for_year, compareoutcome, prog_budget_alloc, constrain_budget):
    '''Objective function for reconciliation process, is used to compare outcomes as well as they use the same logic
       Uses functionality from model.py
        
        Params:
            new_attributes          Project object to run simulations for reconciliation process (type: Python object)
            proj                    Python object, i.e. the country project under consideration           
            parset                  Parameter set to match/reconcile  (type: python object)
            progset                 Program set to match/reconcile  (type: python object)
            impact_pars             Impact pars to test (type: odict)
            results                 Result set to access population compartment sizes (from model runs)
            attribute_dict          Attribute dictionary for comparison purposes i.e. to convert optimization proposed list into a readable dictionary
            reconcile_for_year      Reconciliation year
            compareoutcome          Flag to identify which parameter to return
            
        Returns:
            impact['net difference']    Return difference in progset and parset for reconciliation to minimise
            impact                      Dictionary of original and reconciled parset/progset impact comparison
    '''
    #Options
    if compareoutcome == False:
        proposed_dict = regenerateAttributesDict(new_attributes, attribute_dict)
        if constrain_budget: 
            constrained_dict, orig_total_budget, proposed_total_budget, constrained_total_budget = rescaleAllocation(proposed_dict, attribute_dict)
            new_dict = constrained_dict
        else: new_dict = proposed_dict
        progset = updateProgset(new_dict, progset, year=reconcile_for_year)
        if constrain_budget:
            print('Total Budget: %g\tProposed Budget: %g\tConstrained Budget: %g\n' %(orig_total_budget, proposed_total_budget, constrained_total_budget))
            for key in attribute_dict.keys():
                print('Program: %s\n  Original Unit Cost: %g\tProposed Unit Cost: %g\tConstrained Unit Cost: %g' %(key, attribute_dict[key]['unit_cost'], proposed_dict[key]['unit_cost'], constrained_dict[key]['unit_cost']))
                print('  Original Budget: %g\tProposed Budget: %g\tConstrained Budget: %g\n' %(attribute_dict[key]['budget'], proposed_dict[key]['budget'], constrained_dict[key]['budget']))
            
    par_attributes = odict()
    prog_attributes = odict()
    parser = FunctionParser(debug=False)
    if prog_budget_alloc is None: prog_budget_alloc = {}

    if proj.settings.tvec_start < reconcile_for_year < proj.settings.tvec_end:
        raise OptimaException('Year for reconciliation (%.0f) is not in the simulation interval from %0.f to %0.f!' % (reconcile_for_year, proj.settings.tvec_start, proj.settings.tvec_end))

    rel_t_idx = reconcile_for_year - proj.settings.tvec_start

    
    for par_label in impact_pars:
        for popkey in results.pop_label_index.keys():
            if popkey not in prog_attributes.keys(): prog_attributes[popkey] = odict()
            if par_label in progset.impacts.keys():
                first_prog = True   # True if program in prog_label loop is the first one in the impact dict list.
                impact_list = []    # Notes for each program what its impact would be before coverage limitations.
                overflow_list = []  # Notes for each program how much greater its funded coverage is than people available to be covered.
                if par_label not in prog_attributes[popkey].keys(): prog_attributes[popkey][par_label] = odict()

                pop = results.m_pops[results.pop_label_index[popkey]]
                if par_label in proj.settings.par_deps:
                    par = pop.getDep(par_label)
                else:
                    par = pop.getLinks(proj.settings.linkpar_specs[par_label]['tag'])[0]
                new_val = par.vals[rel_t_idx]

                # get all the relevant programs for the current parameter in the current population
                rel_prog_labels = filter(lambda x: popkey in progset.getProg(x).target_pops, progset.impacts[par_label])
                # add programs from programs which have a global impact but are not covered in this population
                rel_prog_labels = buildRelevantProgs(par_label, popkey, rel_prog_labels, progset)

                for prog_label in rel_prog_labels:
                    prog = progset.getProg(prog_label)
                    prog_type = prog.prog_type

                    # if no budget allocation is given, fill in the allocation with default values
                    if prog_label not in prog_budget_alloc.keys():
                        prog_budget = prog.getDefaultBudget()
                        prog_budget_alloc[prog_label] = prog_budget
                    else:
                        prog_budget = prog_budget_alloc[prog_label]

                    (special, scale_pars) = progset.getProg(prog_label).flag
                    net_impact = prog.getImpact(prog_budget, impact_label=par_label, parser=parser, years=[reconcile_for_year])

                    if 'tag' in proj.settings.linkpar_specs[par_label]:
                        tag = proj.settings.linkpar_specs[par_label]['tag']
                        source_element_size = pop.comps[par.index_from[1]].popsize[-1]
                        source_set_size = 0
                        for from_pop_label in prog.target_pops:
                            from_pop = results.m_pops[results.pop_label_index[from_pop_label]]
                            alt_par = from_pop.getLinks(tag)[0]
                            source_set_size += from_pop.comps[alt_par.index_from[1]].popsize[-1]
                        # Coverage is also split across the source compartments of grouped impact parameters, as specified in the cascade sheet.
                        if 'group' in proj.settings.progtype_specs[prog_type]['impact_pars'][par_label]:
                            group_label = proj.settings.progtype_specs[prog_type]['impact_pars'][par_label]['group']
                            for alt_par_label in proj.settings.progtype_specs[prog_type]['impact_par_groups'][group_label]:
                                if not alt_par_label == par_label:
                                    alt_pars = pop.getLinks(proj.settings.linkpar_specs[alt_par_label]['tag'])
                                    for from_pop_label in prog.target_pops:
                                        from_pop = results.m_pops[results.pop_label_index[from_pop_label]]
                                        source_set_size += from_pop.comps[alt_pars[0].index_from[1]].popsize[-1]


                        if par.val_format == 'fraction':
                            if prog.cov_format == 'fraction':
                                impact = float(net_impact)
                            elif prog.cov_format == 'number':
                                if source_element_size <= settings.TOLERANCE:
                                    impact = 0.0
                                else:
                                    impact = net_impact / source_set_size
                        elif par.val_format == 'number':
                            if prog.cov_format == 'fraction':
                                impact = net_impact * source_element_size
                            elif prog.cov_format == 'number':
                                if source_element_size <= settings.TOLERANCE:
                                    impact = 0.0
                                else:
                                    impact = net_impact * source_element_size / source_set_size

                        if special != 'supp':
                            # Calculate how excessive program coverage is for the provided budgets.
                            # If coverage is a fraction, excess is compared to unity.
                            # If coverage is a number, excess is compared to the total number of people available for coverage.
                            overflow_factor = 0.
                            net_cov = prog.getCoverage(prog_budget)
                            if prog.cov_format == 'fraction':
                                overflow_factor = net_cov
                            elif prog.cov_format == 'number':
                                if float(source_set_size) <= settings.TOLERANCE:
                                    overflow_factor = np.inf
                                else:
                                    overflow_factor = net_cov / float(source_set_size)
                            overflow_list.append(overflow_factor)
                    
                    else:
                        # if coverage format is a number, it is difficult to decide what to do. not doing anything seems alright
                        if special != 'supp' and prog.cov_format == 'fraction':
                            overflow_list.append(prog.getCoverage(prog_budget))
                        impact = net_impact

                    if special == 'scale_prop' and par_label in scale_pars:
                        impact = processScalePropsTagReconcile(par_label, results.m_pops, results.pop_label_index, progset, prog_budget_alloc, parser, rel_t_idx)
                    # elif special == 'split':
                    #     if 'group' in settings.progtype_specs[prog_type]['impact_pars'][par_label]:
                    #         assert (0. <= source_element_size / source_set_size <= 1.)
                    #         impact *= source_element_size / source_set_size
                    # TODO might be reasonable to move this elif part up because impact is computed only here
                    elif special == 'supp':
                        # only apply impact if the paramter is not influenced by GPs; skip otherwise:
                        # if the intersection of impacting programs, GPs and available programs is empty or no global impact parameter is referenced, apply changes
                        if not set(progset.impacts[par_label]).intersection(progset.GPs) and not par_label in progset.GP_pars:
                            impact = processSuppTagReconcile(par_label, prog, progset, prog_budget_alloc, reconcile_for_year)
                            # necessary work-around to prevent supp-programs to change values when they do not have any influence
                            if impact == 0.: continue
                        else:
                            continue

                    if special != '' and special != 'replace':
                        if special == 'scale' and par_label in scale_pars:
                            # !currently only fractions supported!
                            link = pop.getLinks(prog.deps[par_label])[0]
                            impact *= pop.comps[link.index_from[1]].popsize[rel_t_idx]

                        else:
                            if par.val_format == 'number' and isinstance(par, Link):
                                impact += par.vals[rel_t_idx]
                            else:
                                impact *= par.vals[rel_t_idx]
                    else:
                        # originally: new_val = 0 and append impact to impact_list; however, the impact_list may be scaled
                        # depending on the overflow. This should not be the case for number of treatments (the only case
                        # applicable here). so instead of appending the value, overwrite new_val and do not append to impact_list
                        new_val = impact
                        continue

                    # new_val += impact
                    impact_list.append(impact)
                    # prog_attributes[popkey][par_label]['Coverage Cap Impact Value'] = np.nan

                if len(impact_list) > 0:
                    prog_attributes[popkey][par_label]['Original Impact Value'] = new_val + np.sum(impact_list)

                # Checks to make sure that the net coverage of all programs targeting a parameters is capped by those that are available to be covered.
                # Otherwise renormalises impacts.
                if len(overflow_list) > 0 and sum(overflow_list) > 1:
                    impact_list = np.multiply(impact_list, 1/sum(overflow_list))

                if len(impact_list) > 0:
                    new_val += np.sum(impact_list)

                if prog_attributes[popkey][par_label]:    
                    prog_attributes[popkey][par_label]['Coverage Cap Impact Value'] = new_val
                    prog_attributes[popkey][par_label]['overflow_list'] = [format(x, '.2f') for x in overflow_list]
                    
                
                    
    ###############################################################################
    ##Cleanup prog_attributes dictionary if empty odicts exist
    # TODO this step may be skipped?
    for popkey in prog_attributes.keys():
        for par_label in prog_attributes[popkey]:
            if (type(prog_attributes[popkey][par_label]) == odict) and (bool(prog_attributes[popkey][par_label]) == False):
                del prog_attributes[popkey][par_label]
    ###############################################################################
    ##Calculate Impact of parset

    # create parameter key to index in parset dict for convenient access
    key2idx = dict(zip(map(lambda x: x.label, parset), range(len(parset))))
    for par_label in impact_pars:
        for popkey in results.pop_label_index.keys():
            pop = results.m_pops[results.pop_label_index[popkey]]
            if not popkey in par_attributes: par_attributes[popkey] = odict()
            par_attributes[popkey][par_label] = odict()
            specs = proj.settings.linkpar_specs[par_label]

            # Calculate the value of a parameter from its function if it exists, otherwise maintain the current value.
            if 'f_stack' in specs.keys():
                f_stack = dcp(specs['f_stack'])
                deps = dcp(specs['deps'])
                for dep_label in deps.keys():
                    if dep_label in proj.settings.par_deps.keys() or dep_label in proj.settings.charac_deps.keys():
                        val = pop.getDep(dep_label).vals[rel_t_idx]
                    else:
                        val = pop.getLinks(proj.settings.linkpar_specs[dep_label]['tag'])[0].vals[rel_t_idx]  # As links are duplicated for the same tag, can pull values from the zeroth one.
                    deps[dep_label] = val
                new_val = parser.evaluateStack(stack=f_stack, deps=deps)
            else:
                new_val = parset[key2idx[par_label]].interpolate(tvec=np.arange(proj.settings.tvec_start, proj.settings.tvec_end + proj.settings.tvec_dt/2, proj.settings.tvec_dt), pop_label=popkey)[-1]

            par_attributes[popkey][par_label]['Impact Value'] = new_val

    ###############################################################################
    ##Cleanup par_attributes dictionary to match prog_attributes
    # TODO this step may be skipped?
    for popkey in prog_attributes.keys():
        for par_label in prog_attributes[popkey]:
            if par_label in par_attributes[popkey].keys():
                continue
            else: del par_attributes[popkey[par_label]]
    ###############################################################################
    ##Create a single comparison dictionary
    impact = odict()
    
    #return sum of the difference squared to avoid negative numbers
    impact['net_difference'] = 0.0
    for popkey in prog_attributes.keys():
        for par_label in prog_attributes[popkey].keys():
            if par_label not in impact.keys(): impact[par_label] = odict()
            if popkey not in impact[par_label].keys(): impact[par_label][popkey] = odict()
            temp_parset_impact = par_attributes[popkey][par_label]['Impact Value']
            temp_progset_impact = prog_attributes[popkey][par_label]['Coverage Cap Impact Value']
            difference = (temp_parset_impact - temp_progset_impact)**2
            impact[par_label][popkey] = {'Impact Difference': difference}
            impact['net_difference'] += difference    
    if compareoutcome == False: 
        return impact['net_difference']
    else:
        #return comparison between progset and parset
        for popkey in prog_attributes.keys():
            for par_label in prog_attributes[popkey].keys():
                #if popkey not in impact.keys(): impact[popkey] = odict()
                #if par_label not in impact[popkey].keys(): impact[popkey][par_label] = odict()
                #impact[popkey][par_label]['parset_impact_value'] = par_attributes[popkey][par_label]['Impact Value']
                #impact[popkey][par_label]['progset_impact_value'] = prog_attributes[popkey][par_label]['Impact Value']
                if par_label not in impact.keys(): impact[par_label] = odict()
                if popkey not in impact[par_label].keys(): impact[par_label][popkey] = odict()
                impact[par_label][popkey]['parset_impact_value'] = par_attributes[popkey][par_label]['Impact Value']
                impact[par_label][popkey]['progset_impact_uncapped'] = prog_attributes[popkey][par_label]['Original Impact Value']
                impact[par_label][popkey]['progset_impact_capped'] = prog_attributes[popkey][par_label]['Coverage Cap Impact Value']
                impact[par_label][popkey]['overflow_list'] = prog_attributes[popkey][par_label]['overflow_list']

        return impact
