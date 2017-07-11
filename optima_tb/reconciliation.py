from optima_tb.utils import odict, OptimaException
import optima_tb.settings as settings
from copy import deepcopy as dcp
from optima_tb.asd import asd
import numpy as np
from optima_tb.parsing import FunctionParser

import logging
logger = logging.getLogger(__name__)




def reconcileFunc(proj, reconcile_for_year, parset_name, progset_name, unitcost_sigma = 0.05, attribute_sigma = 0.20, impact_pars = None, budget_allocation = None, orig_tvec_end = None):
        """
        Reconciles progset to identified parset, the objective being to match the parameters as closely as possible with identified standard deviation sigma
        
        Params:
            proj                    Project object to run simulations for reconciliation process (type: Python object)
            reconcile_for_year      Year for which reconciliation needs to be done to ensure continuity of parset/progset (type: int)
            parset_name             Parameter set name to match/reconcile  (type: string)
            progset_name            Program set name to match/reconcile  (type: string)
            unitcost_sigma          Standard deviation allowable for Unit Cost (type: float)
            attribute_sigma         Standard deviation allowable for attributes identified in impact_pars (type: float)
            impact_pars             Impact pars to be reconciled (type: list or None)
            budget_allocation       Dictionary of programs with new budget allocations (type: dict)
            
        Returns:
            progset                 Updated progset with reconciled values
            impact                  dictionary of original and reconciled parset/progset impact comparison
            
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
        
        results = proj.runSim(parset_name=parset_name)
        
        #Get original comparison between progset and parset
        impact['original'] = compareOutcomesFunc(proj=proj, parset_name=parset_name, progset_name=progset_name, year=reconcile_for_year, compareoutcome=True, display=False)
        
        #Convert into an optimisable list
        attribute_dict = createAttributeDict(settings=proj.settings, progset=progset)
        attribute_list, unitcost_index = createAttributeList(attribute_dict=attribute_dict)
        #Setup min-max bounds for optimisation
        xmin, xmax = dcp(attribute_list), dcp(attribute_list)
        
        #Setup xmin and xmax conditions, unit_cost indexes are used in case unit_costs have a difference standard deviation value
        for index in range(len(attribute_list)):
            if index in unitcost_index:
                xmin[index] *= (1-unitcost_sigma)
                xmax[index] *= (1+unitcost_sigma)
            else:
                xmin[index] *= (1-attribute_sigma)
                xmax[index] *= (1+attribute_sigma)
                if xmax[index] > 1.: xmax[index] = 1.
            
            if xmin[index] <= 0.: xmin[index] = settings.TOLERANCE
        
        #Run optimisation
        args = {'proj': proj, 'parset': parset, 'progset': progset, 'parset_name': parset_name,
                'impact_pars': impact_pars, 'results': results, 'attribute_dict': attribute_dict, 
                'reconcile_for_year': reconcile_for_year, 'compareoutcome': False, 'prog_budget_alloc': budget_allocation}
        
        optim_args = {
                     'stepsize': proj.settings.autofit_params['stepsize'], 
                     'maxiters': 300,#proj.settings.autofit_params['MaxIter'],
                     'maxtime': proj.settings.autofit_params['timelimit'],
                     'sinc': proj.settings.autofit_params['sinc'],
                     'sdec': proj.settings.autofit_params['sdec'], 
                     'fulloutput': False,
                     'reltol': None
                     }
      
        
        best_attribute_list = asd(reconciliationMetric, attribute_list, args, xmin = xmin, xmax = xmax, **optim_args)
        best_attribute_dict = regenerateAttributesDict(attribute_list = best_attribute_list, orig_attribute_dict = attribute_dict)
        progset = updateProgset(new_pars_dict=best_attribute_dict, progset=progset)
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
                    outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\n'.format(popkey, impact['original'][par_label][popkey]['parset_impact_value'], impact['original'][par_label][popkey]['progset_impact_value'], impact['reconciled'][par_label][popkey]['progset_impact_value'])
                outcome += '\n'
        print outcome
        #Reset back to original runSim durations
        proj.setYear([2000, orig_tvec_end], False)
        return progset

def compareOutcomesFunc(proj, year, parset_name=None, progset_name=None, budget_allocation=None, compareoutcome=None, display=True):
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
    results = proj.runSim(parset_name=parset_name)
    #Declare impact parameters to use for display (use all)
    impact_pars = progset.impacts.keys()
    #use reconcilemetric to get desired result
    impact = reconciliationMetric(new_attributes=[], proj=proj, parset=parset, progset=progset, 
                                  parset_name=parset_name, impact_pars=impact_pars, 
                                  results=results, attribute_dict={}, reconcile_for_year=year, 
                                  compareoutcome=compareoutcome, prog_budget_alloc=budget_allocation)
    #display output
    if display:
        print('Comparing outcomes for year: %i' %year)
        parset_value  = 'parset_impact_value'
        progset_value_uncapped = 'progset_impact_value_uncapped'
        progset_value_capped = 'progset_impact_value_capped'
        outcome = '\n\t\t\t%s\t%s\t%s\n' %(parset_value, progset_value_uncapped, progset_value_capped)
        for par_label in impact.keys():
            if par_label == 'net_difference': continue
            else:
                outcome += '%s\n' %(par_label)
                for popkey in impact[par_label]:
                    try: outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\n'.format(popkey, impact[par_label][popkey][parset_value], impact[par_label][popkey][progset_value_uncapped], impact[par_label][popkey][progset_value_capped])
                    except:
                        try: outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\n'.format(popkey, impact[par_label][popkey][parset_value], impact[par_label][popkey][progset_value_uncapped][0], impact[par_label][popkey][progset_value_capped][0])
                        except:
                            try: outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\n'.format(popkey, impact[par_label][popkey][parset_value], impact[par_label][popkey][progset_value_uncapped][0], impact[par_label][popkey][progset_value_capped])
                            except:
                                outcome += '\t{:<10}\t{:10.2f}\t\t{:10.2f}\t\t{:10.2f}\n'.format(popkey, impact[par_label][popkey][parset_value], impact[par_label][popkey][progset_value_uncapped], impact[par_label][popkey][progset_value_capped][0])
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
        mark_for_delete = False
        if prog_label not in attributes_dict.keys(): attributes_dict[prog_label] = odict()
        index = progset.prog_ids[prog_label]
        try:    
            attributes_dict[prog_label]['unit_cost'] = progset.progs[index].func_specs['pars']['unit_cost']
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
            unitcost_index          Index of locations where unit_costs exist in the list
    '''
    attribute_list = []
    unitcost_index = []
    index = 0
    for prog_label in attribute_dict.keys():
        for par in attribute_dict[prog_label]:
            attribute_list.append(attribute_dict[prog_label][par])
            if par == 'unit_cost':
                unitcost_index.append(index)
                index += 1
            else: index += 1
    return attribute_list, unitcost_index

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

def updateProgset(new_pars_dict, progset):
    '''Update the progset with the new values as obtained from reconciliation process, only updates the last known value
    
       Params:
            new_pars_dict          Dictionary form of the optimized list which is used to update progset
            progset                Program set that is to be updated with the new values
            
       Returns:
            attribute_dict         Dictionary form of the optimized list which is used to update progset
    '''
    for prog_label in progset.prog_ids.keys():
        if prog_label not in new_pars_dict.keys(): continue
        else:
            index = progset.prog_ids[prog_label]
            for attribute in new_pars_dict[prog_label].keys():
                if attribute in progset.progs[index].attributes.keys():
                    progset.progs[index].attributes[attribute][-1] = new_pars_dict[prog_label][attribute]
                else:
                    continue
            progset.progs[index].func_specs['pars']['unit_cost'] = new_pars_dict[prog_label]['unit_cost']
    return progset

def reconciliationMetric(new_attributes, proj, parset, progset, parset_name, impact_pars, results, attribute_dict, reconcile_for_year, compareoutcome, prog_budget_alloc = None):
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
        new_dict = regenerateAttributesDict(new_attributes, attribute_dict)
        progset = updateProgset(new_dict, progset)
    par_attributes = odict()
    prog_attributes = odict()
    parser = FunctionParser(debug=False)
    
    for par_label in (impact_pars):
        for popkey in results.pop_label_index.keys():
            if popkey not in prog_attributes.keys(): prog_attributes[popkey] = odict()
            if par_label in progset.impacts.keys():
                first_prog = True   # True if program in prog_label loop is the first one in the impact dict list.
                impact_list = []    # Notes for each program what its impact would be before coverage limitations.
                overflow_list = []  # Notes for each program how much greater its funded coverage is than people available to be covered.
                if par_label not in prog_attributes[popkey].keys(): prog_attributes[popkey][par_label] = odict()
                for prog_label in progset.impacts[par_label]:
                    prog = progset.getProg(prog_label)
                    prog_type = prog.prog_type
                    # Make sure the population in the loop is a target of this program.x
                    if popkey not in prog.target_pops:
                        continue
                    if prog_budget_alloc is None or prog_label not in prog_budget_alloc.keys():
                        prog_budget = prog.getDefaultBudget()  
                    else:
                        prog_budget = prog_budget_alloc[prog_label]
                    
                    pop = results.m_pops[results.pop_label_index[popkey]]
                    tag = proj.settings.linkpar_specs[par_label]['tag']
                    par = pop.getLinks(tag)[0]
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
                    #print('Program Label: %s, Parameter: %s, Source set size: %f, source_element_size: %f' % (prog_label, par_label, source_set_size, source_element_size))
                    net_impact = prog.getImpact(prog_budget, impact_label = par_label, parser = parser, years = [reconcile_for_year])
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
                    
                    if first_prog: 
                        new_val = 0.
                        temp_val = 0.
                        first_prog = False
                    temp_val += impact
                    impact_list.append(impact)
                    prog_attributes[popkey][par_label]['Original Impact Value'] = temp_val
                    prog_attributes[popkey][par_label]['Coverage Cap Impact Value'] = np.nan
                
                # Checks to make sure that the net coverage of all programs targeting a parameters is capped by those that are available to be covered.
                # Otherwise renormalises impacts.
                if len(overflow_list) > 0 and sum(overflow_list) > 1:
                    impact_list = np.multiply(impact_list, 1/sum(overflow_list))
                new_val += np.sum(impact_list)
#                print('Pop key: %s, Par Label: %s, New Val: %s' % (popkey, par_label, new_val))
                if prog_attributes[popkey][par_label]:    
                    prog_attributes[popkey][par_label]['Coverage Cap Impact Value'] = new_val
                    
                
                    
    ###############################################################################
    ##Cleanup prog_attributes dictionary if empty odicts exist
    for popkey in prog_attributes.keys():
        for par_label in prog_attributes[popkey]:
            if (type(prog_attributes[popkey][par_label]) == odict) and (bool(prog_attributes[popkey][par_label]) == False):
                del prog_attributes[popkey][par_label]
    ###############################################################################
    ##Calculate Impact of parset
    for par_label in (impact_pars):
        for popkey in results.pop_label_index.keys():     
            if popkey not in par_attributes.keys(): par_attributes[popkey] = odict()      
            for index in range((len(parset))):
                if (par_label == parset[index].label) and (par_label not in proj.settings.par_funcs.keys()):
                    if par_label not in par_attributes[popkey].keys(): par_attributes[popkey][par_label] = odict()      
                    par_attributes[popkey][par_label] = {'Impact Value': dcp(parset[index].interpolate(tvec=np.arange(proj.settings.tvec_start, proj.settings.tvec_end + proj.settings.tvec_dt/2, proj.settings.tvec_dt), pop_label=popkey))}
                    par_attributes[popkey][par_label]['Impact Value'] = par_attributes[popkey][par_label]['Impact Value'][-1]                
                elif (par_label == parset[index].label) and (par_label in proj.settings.linkpar_specs.keys()):
                    if par_label not in par_attributes[popkey].keys(): par_attributes[popkey][par_label] = odict()                
                    if 'f_stack' in proj.settings.linkpar_specs[par_label].keys():
                        f_stack = dcp(proj.settings.linkpar_specs[par_label]['f_stack'])
                        for dependency in proj.settings.linkpar_specs[par_label]['deps']:
                            if dependency in proj.parsets[parset_name].par_ids['cascade'].keys():
                                val = proj.parsets[parset_name].par_ids['cascade'][dependency]
                                attribs = dcp(proj.parsets[parset_name].pars['cascade'][val])
                        newattrib = odict()
                        newattrib[attribs.label] = odict()
                        for key in attribs.y:
                            if key == popkey: newattrib[attribs.label][key] = attribs.y[key][-1]
                        finalattrib = odict()
                        for i, pop_label in enumerate(newattrib[attribs.label].keys()):
                                finalattrib[attribs.label] = newattrib[attribs.label][pop_label]
                        par_attributes[popkey][par_label] = {'Impact Value': parser.evaluateStack(stack = f_stack, deps = finalattrib)}
    ###############################################################################
    ##Cleanup par_attributes dictionary to match prog_attributes
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
                impact[par_label][popkey]['progset_impact_value_uncapped'] = prog_attributes[popkey][par_label]['Original Impact Value']
                impact[par_label][popkey]['progset_impact_value_capped'] = prog_attributes[popkey][par_label]['Coverage Cap Impact Value']
        return impact
