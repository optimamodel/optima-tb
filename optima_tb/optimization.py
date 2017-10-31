from optima_tb.asd import asd
from optima_tb.model import runModel
from optima_tb.utils import odict, OptimaException, tic, toc
from optima_tb.defaults import defaultOptimOptions
import optima_tb.settings as project_settings

import logging
import logging.config
logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

import numpy as np
from copy import deepcopy as dcp

def constrainAllocation(alloc, settings, options, algorithm_refs, attempt=0):

#    logging.info('Launching attempt %i to constrain allocation.' % attempt)

    alloc = dcp(np.array(alloc, dtype=np.float64))    # Converting to np array just in case.

    # Convert negative allocation values to zeros.
    alloc[alloc < 0.0] = 0.0

#    print('Attempt: %i' % attempt)

    # Handle values that go beyond limits.
    hit_lower = (alloc == np.nan)  # Track the indices of programs that hit the lower limit. Pre-allocate as False.
    hit_upper = (alloc == np.nan)  # Track the indices of programs that hit the upper limit. Pre-allocate as False.
    for k in xrange(len(alloc)):
#        print 'Old'
#        print alloc[k]
        prog_key = algorithm_refs['alloc_ids']['id_progs'][k]
#        print options['constraints']['limits'][prog_key]['vals']
#        print options['orig_alloc'][prog_key]
        if options['constraints']['limits'][prog_key]['rel']:

            check = options['constraints']['limits'][prog_key]['vals'][1] * options['orig_alloc'][prog_key]
            if check is np.nan: check = options['constraints']['limits'][prog_key]['vals'][1]   # Avoids infinity times zero case.
            # If new budget is less than the relative lower limit times original budget...
            if alloc[k] <= options['constraints']['limits'][prog_key]['vals'][0] * options['orig_alloc'][prog_key]:
                alloc[k] = options['constraints']['limits'][prog_key]['vals'][0] * options['orig_alloc'][prog_key]
                hit_lower[k] = True
#                print options['constraints']['limits'][prog_key]['vals'][0] * options['orig_alloc'][prog_key]
            # If new budget is more than the relative upper limit times original budget...
            elif alloc[k] >= check:
                alloc[k] = check
                hit_upper[k] = True
#                print check
        else:
            # If new budget is less than the absolute lower limit...
            if alloc[k] <= options['constraints']['limits'][prog_key]['vals'][0]:
                alloc[k] = options['constraints']['limits'][prog_key]['vals'][0]
                hit_lower[k] = True
            # If new budget is more than the absolute upper limit...
            elif alloc[k] >= options['constraints']['limits'][prog_key]['vals'][1]:
                alloc[k] = options['constraints']['limits'][prog_key]['vals'][1]
                hit_upper[k] = True
#        print 'New'
#        print alloc[k]

    if 'total' in options['constraints']:
#        print hit_upper
#        print hit_lower
        sum_current = sum(alloc)

        # Ensure that there is never a situation where no budget at all is rescalable.
        if sum_current == 0:
            logger.warn('Allocation had a sum of zero prior to rescaling during the optimization process. All budgets incremented by 1.')
            alloc += 1
            sum_current = sum(alloc)

#        print alloc
        cannot_change = dcp(hit_upper)  # Just to define it as something.
        if options['constraints']['total'] > sum_current:       # Need to scale up.
            cannot_change = dcp(hit_upper)                          # Cannot change values that have already hit upper limit.
        elif options['constraints']['total'] < sum_current:     # Need to scale down.
            cannot_change = dcp(hit_lower)                          # Cannot change values that have already hit lower limit.

        # Ensure that there is never a situation where no modifiable budget is rescalable.
        if sum(alloc[~cannot_change]) == 0:
            logger.warn('Modifiable components of allocation had a sum of zero prior to rescaling during the optimization process. All modifiable budgets incremented by 1.')
            alloc[~cannot_change] += 1
            sum_current = sum(alloc)

        sum_stuck = sum(alloc[cannot_change])   # The total budget amount that is stuck at its limit.
#        print cannot_change
#        print alloc[cannot_change]
        if sum_current - sum_stuck == 0.0:
            logger.warn("An optimization iteration has pushed an allocation entirely to its constraint limits. Not rescaling.")
        else:
            alloc[~cannot_change] *= (options['constraints']['total'] - sum_stuck) / (sum_current - sum_stuck)
#        print (options['constraints']['total']-sum_stuck)/(sum_current-sum_stuck)
#        print alloc
        # Recursively constrain until the budget total rescale does nothing, or the recursive limit gets hit.
        if not abs(sum(alloc) - sum_current) < project_settings.TOLERANCE:
#            logger.info('Allocation not constrained on attempt %i: %f > %s' % (attempt, abs(sum(alloc) - sum_current), project_settings.TOLERANCE))
            if attempt < settings.recursion_limit:
                alloc = constrainAllocation(alloc=alloc, settings=settings, options=options, algorithm_refs=algorithm_refs, attempt=attempt + 1)
            else:
                logger.warn("Tried to constrain an allocation but failed before recursion limit. Reverting to previous iteration allocation.")
                alloc = dcp(algorithm_refs['previous_alloc'])
#        else:
#            logger.info('Budget successfully constrained after %i attempts' % attempt)


    return dcp(alloc)

def calculateObjective(alloc, settings, parset, progset, options, algorithm_refs):
    '''
    Calculates the objective function value for a certain allocation.
    The algorithm_refs dictionary is for storing algorithm-useful values during the process, e.g. the previous allocation tested.
    '''

    logging.debug('Unconstrained alloc...')
    logging.debug(alloc)
    logging.debug(sum(alloc))

    alloc = constrainAllocation(alloc=alloc, settings=settings, options=options, algorithm_refs=algorithm_refs)

    logging.debug('Constrained alloc...')
    logging.debug(alloc)
    logging.debug(sum(alloc))

    if algorithm_refs is None: algorithm_refs = {}

    # Makes sure that time is not wasted on an allocation that was run in the previous step.
    # Avoids useless iteration tests where a program is given negative funding, then constrained to be zero once more.
    # This paragraph is for checking identical allocations.
    identical_alloc = False
    if 'previous_alloc' in algorithm_refs:
        if not algorithm_refs['previous_alloc'] is None:
            identical_alloc = True
            for k in xrange(len(alloc)):
                if not alloc[k] == algorithm_refs['previous_alloc'][k]:
                    identical_alloc = False
                    break
        algorithm_refs['previous_alloc'] = dcp(alloc)   # Store the new allocation for the next iteration.

    # This if-else clause prevents the actual expensive model run.
    if identical_alloc and 'previous_results' in algorithm_refs and not algorithm_refs['previous_results'] is None:
        results = algorithm_refs['previous_results']
        logger.warn("Optimization iteration is testing the same allocation. Skipping model run.")
    else:
        options_iter = dcp(options)

        if 'alloc_ids' in algorithm_refs and 'id_progs' in algorithm_refs['alloc_ids']:
            for k in xrange(len(alloc)):
                options_iter['init_alloc'][algorithm_refs['alloc_ids']['id_progs'][k]] = alloc[k]
        else:
#            logger.warn("No allocation-id program-label conversion dictionary can be found during optimization. \nMaking a potentially dangerous assumption that optimization alloc is in the same order as init alloc in the options dictionary.")
#            for k in xrange(len(alloc)):
#                options_iter['init_alloc'][k] = alloc[k]
            raise OptimaException('ERROR: No allocation-id program-label conversion dictionary can be found during optimization. Cannot push alloc back into an options dictionary.')

        t = tic()
        results = runModel(settings=settings, parset=parset, progset=progset, options=options_iter)
        logging.info(toc(t))

    if 'previous_results' in algorithm_refs:
        algorithm_refs['previous_results'] = dcp(results)   # Store the new results for the next iteration.

#    index_start = np.where(results.sim_settings['tvec'] >= options['progs_start'])[0][0]
#    logging.debug( index_start)
#    logging.debug( results.sim_settings['tvec'][index_start:])
    objective = 0.0

    for objective_label in options['objectives']:
#        charac_label = options['objective_weight'].keys()[0]
        weight = 1.0
        if 'weight' in options['objectives'][objective_label]: weight = options['objectives'][objective_label]['weight']


#        objective = results.getValuesAt(label = objective_label, year_init = options['progs_start'])

        # Handle an objective that wants to know the value in a particular year, e.g. a snapshot of a characteristic.
        if 'year' in options['objectives'][objective_label]:
#            index_start = np.where(results.sim_settings['tvec'] >= options['objectives'][objective_label]['year'])[0][0]
#            for pop_label in results.outputs[objective_label].keys():
#                objective += results.outputs[objective_label][pop_label][index_start] * weight
#            index_start = np.where(results.sim_settings['tvec'] >= options['progs_start'])[0][0]
#             logging.debug("getValuesAt")
#             logging.debug(results.getValuesAt(label=objective_label, year_init=options['objectives'][objective_label]['year'])[0][0] * weight)
#             logging.debug("getValueAt")
#             logging.debug(results.getValueAt(label=objective_label, year_init=options['objectives'][objective_label]['year']) * weight)
#             logging.info("Obj label = %s, eval = %.2f" % (objective_label, results.getValuesAt(label=objective_label, year_init=options['objectives'][objective_label]['year'])[0][0] * weight))
            objective += results.getValuesAt(label=objective_label, year_init=options['objectives'][objective_label]['year'])[0][0] * weight
        # Handle the default objective that wants to know a cumulative value. Timestep scaling is done inside Results.getValueAt().
        else:
#            for pop_label in results.outputs[objective_label].keys():
#                objective += sum(results.outputs[objective_label][pop_label][index_start:]) * results.dt * weight
            objective += results.getValuesAt(label=objective_label, year_init=options['progs_start'], year_end=settings.tvec_end, integrated=True)[0][0] * weight

#     logging.info(objective)

    return objective


def optimizeFunc(settings, parset, progset, options=None, outputqueue=None, thread=None, randseed=None, **optimization_params):

    if options is None:
        logger.warn("An options dictionary was not supplied for optimisation. A default one will be constructed.")
        options = defaultOptimOptions(settings=settings, progset=progset)

    if 'progs_start' not in options:
        options['progs_start'] = 2015.0
        logger.warn("Programs will start overwriting calibrated parameter values in the year: %f" % options['progs_start'])

    if 'init_alloc' not in options:
        options['init_alloc'] = {}
    if 'constraints' not in options:
        options['constraints'] = {'limits':odict()}

    for prog in progset.progs:
        if prog.label not in options['init_alloc']:
            # If this options is off, only a limited subset will be optimized.
            if 'saturate_with_default_budgets' in options and options['saturate_with_default_budgets'] is True:
                options['init_alloc'][prog.label] = prog.getDefaultBudget()

        # If saturation is on, make sure the original allocation is also saturated appropriately.
        # Assumes saturation tag does not change between optimisations.
        if 'orig_alloc' in options:
            if prog.label not in options['orig_alloc']:
                if 'saturate_with_default_budgets' in options and options['saturate_with_default_budgets'] is True:
                    options['orig_alloc'][prog.label] = prog.getDefaultBudget()

        # For programs chosen to be optimised, make sure proper constraints exist.
        if prog.label in options['init_alloc']:
            if prog.label not in options['constraints']['limits']:
                    options['constraints']['limits'][prog.label] = {'vals':[0.0, np.inf], 'rel':True}
                    if prog.func_specs['type'] == 'cost_only':
                        options['constraints']['limits'][prog.label]['vals'] = [1.0, 1.0]
            else:
                if 'vals' not in options['constraints']['limits'][prog.label]:
                    raise OptimaException('ERROR: Limit constraints specified for program "%s", but with no vals defined.' % prog.label)
                elif len(options['constraints']['limits'][prog.label]['vals']) != 2:
                    raise OptimaException('ERROR: Limit constraints for program "%s" must contain a two-element list keyed by "vals", specifying min and max limits.' % prog.label)

#    # Convert alloc into an ordered dictionary so that keys and values are ordered during optimisation.
#    options['init_alloc'] = odict(options['init_alloc'])

    # If user has not supplied constraints, then make the optimisation a redistribution of default budgets.
    if 'total' not in options['constraints']:
        options['constraints']['total'] = sum(options['init_alloc'].values())

    if 'objectives' not in options:
        options['objectives'] = {settings.charac_pop_count : {'weight':-1, 'year':2030.0}}

#     logging.debug( options)

#    logging.debug( options)
#    logging.debug( total_budget)

#    results = runModel(settings = settings, parset = parset, progset = progset, options = options)
#
#    index_start = np.where(results.sim_settings['tvec'] >= options['progs_start'])[0][0]
#    logging.debug( index_start)
#    logging.debug( results.sim_settings['tvec'][index_start:])
#    charac_label = options['objective_weight'].keys()[0]
#    objective = 0
#    for pop_label in results.outputs[charac_label].keys():
#        objective += sum(results.outputs[charac_label][pop_label][index_start:])*results.dt*options['outcome_weight'][charac_label]
#    logging.debug( objective)

    algorithm_refs = {'previous_alloc':None, 'previous_results':None, 'alloc_ids':{'id_progs':{}, 'prog_ids':{}}}

    # Flattens out the initial allocation into a list, but makes sure to note which program labels link to which allocation indices.
    alloc = []
#    algorithm_refs['alloc_ids'] = {}
    k = 0
    for prog_key in options['init_alloc'].keys():
        alloc.append(options['init_alloc'][prog_key])
        algorithm_refs['alloc_ids']['id_progs'][k] = prog_key
        algorithm_refs['alloc_ids']['prog_ids'][prog_key] = k
        k += 1
    alloc = dcp(np.array(alloc))

    # Create an original allocation backup of initial allocation if it does not exist.
    # For standard runs 'orig_alloc' should be identical to 'init_alloc', particularly as optimisation iterations work on duplicated options dictionaries, not the original.
    # For block-optimisations, where options dicts are passed back in with an updated 'init_alloc', 'orig_alloc' should remain fixed and reference the allocation at the very start of the process.
    if 'orig_alloc' not in options:
        options['orig_alloc'] = dcp(options['init_alloc'])

#     logging.debug( alloc )
#     logging.debug( algorithm_refs )

    args = {'settings':settings,
            'parset':parset,
            'progset':progset,
            'options':options,
            'algorithm_refs':algorithm_refs}

    algorithm_refs['previous_alloc'] = dcp(alloc)

    alloc = dcp(constrainAllocation(alloc=alloc, settings=settings, options=options, algorithm_refs=algorithm_refs))
    alloc_new, obj_vals, exit_reason = asd(calculateObjective, alloc, args=args, randseed=randseed, **optimization_params)
    alloc_new = dcp(constrainAllocation(alloc=alloc_new, settings=settings, options=options, algorithm_refs=algorithm_refs))

    # Makes sure allocation is returned in dictionary format.
    alloc_opt = {}
    if 'alloc_ids' in algorithm_refs and 'id_progs' in algorithm_refs['alloc_ids']:
        for k in xrange(len(alloc_new)):
            alloc_opt[algorithm_refs['alloc_ids']['id_progs'][k]] = alloc_new[k]
    else:
        raise OptimaException('ERROR: No allocation-id program-label conversion dictionary can be found during optimization. Cannot convert alloc back into options dictionary format.')

    results = (alloc_opt, obj_vals, exit_reason)

    # For batch, we add this to the output queue
    if outputqueue is not None:
        outputqueue.put([thread, results])

    return results


def parallelOptimizeFunc(settings, parset, progset, options=None, num_threads=4,
                         block_iter=10, max_blocks=10, doplot=False,
                         fullfval=False, randseed=None, **parallel_optimization_params):
    '''
    Same as optimizeFunc, excepts runs in multiple threads in small blocks.
    
    Inputs:
        settings, parset, progset, options = standard optimizeFunc() inputs
        num_threads = the number of parallel threads to run
        block_iter = the number of iterations per block
        max_iter = optional alternate way of specifying the number of iterations
        max_blocks = the total number of blocks to run
        doplot = whether or not to plot the function values (fvals) at every iteration
        fullfval = whether or not to return the full array of fvals
        
    
    In total,this function will run num_threads*block_iter*max_blocks iterations,
    but in block_iter*max_blocks time. What happens is it runs a block of <block_iter>
    iterations on each core, and at the end of the block, picks the best one, and uses
    that as the starting point for the next block. For example, with default values of
    num_threads = 4, block_iter = 10, max_blocks = 10, you can run 400 iterations in
    the time it would take 100 to run normally.
    
    NOTE 1: this switching between blocks eliminates the learning that happens during ASD,
    which will slow things down. However, it would be a fair bit of work to have these
    outputted and then passed back in as well. For another day :)
    
    NOTE 2: If you set max_blocks = 1 and block_iter = max_iter, then this will just run 
    a bunch of ASD instances in parallel, which might also be useful!
        
    
    Usage: you can run it exactly the same way as proj.optimize():
        
        
    
    '''
    # Import dependencies here so no biggiei if they fail
    from multiprocessing import Process, Queue
    from time import time

    # Optionally plot
    if doplot:
        from pylab import figure, plot, clf
        figure()

    # Handle inputs
    if options is None:
        options = defaultOptimOptions(settings=settings, progset=progset)

    if block_iter is None and max_iters is not None and max_blocks is not None:
        logging.info("Parallel optimization: block_iter was not specified; calculating it from max_ters and max_blocks")
        block_iter = np.floor(max_iters / float(max_blocks)) # Calculate block_iter from max_iter, if supplied

    # Crucial step that ensures an original allocation is always stored, regardless of initial allocation changes per block.
    if 'orig_alloc' not in options:
        options['orig_alloc'] = dcp(options['init_alloc'])

    total_iters = block_iter * max_blocks
    fvalarray = np.zeros((num_threads, total_iters)) + np.nan

    msg = "Starting a parallel optimization with %i threads for %i iterations each for %i blocks" % (num_threads, block_iter, max_blocks)
    logger.info(msg)

    # Loop over the optimization blocks
    for block in range(max_blocks):

        # Set up the parallel process
        outputqueue = Queue()
        outputlist = np.empty(num_threads, dtype=object)
        processes = []

        # Loop over the threads, starting the processes
        for thread in range(num_threads):
            if randseed is None:
                randseed = (block + 1) * int((time() - np.floor(time())) * 1e7) # Get a random number based on both the time and the thread
            args = (settings, parset, progset, options, outputqueue, thread, randseed)
            prc = Process(target=optimizeFunc, args=args, kwargs=parallel_optimization_params)
            prc.start()
            processes.append(prc)

        # Tidy up: close the threads and gather the results
        for i in range(num_threads):
            thread, result = outputqueue.get()  # This is needed or else the process never finishes
            outputlist[thread] = result
        for prc in processes:
            prc.join() # Wait for them to finish

        # Figure out which one did best
        bestfvalval = np.inf
        bestfvalind = None

        for i in range(num_threads):
            fvalarray[i, block * block_iter:(block + 1) * block_iter] = outputlist[i][1][-block_iter:] # in case the first timestep is included
            thisbestval = outputlist[i][1][-1]
            if thisbestval < bestfvalval:
                bestfvalval = thisbestval
                bestfvalind = i

        options['init_alloc'] = outputlist[bestfvalind][0] # Update the budget and use it as the input for the next block -- this is key!

        if doplot:
            clf()
            plot(np.transpose(fvalarray))

    results = (options['init_alloc'], fvalarray[bestfvalind, :], 'Parallel optimization')

    if fullfval: results[1] = fvalarray # Use the full fvalarray instead of just the best one

    return results
