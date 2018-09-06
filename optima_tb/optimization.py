from optima_tb.double_asd import asd
from optima_tb.model import runModel
from optima_tb.utils import odict, OptimaException, tic, toc
from optima_tb.defaults import defaultOptimOptions
import optima_tb.settings as project_settings

import logging
logger = logging.getLogger(__name__)

import numpy as np
from copy import deepcopy as dcp

def _constrainAllocation(alloc, settings, options, attempt=0):

    logging.debug('Launching attempt %i to constrain allocation.' % attempt)

    alloc = dcp(np.array(alloc, dtype=np.float64))    # Converting to np array just in case.
    prog_labels = options['init_alloc'].keys()

    # Convert negative allocation values to zeros.
    alloc[alloc < 0.0] = 0.0

    # Handle values that go beyond limits.
    hit_lower = (alloc == np.nan)  # Track the indices of programs that hit the lower limit. Pre-allocate as False.
    hit_upper = (alloc == np.nan)  # Track the indices of programs that hit the upper limit. Pre-allocate as False.

    for k in xrange(len(alloc)):

        prog_key = prog_labels[k]

        if options['constraints']['limits'][prog_key]['rel']:
            raise OptimaException('This branch should not be countered because `_complete_options` should set all limits to absolute')
        else:
            # If new budget is less than the absolute lower limit...
            if alloc[k] <= options['constraints']['limits'][prog_key]['vals'][0]:
                alloc[k] = options['constraints']['limits'][prog_key]['vals'][0]
                hit_lower[k] = True
            # If new budget is more than the absolute upper limit...
            elif alloc[k] >= options['constraints']['limits'][prog_key]['vals'][1]:
                alloc[k] = options['constraints']['limits'][prog_key]['vals'][1]
                hit_upper[k] = True

    if 'total' in options['constraints'] and options['constraints']['total']:

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
                alloc = _constrainAllocation(alloc=alloc, settings=settings, options=options, attempt=attempt + 1)
            else:
                logger.warn("Tried to constrain an allocation but failed before recursion limit. Reverting to previous iteration allocation.")
                alloc = dcp(algorithm_refs['previous_alloc'])
        else:
            logger.debug('Budget successfully constrained after %i attempts' % attempt)


    return dcp(alloc)

def _calculateObjective(alloc, settings, parset, progset, options, previous_results):
    '''
    Calculates the objective function value for a certain allocation.
    The algorithm_refs dictionary is for storing algorithm-useful values during the process, e.g. the previous allocation tested.
    '''

    # First, constrain the alloc
    alloc = _constrainAllocation(alloc=alloc, settings=settings, options=options)

    # Next, if this alloc has been tried before, return the previous objective value
    if tuple(alloc) in previous_results:
        logger.info("Optimization iteration is testing the same allocation. Skipping model run.")
        return previous_results[tuple(alloc)]

    options_iter = dcp(options)
    for prog_key,val in zip(options_iter['init_alloc'].keys(),alloc):
        options_iter['init_alloc'][prog_key] = val
#        logger.info('%s - $%f' % (prog_key,val))

    results = runModel(settings=settings, parset=parset, progset=progset, options=options_iter)

    objective = 0.0
    for objective_label in options['objectives']:
        if 'weight' in options['objectives'][objective_label]:
            weight = options['objectives'][objective_label]['weight']
        else:
            weight = 1.0

        # Handle an objective that wants to know the value in a particular year, e.g. a snapshot of a characteristic.
        if 'year' in options['objectives'][objective_label]:
            objective += results.getValuesAt(label=objective_label, year_init=options['objectives'][objective_label]['year'])[0][0] * weight
        # Handle the default objective that wants to know a cumulative value. Timestep scaling is done inside Results.getValueAt().
        else:
            objective += results.getValuesAt(label=objective_label, year_init=options['progs_start'], year_end=settings.tvec_end, integrated=True)[0][0] * weight

    previous_results[tuple(alloc)] = objective

    return objective

def _complete_options(options, progset, settings):
    # Fills in options['init_alloc'] and fills in options['constraints']['limits'] making them
    # all absolute. Both of these dicts are odicts so that can be mapped to arrays directly
    #
    # This should be done once, prior to doing anything

    options = dcp(options)

    # If already completed (because it was run in parallelOptimizeFunc) return early
    if 'completed' in options and options['completed']:
        return options

    init_alloc = odict()
    limits = odict()

    # Add default budgets to ['init_alloc'] if required
    for prog in progset.progs:
        if prog.label not in options['init_alloc'] and 'saturate_with_default_budgets' in options and options['saturate_with_default_budgets'] is True:
            init_alloc[prog.label] = prog.getDefaultBudget()
        else:
            init_alloc[prog.label] = options['init_alloc'][prog.label]

        # Now process constraints - note checking whether the program is in 'init_alloc' guarantees that
        # every program has a constraint and there are no superfluous constraints
        if prog.label in init_alloc:

            # Default constraints
            if prog.func_specs['type'] == 'cost_only':
                limits[prog.label] = {'vals': (init_alloc[prog.label], init_alloc[prog.label]), 'rel': False}  # Upper and lower bounds equal current value i.e. unoptimizable
            else:
                limits[prog.label] = {'vals': (0.0,np.inf), 'rel': False} # Default is unlimited range

            # If previous constraints were specified
            if prog.label in options['constraints']['limits']:
                this_limit = options['constraints']['limits'][prog.label]
                if this_limit['rel']:
                    limits[prog.label]['vals'] = (init_alloc[prog.label]*this_limit['vals'][0], np.inf if np.isinf(this_limit['vals'][1]) else init_alloc[prog.label]*this_limit['vals'][1])
                else:
                    limits[prog.label]['vals'] = tuple(this_limit['vals'])

    options['init_alloc'] = init_alloc
    options['constraints']['limits'] = limits

    # If user has not supplied constraints, then make the optimisation a redistribution of default budgets.
    # If the user specifies options['constraints']['total']=False then constrainAllocation will skip it
    if 'total' not in options['constraints']:
        options['constraints']['total'] = sum(options['init_alloc'].values())

    if 'objectives' not in options:
        options['objectives'] = {settings.charac_pop_count : {'weight':-1, 'year':2030.0}}

    if 'progs_start' not in options:
        options['progs_start'] = 2015.0
        logger.info("Programs will start overwriting calibrated parameter values in the year: %f" % options['progs_start'])

    options['completed'] = True
    return options

def optimizeFunc(settings, parset, progset, options=None, outputqueue=None, thread=None, randseed=None, **optimization_params):
    # The options dict critically but optionally specifies
    # - init_alloc : the x0 values for ASD per program, defaulting depends on saturate_with_default_budgets
    # - orig_alloc : if not specified, it's the same as init_alloc. This is the very first allocation that the constraints are relative to

    if options is None:
        logger.info("An options dictionary was not supplied for optimisation. A default one will be constructed.")
        options = defaultOptimOptions(settings=settings, progset=progset)

    options = _complete_options(options, progset, settings)

    # Assemble the arrays
    x0 = [x for x in options['init_alloc'].values()]
    x0 = dcp(_constrainAllocation(alloc=x0, settings=settings, options=options))

    # Set limits
    optimization_params['xmin'] = [x['vals'][0] for x in options['constraints']['limits'].values()]
    optimization_params['xmax'] = [x['vals'][1] for x in options['constraints']['limits'].values()]

    args = {'settings':settings,
            'parset':parset,
            'progset':progset,
            'options':options,
            'previous_results':dict()}

    alloc_new, obj_vals, exit_reason = asd(_calculateObjective, x0, args=args, randseed=randseed, **optimization_params)
    alloc_new = dcp(_constrainAllocation(alloc=alloc_new, settings=settings, options=options))

    # Makes sure allocation is returned in dictionary format.
    alloc_opt = odict()
    for i,prog_key in enumerate(options['init_alloc'].keys()):
        alloc_opt[prog_key] = alloc_new[i]

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

    options = _complete_options(options, progset, settings)

    if block_iter is None and max_iters is not None and max_blocks is not None:
        logging.info("Parallel optimization: block_iter was not specified; calculating it from max_ters and max_blocks")
        block_iter = np.floor(max_iters / float(max_blocks)) # Calculate block_iter from max_iter, if supplied

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
            thisbestval = outputlist[i][1][-1]
            if thisbestval < bestfvalval:
                bestfvalval = thisbestval
                bestfvalind = i
            tmp_fval = outputlist[i][1][-block_iter:]
            if len(tmp_fval) < block_iter:
                # we could use other values, i.e. nan, but if a thread ends because it has stopped decreasing
                # then we should use the last known value
                logging.warn("Padding feval output with last value")
                tmp_fval = np.concatenate((tmp_fval , np.array([thisbestval] * (block_iter - len(tmp_fval)))))
            fvalarray[i, block * block_iter:(block + 1) * block_iter] = tmp_fval # in case the first timestep is included



        options['init_alloc'] = outputlist[bestfvalind][0] # Update the budget and use it as the input for the next block -- this is key!

        if doplot:
            clf()
            plot(np.transpose(fvalarray))

    results = (options['init_alloc'], fvalarray[bestfvalind, :], 'Parallel optimization')

    if fullfval: results[1] = fvalarray # Use the full fvalarray instead of just the best one

    return results
