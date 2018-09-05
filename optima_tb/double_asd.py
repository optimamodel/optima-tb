import logging
logger = logging.getLogger(__name__)

def asd(function, x, args=None, stepsize=0.1, sinc=2, sdec=2, pinc=2, pdec=2,
    pinitial=None, sinitial=None, absinitial=None, xmin=None, xmax=None,
    maxiters=None, maxtime=None, abstol=None, reltol=1e-3, stalliters=None,
    stoppingfunc=None, randseed=None, label=None, fulloutput=True, verbose=None, 
    minimprove=0., minstep=1e-6, numparmod = 2, **kwargs):
    """
    Optimization using adaptive stochastic descent (ASD).
    
    output = asd(func,x0) starts at x0 and attempts to find a 
    local minimizer x of the function func. func accepts input x and returns a scalar 
    function value evaluated  at x. x0 can be a scalar, list, or Numpy array of 
    any size. 
    
    If fulloutput is False, then asd() returns x only. If it is true, then it returns
    a tuple with the following items:
        x       -- The parameter set that minimizes the objective function
        fval    -- The value of the objective function at each iteration (use fval[-1] for final)
        details -- A dictionary with additional output: exitreason, why the algorithm terminated;
                   probabilities, the probability of each step; and stepsizes, the size of each
                   step for each parameter
    
    asd() has the following options that can be set using keyword arguments. Their
    names and default values are as follows:
    
      stepsize       0.1     Initial step size as a fraction of each parameter
      sinc           2       Step size learning rate (increase)
      sdec           2       Step size learning rate (decrease)
      pinc           2       Parameter selection learning rate (increase)
      pdec           2       Parameter selection learning rate (decrease)
      pinitial       None    Set initial parameter selection probabilities
      sinitial       None    Set initial step sizes; if empty, calculated from stepsize instead
      minimprove     0       Any improvement in the objective function less than this will be treated as no improvement
      minstep        1e-6    Stop adjusting a parameter once the step size drops below this value
      numparmod      1       Number of parameters to modify each step (recommended not high!)
      xmin           None    Min value allowed for each parameter  
      xmax           None    Max value allowed for each parameter 
      maxiters       1000    Maximum number of iterations (1 iteration = 1 function evaluation)
      maxtime        3600    Maximum time allowed, in seconds
      abstol         1e-6    Minimum absolute change in objective function
      reltol         5e-3    Minimum relative change in objective function
      stalliters     50      Number of iterations over which to calculate TolFun
      stoppingfunc   None    External method that can be used to stop the calculation from the outside.
      randseed       None    The random seed to use
      fulloutput     True    Whether or not to return the full output
      verbose        2       How much information to print during the run
      label          None    A label to use to annotate the output
  
    Example:
        from asd import asd
        from numpy.linalg import norm
        x, fval, details = asd(norm, [1, 2, 3])
    
    Version: 2017may17 by Cliff Kerr (cliff@thekerrlab.com)
    """
#     logging.info("asd.py: Unused kwargs", kwargs)

    from numpy import array, shape, reshape, ones, zeros, mean, cumsum, mod, concatenate, floor, flatnonzero, isnan, inf
    from numpy.random import random, seed
    from copy import deepcopy # For arrays, even y = x[:] doesn't copy properly
    from time import time
    if randseed is not None:
        seed(int(randseed)) # Don't reset it if not supplied
        logger.debug('Launching ASD with random seed is %i; sample: %f' % (randseed, random()))

    def consistentshape(userinput, origshape=False):
        """
        Make sure inputs have the right shape and data type.
        """
        output = reshape(array(userinput, dtype='float'), -1)
        if origshape: return output, shape(userinput)
        else:         return output

    # Handle inputs and set defaults
    if maxtime  is None: maxtime = 3600
    if maxiters is None: maxiters = 1000
    maxrangeiters = 1000 # Number of times to try generating a new parameter
    x, origshape = consistentshape(x, origshape=True) # Turn it into a vector but keep the original shape (not necessarily class, though)
    nparams = len(x) # Number of parameters

    nps = range(numparmod)  #The number of things parameters that may need to be changed each time step.

    # Set initial parameter selection probabilities -- uniform by default
    if pinitial is None: probabilities = ones(2 * nparams)
    else:                probabilities = consistentshape(pinitial)

    # Handle step sizes
    if sinitial is None:
        stepsizes = abs(stepsize * x)
        stepsizes = concatenate((stepsizes, stepsizes)) # need to duplicate since two for each parameter
    else:
        stepsizes = consistentshape(sinitial)

    # Handle x limits
    xmin = zeros(nparams) - inf if xmin is None else consistentshape(xmin)
    xmax = zeros(nparams) + inf if xmax is None else consistentshape(xmax)

    # Final input checking
    if sum(isnan(x)):
        errormsg = 'At least one value in the vector of starting points is NaN:\n%s' % x
        raise Exception(errormsg)
    if label is None: label = ''
    if stalliters is None: stalliters = 10 * nparams # By default, try 10 times per parameter on average
    stalliters = int(stalliters)
    maxiters = int(maxiters)

    # Initialization
    if all(stepsizes == 0): stepsizes += stepsize # Handle the case where all step sizes are 0
    if any(stepsizes == 0): stepsizes[stepsizes == 0] = mean(stepsizes[stepsizes != 0]) # Replace step sizes of zeros with the mean of non-zero entries
    if args is None: args = {} # Reset if no function arguments supplied
    fval = function(x, **args) # Calculate initial value of the objective function
    fvalorig = fval # Store the original value of the objective function, since fval is overwritten on each step
    xorig = deepcopy(x) # Keep the original x, just in case

    # Initialize history
    abserrorhistory = zeros(stalliters) # Store previous error changes
    relerrorhistory = zeros(stalliters) # Store previous error changes
    fvals = zeros(maxiters + 1) # Store all objective function values
    allsteps = zeros((maxiters + 1, nparams)) # Store all parameters
    fvals[0] = fvalorig # Store initial function output
    allsteps[0, :] = xorig # Store initial input vector

    # Loop
    count = 0 # Keep track of how many iterations have occurred
    start = time() # Keep track of when we begin looping
    offset = ' ' * 4 # Offset the print statements
    exitreason = 'Unknown exit reason' # Catch everything else
    while True:
        count += 1 # Increment the count

        # Calculate next parameters
        if sum(probabilities) == 0:
            exitreason = 'Ran out of parameters to select'
            break

        probabilities = probabilities / sum(probabilities) # Normalize probabilities
        cumprobs = cumsum(probabilities) # Calculate the cumulative distribution
        inrange = False
        for r in range(maxrangeiters): # Try to find parameters within range
            choices = [flatnonzero(cumprobs > random())[0] for _ in nps] # Choose a parameter and upper/lower at random
            pars = mod(choices, nparams) # Which parameter was chosen
            pms = list(floor(array((choices)) / nparams)) # Plus or minus
            newvals = [x[pars[n]] + ((-1) ** pms[n]) * stepsizes[choices[n]] for n in nps] # Calculate the new vector
            if array([newvals[n] >= xmin[pars[n]] and newvals[n] <= xmax[pars[n]] for n in nps]).all(): # Make sure it's in range
                inrange = True
                break
            else:
                for n in nps:
                    if stepsizes[choices[n]] < minstep:
                        probabilities[choices[n]] = 0.0 # We have decreased the stepsize too far for this parameter to matter, stop changing it
                else:
                    stepsizes[choices[n]] = stepsizes[choices[n]] / sdec # Decrease size of step for next time


        if not inrange:
            logger.warning('======== Can\'t find parameters within range after %i tries, terminating ========' % maxrangeiters)
            break

        # Calculate the new value
        xnew = deepcopy(x) # Initialize the new parameter set
        for n in nps:
            xnew[pars[n]] = newvals[n] # Update the new parameter set
        fvalnew = function(xnew, **args) # Calculate the objective function for the new parameter set
        abserrorhistory[mod(count, stalliters)] = max(0, fval - fvalnew) # Keep track of improvements in the error
        relerrorhistory[mod(count, stalliters)] = max(0, fval / float(fvalnew) - 1.0) # Keep track of improvements in the error
        logger.info(offset + 'step=%i choices=%s, pars=%s, pms=%s, origval=%s, newval=%s, inrange=%s' % (count, choices, pars, pms, [x[par] for par in pars], [xnew[par] for par in pars], inrange))

        # Check if this step was an improvement
        fvalold = fval # Store old fval
        if fvalnew < fvalold - minimprove: # New parameter set is better than previous one
            for n in nps:  #assign the desire to do more to all parameters changed
                probabilities[choices[n]] = probabilities[choices[n]] * pinc # Increase probability of picking this parameter again
                stepsizes[choices[n]] = stepsizes[choices[n]] * sinc # Increase size of step for next time
            x = xnew # Reset current parameters
            fval = fvalnew # Reset current error
            flag = '++' # Marks an improvement
        elif fvalnew >= fvalold - minimprove: # New parameter set is the same or worse than the previous one
            for n in nps:            #assign the desire to do less to all parameters changed    
                probabilities[choices[n]] = probabilities[choices[n]] / pdec # Decrease probability of picking this parameter again
                stepsizes[choices[n]] = stepsizes[choices[n]] / sdec # Decrease size of step for next time
            flag = '--' # Marks no change
        else:
            exitreason = 'Objective function returned NaN'
            break
        logger.info(offset + label + ' step %i (%0.1f s) %s (orig: %s | best:%s | new:%s | diff:%s)' % ((count, time() - start, flag) + multisigfig([fvalorig, fvalold, fvalnew, fvalnew - fvalold])))

        # Store output information
        fvals[count] = fval # Store objective function evaluations
        allsteps[count, :] = x # Store parameters

        # Stopping criteria
        if count >= maxiters: # Stop if the iteration limit is exceeded
            exitreason = 'Maximum iterations reached'
            break
        if (time() - start) > maxtime:
            exitreason = 'Time limit reached (%s > %s)' % multisigfig([(time() - start), maxtime])
            break
        if (count > stalliters) and (abs(mean(abserrorhistory)) < abstol): # Stop if improvement is too small
            exitreason = 'Absolute improvement too small (%s < %s)' % multisigfig([mean(abserrorhistory), abstol])
            break
        if (count > stalliters) and (sum(relerrorhistory) < reltol): # Stop if improvement is too small
            exitreason = 'Relative improvement too small (%s < %s)' % multisigfig([mean(relerrorhistory), reltol])
            break
        if stoppingfunc and stoppingfunc():
            exitreason = 'Stopping function called'
            break

    # Return
    x = reshape(x, origshape) # Parameters
    fvals = fvals[:count + 1] # Function evaluations
    logger.info('=== %s %s (%i steps, orig: %s | best: %s | ratio: %s) ===' % ((label, exitreason, count) + multisigfig([fvals[0], fvals[-1], fvals[-1] / fvals[0]])))
    if fulloutput:
        details = dict()
        details['exitreason'] = exitreason
        details['probabilities'] = probabilities
        details['stepsizes'] = stepsizes
        details['allsteps'] = allsteps[:count + 1, :]
        return (x, fvals, details)
    else:
        return (x, None, None) # so that it has the same shape regardless of fulloutput




def multisigfig(X, sigfigs=5):
    """ Return a string representation of variable x with sigfigs number of significant figures """
    from numpy import log10, floor

    output = []
    try:
        n = len(X)
        islist = True
    except:
        X = [X]
        n = 1
        islist = False
    for i in range(n):
        x = X[i]
        try:
            if x == 0:
                output.append('0')
            else:
                magnitude = floor(log10(abs(x)))
                factor = 10 ** (sigfigs - magnitude - 1)
                x = round(x * factor) / float(factor)
                digits = int(abs(magnitude) + max(0, sigfigs - max(0, magnitude) - 1) + 1 + (x < 0) + (abs(x) < 1)) # one because, one for decimal, one for minus
                decimals = int(max(0, -magnitude + sigfigs - 1))
                strformat = '%' + '%i.%i' % (digits, decimals) + 'f'
                string = strformat % x
                output.append(string)
        except:
            output.append(str(x))
    if islist:
        return tuple(output)
    else:
        return output[0]
