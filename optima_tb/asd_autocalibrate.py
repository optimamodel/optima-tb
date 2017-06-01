import logging
logger = logging.getLogger(__name__)

from optima_tb.utils import OptimaException



def asd(function, init_params, init_compartments=[], args=None, stepsize=0.1, xmin=None, xmax=None, xnames=None,
        sinc=2, sdec=2, pinc=2, pdec=2, pinitial=None, sinitial=None, absinitial=None, MaxRangeIter=1000,
        MaxFunEvals=None, MaxIter=1e3, AbsTolFun=1e-6, RelTolFun=1e-2, TolX=None, StallIterLimit=100,
        fulloutput=False, maxarraysize=1e6, timelimit=3600, stoppingfunc=None, randseed=None,useYFactor=False,**kwargs):
    """
    Optimization using the adaptive stochastic descent algorithm.
    
    X, FVAL, EXITFLAG, OUTPUT = asd(FUN,X0) starts at X0 and attempts to find a 
    local minimizer X of the function FUN. FUN accepts input X and returns a scalar 
    function value F evaluated  at X. X0 can be a scalar, list, or Numpy array of 
    any size. The outputs are:
        X        -- The parameter set that minimizes the objective function
        FVAL     -- The value of the objective function at X
        EXITFLAG -- The exit condition of the algorithm possibilities are:
                     0 -- Maximum number of function evaluations or iterations reached
                     1 -- Step size below threshold
                     2 -- Improvement in objective function below minimum threshold
                     3 -- Maximum number of iterations to calculate new parameter when out of range reached
                     4 -- Time limit exceeded
                     5 -- Stopping function criteria met
                    -1 -- Algorithm terminated for other reasons
          OUTPUT -- An object with the following attributes:
            iterations -- Number of iterations
            funcCount  -- Number of function evaluations
            fval       -- Value of objective function at each iteration
            x          -- Vector of parameter values at each iteration
    
    asd() has the following options that can be set using keyword arguments. Their
    names and default values are as follows:
    
    # TODO : update documentation <-----------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
      stepsize       {0.1}      -- Initial step size as a fraction of each parameter
      xmin           {[]}       -- Min value allowed for each parameter  
      xmax           {[]}       -- Max value allowed for each parameter 
      sinc           {2}        -- Step size learning rate (increase)
      sdec           {2}        -- Step size learning rate (decrease)
      pinc           {2}        -- Parameter selection learning rate (increase)
      pdec           {2}        -- Parameter selection learning rate (decrease)
      pinitial       {ones(2N)} -- Set initial parameter selection probabilities
      sinitial       {[]}       -- Set initial step sizes; if empty, calculated from stepsize instead
      MaxRangeIter   {1000}     -- Maximum number of iterations to calculate new parameter when out of range
      MaxFunEvals    {N*1e3}    -- Maximum number of function evaluations
      MaxIter        {1e3}      -- Maximum number of iterations (1 iteration = 1 function evaluation)
      AbsTolFun      {1e-3}     -- Minimum absolute change in objective function
      RelTolFun      {1e-2}     -- Minimum relative change in objective function
      TolX           {N*1e-6}   -- Minimum change in parameters
      StallIterLimit {100}      -- Number of iterations over which to calculate TolFun
      fulloutput     {True}     -- Whether or not to output the parameters and errors at each iteration
      maxarraysize   {1e6}      -- Limit on MaxIter and StallIterLimit to ensure arrays don't get too big
      timelimit      {3600}     -- Maximum time allowed, in seconds
      stoppingfunc   {None}     -- External method that can be used to stop the calculation from the outside.
      randseed       {None}     -- The random seed to use
  
    
    Example:
        from asd import asd
        from numpy import norm
        x, fval, exitflag, output = asd(norm, [1, 2, 3])
    
    
    Version: 2016feb11 by Cliff Kerr (cliff@thekerrlab.com) in optimamodel / optima-HIV
    Modified:   2016dec05 by Sarah Jarvis, when added to Optima-TB
    """
    
    from numpy import array, shape, reshape, ones, zeros, size, mean, cumsum, mod, hstack, floor, flatnonzero, isnan, sum
    from numpy.random import random, seed
    from copy import deepcopy # For arrays, even y = x[:] doesn't copy properly
    from time import time
    seed(randseed)
    
    # Initialization of required variables
    x = init_params+init_compartments
    n_init_params = len(init_params)
    nparams = len(x)
    logger.debug("nparams for asd = %g"%nparams)
    p = ones(2*nparams)  # Set initial parameter selection probabilities -- uniform by default
    p = p/sum(p) # Normalize probabilities
    steps = ones(2*nparams)*stepsize
    #print "steps = ", steps
    #print "pi = ", p
    # helper variables
    MaxFunEvals = 1000*nparams if MaxFunEvals == None else MaxFunEvals # Maximum number of function evaluations
    TolX = 1e-6 if TolX == None else TolX  # Minimum change in parameters
    MaxIter = min(MaxIter, maxarraysize);
    StallIterLimit = min(StallIterLimit, maxarraysize); # Don't by default let users create arrays larger than this -- slow and pointless
    
        
    # Setup params related to first pass of evaluating parameters
    fval = sum(function(x,init_compartments))# Calculate initial value of the objective function
    fval = fval.sum()
    fvalorig = fval # Store the original value of the objective function, since fval is overwritten on each step
    logger.info("Initial score for fit : %g"%fvalorig)
    
    count = 0 # Keep track of how many iterations have occurred
    
    exitflag = -1 # Set default exit flag
    abserrorhistory = zeros(int(StallIterLimit)) # Store previous error changes
    relerrorhistory = zeros(int(StallIterLimit)) # Store previous error changes
    
    if fulloutput: # Include additional output structure
        fulloutputfval = zeros(int(MaxIter)) # Store all objective function values
        fulloutputx = zeros((int(MaxIter),int(nparams))) # Store all parameters
    
    ## Loop
    start = time()
    offset = ' '*4 # Offset the print statements
    
    xlabel = ''
    
    while True:
        logger.info(offset+'Iteration %i; elapsed %0.1f s; objective: --' % (count+1, time()-start))
        # Calculate next step
        count += 1 # On each iteration there are two function evaluations
        p = p/sum(p) # Normalize probabilities
        cumprobs = cumsum(p) # Calculate the cumulative distribution
        inrange = 0
        inner_count = 0
        while not inrange:
            inner_count += 1
            choice = flatnonzero(cumprobs > random())[0] # Choose a parameter and upper/lower at random
            par = mod(choice,nparams) # Which parameter was chosen
            

            
            pm = floor((choice)/nparams) # Plus or minus
            newval = x[par] + ((-1)**pm)*steps[choice] # Calculate the new parameter set
            """    
            if par < n_init_params:
                newval = x[par] + ((-1)**pm)*steps[choice] # Calculate the new parameter set
            else: # it's a compartment value
                newval = x[par]
                if x[par] == 0:
                    newval = 1.
                newval *= (1+((-1)**pm)*steps[choice] )
            """
            if inner_count > MaxRangeIter: # f stuck due to x range limits, exit after 1000 iterations
                newval = x[par]
                #exitflag = -1
                inrange = 1
            elif (xmax is None) and (xmin is None):
                inrange = 1
            elif (xmin is None) and (xmax is not None) and (newval <= xmax[par]):
                inrange = 1
            elif (xmax is None) and (xmin is not None) and (newval >= xmin[par]):
                inrange = 1
            elif (xmax is not None) and (xmin is not None) and (xmax[par] is None) and (newval >= xmin[par]):
                inrange = 1
            elif (xmax is not None) and (xmin is not None) and (newval <= xmax[par]) and (newval >= xmin[par]):
                inrange = 1
            else:
                p[choice] = p[choice]/pdec # decrease probability of picking this parameter again
                steps[choice] = steps[choice]/sdec # decrease size of step for next time

        # Set up copies
        xnew = deepcopy(x) # Initialize the new parameter set
        xnew[par] = newval # Update the new parameter set
        param_new = xnew[:n_init_params]
        compartments_new = xnew[n_init_params:]
        
        # TODO: use updated initcompartments
        # Take the parameter set x and run it on the model:
        try:
            fvalnew = function(param_new,compartments_new) # Calculate the objective function for the new parameter set
            fvalnew = sum(fvalnew)
            fvalnew = fvalnew.sum()
        except OptimaException as optima_error:
            logger.debug("Encountered an OptimaException when running autofit:\n%s\n"%(str(optima_error)))
            continue
        
        abserrorhistory[mod(count,StallIterLimit)] = max(0, fval-fvalnew) # Keep track of improvements in the error
        relerrorhistory[mod(count,StallIterLimit)] = max(0, fval/float(fvalnew)-1.0) # Keep track of improvements in the error  
        if xnames is not None: 
            xlabel = ' (xlabel=%s)'%xnames[par]
        
        logger.info(offset+'step=%i choice=%s%s, par=%s, pm=%s, origval=%s, newval=%s, inrange=%s' % (count, choice, xlabel,par, pm, x[par], xnew[par], inrange))

        # Check if this step was an improvement
        fvalold = sum(fval) # Store old fval
        
        if fvalnew < fvalold: # New parameter set is better than previous one
            p[choice] = p[choice]*pinc # Increase probability of picking this parameter again
            steps[choice] = steps[choice]*sinc # Increase size of step for next time
            x = xnew # Reset current parameters
            fval = fvalnew # Reset current error
            flag = 'Improvement'
        elif fvalnew >= fvalold: # New parameter set is the same or worse than the previous one
            p[choice] = p[choice]/pdec # Decrease probability of picking this parameter again
            steps[choice] = steps[choice]/sdec # Decrease size of step for next time
            flag = '  No change'
        else:
            exitflag = -1
            logger.info('======== Objective function returned NaN, terminating ========')
            break
        logger.info(offset + 'Step %i (%0.1f s): %s (orig: %s | best:%s | new:%s | diff:%s | ratio:%0.5f)' % ((count, time()-start, flag)+multisigfig([fvalorig, fvalold, fvalnew, fvalnew-fvalold]) + (fvalnew/fvalold,)))
        
        # Optionally store output information
        if fulloutput: # Include additional output structure
            fulloutputfval[count-1] = fval # Store objective function evaluations
            fulloutputx[count-1,:] = x # Store parameters
        
        # Stopping criteria
        if (count+1) >= MaxFunEvals: # Stop if the function evaluation limit is exceeded
            exitflag = 0 
            logger.info('======== Maximum function evaluations reached (%i >= %i), terminating ========' % ((count+1), MaxFunEvals))
            break
        if count >= MaxIter: # Stop if the iteration limit is exceeded
            exitflag = 0 
            logger.info('======== Maximum iterations reached (%i >= %i), terminating ========' % (count, MaxIter))
            break 
        if mean(steps) < TolX: # Stop if the step sizes are too small
            exitflag = 1 
            logger.info('======== Step sizes too small (%f < %f), terminating ========' % (mean(steps), TolX))
            break
        if (count > StallIterLimit) and (abs(mean(abserrorhistory)) < AbsTolFun): # Stop if improvement is too small
            exitflag = 2 
            logger.info('======== Absolute improvement too small (%f < %f), terminating ========' % (mean(abserrorhistory), AbsTolFun))
            break
        if (count > StallIterLimit) and (mean(relerrorhistory) < (RelTolFun/StallIterLimit)): # Stop if improvement is too small
            exitflag = 2 
            logger.info('======== Relative improvement too small (%f < %f), terminating ========' % (mean(relerrorhistory), RelTolFun))
            break
        if inner_count > MaxRangeIter: 
            exitflag = 3
            logger.info('======== Can\'t find parameters within range (%i > %i), terminating ========' % (inner_count, MaxRangeIter))
            break
        if timelimit is not None and (time()-start)>timelimit:
            exitflag = 4
            logger.info('======== Time limit reached (%f > %f), terminating ========' % ((time()-start), timelimit))
            break
        if stoppingfunc and stoppingfunc():
            exitflag = 5
            logger.info('======== Stopping function called, terminating ========')
            break

    # Create additional output
    class makeoutput:
        iterations = count # Number of iterations
        funcCount = count+1 # Number of function evaluations
        if fulloutput: # Include additional output structure
            fval = fulloutputfval[:count] # Function evaluations
            x = fulloutputx[:count,:] # Parameters
            
    output = makeoutput()
    
    """
    # For debugging: check on variables 
    print p
    print steps
    print abserrorhistory
    print relerrorhistory
    """
    return x, fval, exitflag, output




def multisigfig(X, sigfigs=5):
    """ Return a string representation of variable x with sigfigs number of significant figures """
    from numpy import log10, floor
    
    output = []
    try: 
        n=len(X)
        islist = True
    except:
        X = [X]
        n = 1
        islist = False
    for i in range(n):
        x = X[i]
        try:
            if x==0:
                output.append('0')
            else:
                magnitude = floor(log10(abs(x)))
                factor = 10**(sigfigs-magnitude-1)
                x = round(x*factor)/float(factor)
                digits = int(abs(magnitude) + max(0, sigfigs - max(0,magnitude) - 1) + 1 + (x<0) + (abs(x)<1)) # one because, one for decimal, one for minus
                decimals = int(max(0,-magnitude+sigfigs-1))
                strformat = '%' + '%i.%i' % (digits, decimals)  + 'f'
                string = strformat % x
                output.append(string)
        except:
            output.append(str(x))
    if islist:
        return tuple(output)
    else:
        return output[0]