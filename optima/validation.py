#%% Imports
from utils import flattenDict, odict, OptimaException
from parsing import FunctionParser
import settings

import logging
logger = logging.getLogger(__name__)

import numpy as np



def isDPopValid(model,msettings,dpopsizes,ti,validationSetting):
    
    logging.info("(t=%g) Checking dpop update is valid "%ti)
    if validationSetting == settings.VALIDATION_AVERT:
        pops = model.pops
        num_pops = len(pops)
        num_comps = len(pops[0].comps) 
    
        for pid in xrange(num_pops):
            for cid in xrange(num_comps):
                did = pid * num_comps + cid
                
                if pops[pid].comps[cid].popsize[ti+1] < 0 :
                    return False
                # check dpop won't result in a negative population compartment
                if pops[pid].comps[cid].popsize[ti+1] + dpopsizes[did] < 0. : 
                    #warning = "population=%g, cid=%g, did=%g, t=%g : value = %g, dpop = %g"%(pid,cid,did,ti,pops[pid].comps[cid].popsize[ti+1],dpopsizes[did])
                    #logging.info("Invalid dpop proposed for %s"%warning)
                    return False
                if pops[pid].comps[cid].popsize[ti+1] == 0. and dpopsizes[did] < 0.:
                    return False
        return True
    
    # for all other settings, we'll either ignore, warn or throw an error, so we don't need to check
    return True



def checkNegativePopulation(model,msettings,dpopsizes,dpop_out,ti,dt,validationSetting):
    """ 
    This validation check is run during the simulation, and either warns, stops or averts instances when 
    the calculated outflow of compartment is greater than the current value of that compartment.
    
    For validation settings for ERROR or WARN: the offending compartment is noted and either an error or
        warning is raised, respectively. 
    
    For validation settings for AVERT: the compartment for a population that has an outflow larger than 
        the current population is noted, and the ratio for rescaling is noted. This is passed back to 
        model._calculateDPops to be applied when calculating the new dpopsizes, as we need to conserve values
        which requires that not only do we change the outflow of a compartment, but that the inflow to the next compartment
        is also updated too.
        
        Note that the rescaling factor is calculated as the (current compartment size) / (the observed outwards dpop amount). 
        Initially, we tried (current compartment size) / (the observed total dpop amount), but this led to errors due to the inflow
        component of dpop . 
    
    For validation settings for IGNORE: no check is performed
    
    """
    pops = model.pops
    num_pops = len(pops)
    num_comps = len(pops[0].comps) 
    
    logging.info("(t=%g) Checking validation of negative populations "%ti)
    
    if validationSetting == settings.VALIDATION_ERROR or validationSetting == settings.VALIDATION_WARN:
        # go through each compartment and check that the value of the corresponding dpop is less than the current compartment's value
        for pid in xrange(num_pops):
            for cid in xrange(num_comps):
                did = pid * num_comps + cid
                #logging.debug("%g %g %g"%(did, pops[pid].comps[cid].popsize[ti+1] , dpopsizes[did]))
                
                # if not, then either thrown an OptimaException or warn, as required
                if pops[pid].comps[cid].popsize[ti+1] + dpopsizes[did] < 0.:
                    
                    warning = "Negative value encountered for: population=%g, cid=%g, t=%g : value = %g, dpop = %g"%(pid,cid,ti,pops[pid].comps[cid].popsize[ti+1],dpopsizes[did])
                    
                    if validationSetting == settings.VALIDATION_ERROR:
                        raise OptimaException("ERROR: "+warning)
                    else:
                        logging.warn(warning)
                    
    
    elif validationSetting == settings.VALIDATION_AVERT:
        
        
        empty_compartment = []
        reduced_compartment = {}
        reset_compartment = []
        
        for pid in xrange(num_pops):
            for cid in xrange(num_comps):
                did = pid * num_comps + cid
    
                warning = "... encountered for: population=%g, cid=%g, did=%g, t=%g : value = %g, dpop = %g"%(pid,cid,did,ti,pops[pid].comps[cid].popsize[ti+1],dpopsizes[did])
                    
                if pops[pid].comps[cid].popsize[ti+1] == 0. and dpopsizes[did] < 0.:
                    logging.debug("--- Compartment is empty; cannot remove from empty compartment: %s"%warning)
                    empty_compartment.append(did)
                    
                elif pops[pid].comps[cid].popsize[ti+1] + dpopsizes[did] < 0.:
                    logging.debug("Have to reduce the corresponding flow for this compartment: %s"%warning)
                    #print did, dpopsizes[did], pops[pid].comps[cid].popsize[ti+1]
                    curr_pop = pops[pid].comps[cid].popsize[ti+1]
                    reduced_compartment[did] = curr_pop / np.abs(dpop_out[did])
                    
                elif pops[pid].comps[cid].popsize[ti+1] < 0 :
                    # TODO: modify appropriately
                    logging.debug("--- Compartment value is negative: %s"%warning)
                    reset_compartment.append(did)
                    
    
        
        # go back and recalculate links so that 
        # - there is no corresponding inflow from an empty compartments 
        # - for reduced compartments, all contributing amounts that contribute from this compartment are proportioned
        #   based on the theoretical sum and the available amount
        
        dp, dp_out =  model._calculateDPops(msettings,ti,dt,reduced_compartment,empty_compartment,reset_compartment)
        
        """
        print empty_compartment , reduced_compartment
        print "changes for reduced compartments"
        for k in reduced_compartment.keys():
            print "k: old dpop    pop value         new dpop value"
            print k, dpopsizes[k] , pops[k/num_comps].comps[k%num_comps].popsize[ti+1], dp[k]
        
        print "changes for empty compartments"
        for k in empty_compartment:
            print "k: old dpop    pop value         new dpop value"
            print k , dpopsizes[k] , pops[k/num_comps].comps[k%num_comps].popsize[ti+1], dp[k]
        """
            
        return dp, dp_out, reset_compartment
        
    elif validationSetting == settings.VALIDATION_IGNORE:
        # do nothing
        pass
    
    else:
        logging.warn("Unknown validation setting for validation.checkNegativePopulation: %g\nIgnoring validation check"%validationSetting)

    return dpopsizes, [], []





def checkSummedOutflowRate():
    """
    This validation check should be performed after building the model, rather than during the simulation.
    For each compartment, it checks that the summed outflow is valid i.e. fractions 
    
    Note that in practice, the inflow to a compartment could be sufficiently large to allow summed 
    outflow rates that are greater than 1.
    For this reason, this check should NEVER be run as an AVERT or ERROR level, but only as IGNORE or WARN.
    
    """
    pass




def checkInitialization():
    """
    This validation check should be run after loading data in from a country's databook.
    It checks that all states have been initialized and should either WARN, ERROR or IGNORE if not.
    
    """
    pass
    