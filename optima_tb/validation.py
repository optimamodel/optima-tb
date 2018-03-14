#%% Imports
from utils import flattenDict, odict, OptimaException
from parsing import FunctionParser
import optima_tb.settings as settings

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

def checkTransitionFraction(transition,validationSettings):
    """
    
    Params:
        transition            float
        validationSetting     optima_tb.settings.ValidationSettings object
    
    """
    validation_level = validationSettings['transition_fraction']
    if transition > 1.: 
        warning = "Transition value observed larger than 1.: %g"%transition
        if validation_level == settings.VALIDATION_ERROR:
            raise OptimaException(warning)
        elif validation_level == settings.VALIDATION_WARN:
            logger.warn(warning)
        transition = 1.
    return transition


def checkSummedOutflowRate():
    """
    This validation check should be performed after building the model, rather than during the simulation.
    For each compartment, it checks that the summed outflow is valid i.e. fractions 
    
    Note that in practice, the inflow to a compartment could be sufficiently large to allow summed 
    outflow rates that are greater than 1.
    For this reason, this check should NEVER be run as an AVERT or ERROR level, but only as IGNORE or WARN.
    
    """
    # TODO: implement checkSummedOutflowRate

    pass




def checkInitialization():
    """
    This validation check should be run after loading data in from a country's databook.
    It checks that all states have been initialized and should either WARN, ERROR or IGNORE if not.
    
    """
    # TODO: implement checkInitialization

    pass


def checkPopulationExplosion():
    """
    Validation to check on whether our populations have exploded in size or not, beyond some tolerance epsilon. This is due to the fact
    that initial populations are used as well as a parameter transitions, but that the population size after this is not enforced or checked
    to determine how the population grows. 
    
    For validation settings for ERROR or WARN: if the summed population change for Population P between time t and t+1 is beyond
        some tolerance value epsilon*, either an error is raised or print a warning.
        *Note that this change is in addition to aging, other migrations or interpopulation transfers.
    
    For validation settings for AVERT: if the summed population change for Population P between time t and t+1 is beyond
        some tolerance value epsilon*, then the values for each compartment are normalized so that the sum total is P(t) +/- epsilon.
        *Note that this change is in addition to aging, other migrations or interpopulation transfers.
    
    For validation settings for IGNORE: no check is performed
    
    """
    # TODO: implement checkPopulationExplosion
    pass




    