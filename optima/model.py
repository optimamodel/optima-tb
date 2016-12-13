#%% Imports

from utils import flattenDict, odict, OptimaException
from validation import checkNegativePopulation
import settings as project_settings
from results import ResultSet

import logging
logger = logging.getLogger(__name__)

import numpy as np
from copy import deepcopy as dcp
from collections import Counter



#%% Abstract classes used in model

class Node(object):
    ''' Lightweight abstract class to represent one object within a network. '''
    def __init__(self, label = 'default', index = 0):
        self.label = label          # Reference name for this object.
        self.index = index          # Index to denote storage position in Model object. Format is left unspecified.
                                    # For a ModelPopulation, this is intended as an integer.
                                    # For a ModelCompartment, this is intended as a tuple with first element denoting index of ModelPopulation.
        self.num_outlinks = 0       # Tracks number of nodes this one is linked to as an initial node.
        self.outlink_ids = []       # List of indices corresponding to each outgoing link.
        
    def makeLinkTo(self, other_node, link_index, link_label):
        '''
        Link this node to another (i.e. create a transition link).
        Must provide an index that, by intention, uniquely represents where this link is positioned in its storage container.
        Must also provide the variable label for the link, in case Settings.linkpar_specs[link_label] need to be referred to.
        '''
        if not isinstance(other_node, Node):
            raise OptimaException('ERROR: Attempting to link compartment to something that is not a compartment.')
        self.num_outlinks += 1
        self.outlink_ids.append(link_index)
        return Link(object_from = self, object_to = other_node, label = link_label)

class Variable(object):
    '''
    Lightweight abstract class to store a variable array of values (presumably corresponding to an external time vector).
    Includes an attribute to describe the format of these values.
    '''
    def __init__(self, label = 'default', val = 0.0):
        self.label = label
        self.vals = np.array([float(val)])   # An abstract array of values.
        if val > 1:
            self.val_format = 'number'
        else:
            self.val_format = 'fraction'
        

class Link(Variable):
    '''
    A more involved version of Variable, representing unidirectional flow between two objects in a network.
    The values stored in this extended version of Variable refer to flow rates.
    If used in ModelPop, the Link refers to two cascade compartments within a single population.
    If used in Model, the Link refers to two distinct population groups.
    In the latter case, intended logic should transfer agents between all (non-dead) corresponding compartments.
    '''
    def __init__(self, object_from, object_to, label = 'default', val = 0.0, scale_factor = 1.):
        Variable.__init__(self, label = label, val = val)
        self.index_from = object_from.index
        self.index_to = object_to.index
        self.label_from = object_from.label
        self.label_to = object_to.label
        self.scale_factor = scale_factor

    def __repr__(self, *args, **kwargs):
        print "self.scale_factor = ", self.scale_factor, type(self.scale_factor)
        return "Link %s: from (%s, %s) to (%s, %s) with scale_factor=%g"%(self.label,self.index_from,self.label_from,self.index_to,self.label_to,self.scale_factor)

#%% Cascade compartment and population classes

class ModelCompartment(Node):
    ''' A class to wrap up data for one compartment within a cascade network. '''
    
    def __init__(self, label = 'default', index = 0, popsize = 0.0):
        Node.__init__(self, label = label, index = index)
        self.popsize = np.array([float(popsize)])   # Number of people in compartment.
        self.tag_birth = False                      # Tag for whether this compartment contains unborn people.
        self.tag_dead = False                       # Tag for whether this compartment contains dead people.
        self.junction = False
    
    def __repr__(self, *args, **kwargs):
        return "%s: %g"%(self.label, self.popsize[0])
    
    def getValue(self, ti):
        """ Get value of population at timestep ti """
        return self.popsize[ti]

class ModelPopulation(Node): 
    '''
    A class to wrap up data for one population within model.
    Each model population must contain a set of compartments with equivalent labels.
    '''
    
    def __init__(self, settings, label = 'default', index = 0):
        Node.__init__(self, label = label, index = index)
        self.comps = list()         # List of cascade compartments that this model population subdivides into.
        self.links = list()         # List of intra-population cascade transitions within this model population.
        self.deps = list()          # List of value arrays for link-dependencies (both characteristic and untagged parameters) required by model.
        self.comp_ids = dict()      # Maps label of a compartment to its position index within compartments list.
        self.link_ids = dict()      # Maps cascade transition tag to indices for all relevant transitions within links list.
        self.dep_ids = dict()       # Maps label of a relevant characteristic/parameter to its position index within dependencies list.
        self.t_index = 0            # Keeps track of array index for current timepoint data within all compartments.
        
        self.genCascade(settings = settings)    # Convert compartmental cascade into lists of compartment and link objects.
    
    def __repr__(self, *args, **kwargs):
        return "".join("%s"%self.comps)
    
    def getModelState(self,ti):
        states = [c.getValue(ti) for c in self.comps]
        return states
    
    def getComp(self, comp_label):
        ''' Allow compartments to be retrieved by label rather than index. Returns a ModelCompartment. '''
        comp_index = self.comp_ids[comp_label]
        return self.comps[comp_index]
        
    def getLinks(self, link_tag):
        ''' Allow links to be retrieved by tag rather than index. Returns a list of Links. '''
        link_index_list = self.link_ids[link_tag]
        link_list = []
        for link_index in link_index_list:
            link_list.append(self.links[link_index])
        return link_list
        
    def getDep(self, dep_label):
        ''' Allow dependencies to be retrieved by label rather than index. Returns a Variable. '''
        dep_index = self.dep_ids[dep_label]
        return self.deps[dep_index]
        
    def genCascade(self, settings):
        '''
        Generate a compartmental cascade as defined in a settings object.
        Fill out the compartment, transition and dependency lists within the model population object.
        Maintaining order as defined in a cascade workbook is crucial due to cross-referencing.
        '''
        for k, label in enumerate(settings.node_specs.keys()):
            self.comps.append(ModelCompartment(label = label, index = (self.index,k)))
            if 'tag_birth' in settings.node_specs[label].keys():
                self.comps[-1].tag_birth = True
            if 'tag_dead' in settings.node_specs[label].keys():
                self.comps[-1].tag_dead = True
            if 'junction' in settings.node_specs[label].keys():
                self.comps[-1].junction = True
            self.comp_ids[label] = k
        
        k = 0
        # Create actual link objects for parameters with tags, a.k.a. transitions.
        for label in settings.linkpar_specs.keys():
            if 'tag' in settings.linkpar_specs[label]:
                tag = settings.linkpar_specs[label]['tag']
                for pair in settings.links[tag]:
                    self.links.append(self.getComp(pair[0]).makeLinkTo(self.getComp(pair[1]),link_index=k,link_label=label))
                    if not tag in self.link_ids.keys():
                        self.link_ids[tag] = []
                    self.link_ids[tag].append(k)
                    k += 1
                    
        k = 0
        # Now create variable objects for dependencies, starting with characteristics.
        for label in settings.charac_specs.keys():
            if 'par_dependency' in settings.charac_specs[label]:
                self.deps.append(Variable(label=label))
                self.dep_ids[label] = k
                k += 1
        
        # Finish dependencies with untagged parameters, i.e. ones that are not considered transition flow rates.
        for label in settings.par_deps.keys():
            self.deps.append(Variable(label=label))
            self.dep_ids[label] = k
            k += 1
                    
            
#    def printCompVars(self, full = False):
#        ''' Loop through all compartments and print out current variable values. '''
#        for comp in self.comps:
#            if not full:
#                print('[Pop: %s][Compartment: %5s][Popsize: %15.4f]' % (self.label, comp.label, comp.popsize[self.t_index]))
#            else:
#                print('[Pop: %s][Compartment: %5s][Popsize...]' % (self.label, comp.label))
#                print(comp.popsize)

#    def printLinkVars(self, full = False):
#        ''' Loop through all links and print out current variable values. '''
#        for link in self.links:
#            if not full:
#                print('[Pop: %s][%5s --> %-5s][Transit. Frac.: %5.4f]' % (self.label, link.label_from, link.label_to, link.vals[self.t_index]))
#            else:
#                print('[Pop: %s][%5s --> %-5s][Transit. Fraction...]' % (self.label, link.label_from, link.label_to))
#                print(link.vals)
                
    def preAllocate(self, sim_settings):
        '''
        Pre-allocate variable arrays in compartments, links and dependent variables for faster processing.
        Array maintains initial value but pre-fills everything else with NaNs.
        Thus errors due to incorrect parset value saturation should be obvious from results.
        '''
        for comp in self.comps:
            init_popsize = comp.popsize[0]
            comp.popsize = np.ones(len(sim_settings['tvec']))*np.nan
            comp.popsize[0] = init_popsize
        for link in self.links:
            init_val = link.vals[0]
            link.vals = np.ones(len(sim_settings['tvec']))*np.nan
            link.vals[0] = init_val
        for dep in self.deps:
            init_val = dep.vals[0]
            dep.vals = np.ones(len(sim_settings['tvec']))*np.nan
            dep.vals[0] = init_val
            
            
            
#%% Model class
            
class Model(object):
    ''' A class to wrap up multiple populations within model and handle cross-population transitions. '''
    
    def __init__(self):
        
        self.pops = list()              # List of population groups that this model subdivides into.
        self.pop_ids = dict()           # Maps label of a population to its position index within populations list.
        
        self.sim_settings = odict()
        
        self.t_index = 0                # Keeps track of array index for current timepoint data within all compartments.
        
        
        
    def getPop(self, pop_label):
        ''' Allow model populations to be retrieved by label rather than index. '''
        pop_index = self.pop_ids[pop_label]
        return self.pops[pop_index]
        
        
    def build(self, settings, parset):
        ''' Build the full model. '''
        
        self.sim_settings['tvec'] = np.arange(settings.tvec_start, settings.tvec_end + settings.tvec_dt/2, settings.tvec_dt)
        
        
        for k, pop_label in enumerate(parset.pop_labels):
            self.pops.append(ModelPopulation(settings = settings, label = pop_label, index = k))
            self.pops[-1].preAllocate(self.sim_settings)     # Memory is allocated, speeding up model. However, values are NaN so as to enforce proper parset value saturation.
            self.pop_ids[pop_label] = k
                    
        # Propagating initial characteristic parset values into ModelPops.
        # NOTE: Extremely involved process, so might be worth extracting the next few paragraphs as a separate method.
        # First interpolate initial value for each relevant characteristic (i.e. one that has an entry point).
        # Maintaining definitional order (i.e. the order characteristics were defined in cascade workbook) is crucial.
        init_dict = odict()
        charac_for_entry = odict()
        t_init = np.array([self.sim_settings['tvec'][0]])
        for charac_label in settings.charac_specs.keys():
            if 'entry_point' in settings.charac_specs[charac_label].keys():
                entry_point = settings.charac_specs[charac_label]['entry_point']
                init_dict[charac_label] = odict()
                charac_for_entry[entry_point] = charac_label
                par = parset.pars['characs'][parset.par_ids['characs'][charac_label]]
                for pop_label in parset.pop_labels:
                    val = par.interpolate(tvec = t_init, pop_label = pop_label)
                    init_dict[charac_label][pop_label] = dcp(val)
        
        # Next, multiply out any denominators that exist. Again, definitional order matters.
        # These should all be other previously-defined entry-point characteristics, according to validation in settings.py.
        for charac_label in init_dict.keys():
            if 'denom' in settings.charac_specs[charac_label].keys():
                denom_label = settings.charac_specs[charac_label]['denom']
                entry_point = settings.charac_specs[charac_label]['entry_point']
                for pop_label in parset.pop_labels:
                    init_dict[charac_label][pop_label] *= init_dict[denom_label][pop_label]
#        print init_dict
#        print charac_for_entry
        
        # Then map each characteristic to those that include its entry point.
        sub_dict = odict()
        for charac_label in init_dict.keys():
            entry_point = settings.charac_specs[charac_label]['entry_point']
            flat_list, dep_list = flattenDict(input_dict = settings.charac_specs, base_key = charac_label, sub_keys = ['includes'])
            flat_list.remove(entry_point)
            for include in flat_list:
                if include in charac_for_entry.keys():
                    if charac_for_entry[include] not in sub_dict.keys(): sub_dict[charac_for_entry[include]] = []
                    sub_dict[charac_for_entry[include]].append(charac_label)
#        print sub_dict
                    
        # Flatten the previous mapping and remove duplicates to determine the smallest sets that include the entry point of each characteristic.
        # This process determines which values should be subtracted from which characteristics.
        # For example, if a is a subgroup of b and c, but b is a subgroup of c, subtracting both a and b from c would b double-counting.
        # In this example, b should be extracted from c and a should be extracted from b alone.
        minus_dict = odict()
        for charac_label in init_dict.keys():
            try:
                sub_list, key_list = flattenDict(input_dict = sub_dict, base_key = charac_label)
            except:
                sub_list = []
                key_list = [charac_label]
            combine_list = key_list + sub_list
            combine_list.remove(charac_label)
            minus_list = [k for k,v in Counter(combine_list).items() if v==1]
            minus_dict[charac_label] = minus_list
#        print minus_dict
        
        # Another loop to set compartment values.
        for charac_label in init_dict.keys():
            entry_point = settings.charac_specs[charac_label]['entry_point']
            for pop_label in parset.pop_labels:
                val = init_dict[charac_label][pop_label]
                self.getPop(pop_label).getComp(entry_point).popsize[0] = val
        
        # A final loop to apply the subtractions determined previously (e.g. susceptibles = total - infected).        
        for minus_label in minus_dict.keys():
            for charac_label in minus_dict[minus_label]:
                entry_point = settings.charac_specs[charac_label]['entry_point']
                for pop_label in parset.pop_labels:
                    val = init_dict[minus_label][pop_label]
                    self.getPop(pop_label).getComp(entry_point).popsize[0] -= val
                    
#        for charac_label in init_dict.keys():
#            entry_point = settings.charac_specs[charac_label]['entry_point']
#            print entry_point
#            print self.getPop(pop_label).getComp(entry_point).popsize[0]
                    
        
        # Propagating cascade parameter parset values into ModelPops. Handle both 'tagged' links and 'untagged' dependencies.
        for par in parset.pars['cascade']:
            
            if 'tag' in settings.linkpar_specs[par.label]:
                tag = settings.linkpar_specs[par.label]['tag']                  # Map parameter label to link tag.
                for pop_label in parset.pop_labels:
                    for link_id in self.getPop(pop_label).link_ids[tag]:        # Map link tag to link id in ModelPop.
                        self.getPop(pop_label).links[link_id].vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_label)
                        self.getPop(pop_label).links[link_id].val_format = par.y_format[pop_label]
                        self.getPop(pop_label).links[link_id].scale_factor = par.y_factor[pop_label]
            else:
                for pop_label in parset.pop_labels:
                    dep_id = self.getPop(pop_label).dep_ids[par.label]          # Map dependency label to dependency id in ModelPop.
                    self.getPop(pop_label).deps[dep_id].vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_label)
                    self.getPop(pop_label).deps[dep_id].val_format = par.y_format[pop_label]
                    self.getPop(pop_label).deps[dep_id].scale_factor = par.y_factor[pop_label]

        
        # Propagating transfer parameter parset values into Model object.
        for trans_type in parset.transfers.keys():
            if parset.transfers[trans_type]:
                for pop_source in parset.transfers[trans_type].keys():
                    par = parset.transfers[trans_type][pop_source]
                    
                    for pop_target in par.y:
                        for comp in self.getPop(pop_source).comps:
                            trans_tag = comp.label + '_' + trans_type + '_to_' + pop_target       # NOTE: Perhaps there is a nicer way to set up transfer tagging.
                            if not comp.tag_birth and not comp.tag_dead and not comp.junction:
                                
                                num_links = len(self.getPop(pop_source).links)
                                link = comp.makeLinkTo(self.getPop(pop_target).getComp(comp.label),link_index=num_links,link_label=trans_tag)
                                link.vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_target)
                                link.val_format = par.y_format[pop_target]
                                link.scale_factor = par.y_factor[pop_target]
                                
                                self.getPop(pop_source).links.append(link)
                                self.getPop(pop_source).link_ids[trans_tag] = [num_links]
        
        # Make sure initially-filled junctions are processed and initial dependencies are calculated.
        self.processJunctions(settings = settings)
        self.updateDependencies(settings = settings)


        # set up sim_settings for later use wrt population tags
        self.sim_settings['tag_birth'] = []
        self.sim_settings['tag_death'] = []
        self.sim_settings['tag_no_plot'] = []
        for node_label in settings.node_specs.keys():
            for tag in ['tag_birth','tag_death','tag_no_plot']:
                if tag in settings.node_specs[node_label]:
                    self.sim_settings[tag] = node_label

            
    def process(self, settings):
        ''' 
        Run the full model. 
        
        Note that updateBirths is called before stepForward: this can be set to be after the other interpopulation transitions, but requires
        some care with how the values are set up. 
        '''
        
        for t in self.sim_settings['tvec'][1:]:
#            self.printModelState(self.t_index)
            self.stepForward(settings = settings, dt = settings.tvec_dt)
            self.processJunctions(settings = settings)
            self.updateDependencies(settings = settings)
        
        return self.pops, self.sim_settings

    def _calculateDPops(self,settings,ti,dt,modified_compartment={},empty_compartment=[],reset_compartment=[]):
        """
        
        This function gets called from model.stepForward() but also from within a validation check. 
        
        Params:
            settings
            ti
            dt
            modified_compartment optional: Dictionary of (did,reduction_ratio) pairs that indicate by what percentage an dpop amount should be reduced.  
            empty_compartment    optional: List of compartment did that are empty
            reset_compartment    optional: List of compartment did that already have a negative value
        
        """
        
        # Preallocate a change-in-popsize array. Requires each population to have the same cascade.
        num_pops = len(self.pops)
        num_comps = len(self.pops[0].comps)         # NOTE: Hard-coded reference to 'zeroth' population. Improve later.
        dpopsize = np.zeros(num_pops*num_comps)
        # Variables for tracking in and out amounts to a compartment, separately
        #dpop_in = np.zeros(num_pops*num_comps)
        dpop_out = np.zeros(num_pops*num_comps)
        
        for pop in self.pops:
            for comp in pop.comps:
                
                # If not pre-allocated, extend compartment variable arrays. Either way copy current value to the next position in the array.
                if not len(comp.popsize) > ti + 1:
                    comp.popsize = np.append(comp.popsize, 0.0)
                if not len(comp.popsize) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                    raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment "%s".' % (comp.label))
                comp.popsize[ti+1] = comp.popsize[ti]
                    
                if comp.label not in settings.junction_labels:      # Junctions collect inflows during this step. They do not process outflows here.
                    for lid in comp.outlink_ids:
                        link = pop.links[lid]
                        
                        # If link values are not pre-allocated to the required time-vector length, extend with the same value as the last index.
                        if not len(link.vals) > ti + 1:
                            link.vals = np.append(link.vals, link.vals[-1])
                        if not len(link.vals) > ti + 1:         # If one extension did not create an index of ti+1, something is seriously wrong...
                            raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment "%s".' % (link.label))
                        
                        did_from = link.index_from[0] * num_comps + link.index_from[1]
                        did_to = link.index_to[0] * num_comps + link.index_to[1]
                        comp_source = self.pops[link.index_from[0]].getComp(link.label_from)
                        
#                        # Evaluate custom function for variable if it exists and overwrite any value that is currently stored for transition.
#                        if link.label in settings.linkpar_specs:
#                            if 'f_stack' in settings.linkpar_specs[link.label]:
#                                f_stack = dcp(settings.linkpar_specs[link.label]['f_stack'])
#                                deps = dcp(settings.linkpar_specs[link.label]['deps'])
#                                for dep_label in deps.keys():
#                                    deps[dep_label] = pop.getDep(dep_label).vals[ti]
#                                link.vals[ti] = settings.parser.evaluateStack(stack = f_stack, deps = deps)
                        
                        
                        convert_amt = 0 
                        
                        transition = link.vals[ti]
                        
                        if link.scale_factor is not None and link.scale_factor != project_settings.DO_NOT_SCALE : # scale factor should be available to be used 
                            transition *= link.scale_factor
                        
                        
                        if link.val_format == 'fraction': 
                            converted_frac = 1 - (1 - transition) ** dt      # A formula for converting from yearly fraction values to the dt equivalent.
                            converted_amt = comp_source.popsize[ti] * converted_frac
                        elif link.val_format == 'number':
                            converted_amt = transition * dt
                        else:
                            raise OptimaException("Unknown link type: %s in model\nObserved for population %s, compartment %s"%(link.val_format,pop.label,comp.label))
                        
                        
                        # If we've observed that this link contributes to the creation of a negative population size, we modify the amounts: 
                        if did_from in empty_compartment:
                            # If the compartment is an empty one, we shouldn't update
                            logging.debug("Empty compartment: no change made for link from=%g,to=%g"%(did_from,did_to))
                            converted_amt = 0
                        
                        elif did_from in modified_compartment.keys():
                            # during validation, we encountered an outflow from a compartment that was more than the population size
                            logging.debug("Modifying flow amount for link (from=%g,to=%g) by %g: prev amount = %g, new amount = %g"%(did_from,did_to,modified_compartment[did_from],converted_amt,modified_compartment[did_from]*converted_amt))
                            converted_amt *= modified_compartment[did_from]
                        
                        elif did_from in reset_compartment:
                            # set converted amount in the other direction. Hmm still not perfect
                            converted_amt = 0
                            
                            
                        dpopsize[did_from] -= converted_amt
                        dpopsize[did_to] += converted_amt
                        
                        # Tracking in and out for a compartment
                        #dpop_in[did_to] += converted_amt
                        dpop_out[did_from] += converted_amt
                            
        return dpopsize, dpop_out
    
    

    def stepForward(self, settings, dt = 1.0):
        '''
        Evolve model characteristics by one timestep (defaulting as 1 year).
        Each application of this method writes calculated values to the next position in popsize arrays, regardless of dt.
        Thus the corresponding time vector associated with variable dt steps must be tracked externally.
        '''
        
        ti = self.t_index
        
        # Requires each population to have the same cascade.
        num_pops = len(self.pops)
        num_comps = len(self.pops[0].comps)         # NOTE: Hard-coded reference to 'zeroth' population. Improve later.
        
        
        # First loop through all pops and compartments to calculate value changes.
        dpopsize, dpop_out = self._calculateDPops(settings,ti,dt)
                  
        # Check calculated rates for proposed dpop value. Note that this should be made an iterative process for 
        # when we avert an incorrect dpop change, as the avert may result in new errors. 
        # TODO: create validation as an iterative process (currently found that it was possible that there was no good solution found
        """
        # SUGGESTION: 
        count = 0 
        while count < some_limit_count: 
            if isDPopValid(self,settings,dpopsize,ti,validation):
                break
            dpopsize,dpop_out,reset_ids = checkNegativePopulation(self,settings,dpopsize,dpop_out,ti,validation)
        """
        validation = settings.validation['negative_population']
        dpopsize,dpop_out,reset_ids = checkNegativePopulation(self,settings,dpopsize,dpop_out,ti,dt,validation)

        # Second loop through all pops and comps to apply value changes at next timestep index.
        for pid in xrange(num_pops):
            for cid in xrange(num_comps):
                did = pid * num_comps + cid
                self.pops[pid].comps[cid].popsize[ti+1] += dpopsize[did]
        
        # Update timestep index.
        self.t_index += 1
        for pop in self.pops:
            pop.t_index = self.t_index       # Keeps ModelPopulations synchronised with Model object.


    def processJunctions(self, settings):
        '''
        For every compartment considered a junction, propagate the contents onwards until all junctions are empty.
        '''
        
        ti = self.t_index
        ti_link = ti - 1
        if ti_link < 0: ti_link = ti    # For the case where junctions are processed immediately after model initialisation.
        final_review = False
        
        review_count = 0
        while not final_review:
            if review_count > settings.recursion_limit: raise OptimaException('ERROR: Processing junctions (i.e. propagating contents onwards) for timestep %i is taking far too long. Infinite loop suspected.' % ti_link)
            final_review = True     # Assume that this is the final sweep through junctions to verify their emptiness.
            for pop in self.pops:
                for junction_label in settings.junction_labels:
                    comp = pop.getComp(junction_label)
                    popsize = comp.popsize[ti]
                    if popsize < 0:     # NOTE: Hacky fix for negative values. Needs work in case of super-negatives.
                        comp.popsize[ti] = 0
                        popsize = 0
                    if popsize > project_settings.TOLERANCE:  # NOTE: Hard-coded tolerance. Bad.
                        final_review = False    # Outflows could propagate into other junctions requiring another review.
                        denom_val = sum(pop.links[lid].vals[ti_link] for lid in comp.outlink_ids)
                        if denom_val == 0: raise OptimaException('ERROR: Proportions for junction "%s" outflows sum to zero, resulting in a nonsensical ratio. There may even be (invalidly) no outgoing transitions for this junction.' % junction_label)
                        for lid in comp.outlink_ids:
                            link = pop.links[lid]
                            
                            comp.popsize[ti] -= popsize * link.vals[ti_link] / denom_val
                            pop.getComp(link.label_to).popsize[ti] += popsize * link.vals[ti_link] / denom_val
            review_count += 1



    def updateDependencies(self, settings):
        '''
        Run through all parameters and characteristics flagged as dependencies for custom-function parameters and evaluate them for the current timestep.
        These dependencies must be calculated in the same order as defined in settings, characteristics before parameters, otherwise references may break.
        Also, parameters in the dependency list do not need to be calculated unless explicitly depending on another parameter.
        '''
        
        ti = self.t_index
        
        for pop in self.pops:
            for dep in pop.deps:
                if dep.label in settings.charac_deps.keys():
                    dep.vals[ti] = 0                                    
                    
                    # Sum up all relevant compartment popsizes (or previously calculated characteristics).
                    for inc_label in settings.charac_specs[dep.label]['includes']:
                        if inc_label in pop.comp_ids.keys():
                            val = pop.getComp(inc_label).popsize[ti]
                        elif inc_label in pop.dep_ids.keys():    # NOTE: This should not select a parameter-type dependency due to settings validation, but can validate here if desired.
                            val = pop.getDep(inc_label).vals[ti]
                        else:
                            raise OptimaException('ERROR: Compartment or characteristic "%s" has not been pre-calculated for use in calculating "%s".' % (inc_label, dep.label))
                        
                        dep.vals[ti] += val
                    
                    # Divide by relevant compartment popsize (or previously calculated characteristic).
                    if 'denom' in settings.charac_specs[dep.label]:
                        den_label = settings.charac_specs[dep.label]['denom']
                        if den_label in pop.dep_ids.keys():  # NOTE: See above note for avoiding parameter-type dependencies.
                            val = pop.getDep(den_label).vals[ti]
                        elif den_label in pop.comp_ids.keys():
                            val = pop.getComp(den_label).popsize[ti]
                        else:
                            raise OptimaException('ERROR: Compartment or characteristic "%s" has not been pre-calculated for use in calculating "%s".' % (inc_label, dep.label))                      
                        
                        dep.vals[ti] /= val
                
#                # If the dependency is a parameter, evaluate its stack.
#                elif dep.label in settings.linkpar_specs.keys():
#                    if 'deps' in settings.linkpar_specs[dep.label]:
#                        if len(settings.linkpar_specs[dep.label]['deps'].keys()) >= 0:
#                            f_stack = dcp(settings.linkpar_specs[dep.label]['f_stack'])
#                            dep_deps = dcp(settings.linkpar_specs[dep.label]['deps'])
#                            for dep_dep_label in dep_deps.keys():
#                                dep_deps[dep_dep_label] = pop.getDep(dep_dep_label).vals[ti]
#                            dep.vals[ti] = settings.parser.evaluateStack(stack = f_stack, deps = dep_deps)

#                else:
#                    raise OptimaException('ERROR: Dependency "%s" does not appear to be either a characteristic or parameter.' % (dep.label))
                    
            for par_label in settings.par_funcs.keys():
                pars = []
                if par_label in settings.par_deps:
                    pars.append(pop.getDep(par_label))
                else:
                    pars = pop.getLinks(settings.linkpar_specs[par_label]['tag'])
                specs = settings.linkpar_specs[par_label]
                f_stack = dcp(specs['f_stack'])
                deps = dcp(specs['deps'])
                for dep_label in deps.keys():
                    if dep_label in settings.par_deps.keys() or dep_label in settings.charac_deps.keys():
                        val = pop.getDep(dep_label).vals[ti]
                    else:
                        val = pop.getLinks(settings.linkpar_specs[dep_label]['tag'])[0].vals[ti]    # As links are duplicated for the same tag, can pull values from the zeroth one.
                    deps[dep_label] = val
                new_val = settings.parser.evaluateStack(stack = f_stack, deps = deps)
                for par in pars:
                    par.vals[ti] = new_val
                

    
    def calculateOutputs(self, settings):
        '''
        Calculate outputs (called cascade characteristics in settings).
        These outputs must be calculated in the same order as defined in settings, otherwise references may break.
        '''
        
        outputs = odict()
        for cid in settings.charac_specs.keys():
            outputs[cid] = odict()
            for pop in self.pops:
                outputs[cid][pop.label] = None
                
                # Sum up all relevant compartment popsizes (or previously calculated characteristics).
                for inc_label in settings.charac_specs[cid]['includes']:
                    if inc_label in self.getPop(pop.label).comp_ids.keys():
                        vals = self.getPop(pop.label).getComp(inc_label).popsize
                    elif inc_label in outputs.keys()[:-1]:
                        vals = outputs[inc_label][pop.label]
                    else:
                        raise OptimaException('ERROR: Compartment or characteristic "%s" has not been pre-calculated for use in calculating "%s".' % (inc_label, cid))
                        
                    if outputs[cid][pop.label] is None:
                        outputs[cid][pop.label] = dcp(vals)
                    else:
                        outputs[cid][pop.label] += vals
                
                # Divide by relevant compartment popsize (or previously calculated characteristic).
                if 'denom' in settings.charac_specs[cid]:
                    den_label = settings.charac_specs[cid]['denom']
                    if den_label in outputs.keys()[:-1]:
                        vals = outputs[den_label][pop.label]
                    elif den_label in self.getPop(pop.label).comp_ids.keys():
                        vals = self.getPop(pop.label).getComp(den_label).popsize
                    else:
                        raise OptimaException('ERROR: Compartment or characteristic "%s" has not been pre-calculated for use in calculating "%s".' % (inc_label, cid))
                    
                    outputs[cid][pop.label] /= (vals+project_settings.TOLERANCE)
                    
        return outputs
    
    def printModelState(self,ti):
        for pop in self.pops:
            logging.info("Population: %s"%pop.label)
            logging.info(pop.getModelState(ti))
            
    
    def isBirthCompartment(self,comp_id,settings):
        """
        Returns True if comp_label is label for a compartment tagged as births. 
        """
        if 'tag_birth' in settings.node_specs[comp_id].keys():
            return True
        return False
    
    def isDeadCompartment(self,comp_id,settings):
        """
        Returns True if comp_label is label for a compartment tagged as dead. 
    
        """
        if 'tag_dead' in settings.node_specs[comp_id].keys():
            return True
        return False
        
    

#%% Model function (simulates epidemic dynamics)

def runModel(settings, parset):
    ''' Processes the TB epidemiological model. '''

    m = Model()
    m.build(settings = settings, parset = parset)
    m.process(settings = settings)
    
    # Intention: to replace these 3 lines
    #outputs = m.calculateOutputs(settings = settings)
    #m_pops = m.pops
    #sim_settings = m.sim_settings
    # with 
    results = ResultSet(m,parset,settings)
    # and
    #return results
    
    # For moment, we'll return all, as further plotting methods require  m_pops, sim_settings, outputs
    # (which are within results, but we'll do it in a piecewise manner at the moment)
    #%% Collect and return raw results    
    return results#, m_pops, sim_settings, outputs


