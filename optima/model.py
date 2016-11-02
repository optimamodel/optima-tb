#%% Imports

from utils import flattenDict, odict, OptimaException

import numpy as np
from copy import deepcopy as dcp



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
        self.val_format = 'fraction'

class Link(Variable):
    '''
    A more involved version of Variable, representing unidirectional flow between two objects in a network.
    The values stored in this extended version of Variable refer to flow rates.
    If used in ModelPop, the Link refers to two cascade compartments within a single population.
    If used in Model, the Link refers to two distinct population groups.
    In the latter case, intended logic should transfer agents between all (non-dead) corresponding compartments.
    '''
    def __init__(self, object_from, object_to, label = 'default', val = 0.0):
        Variable.__init__(self, label = label, val = val)
        self.index_from = object_from.index
        self.index_to = object_to.index
        self.label_from = object_from.label
        self.label_to = object_to.label



#%% Cascade compartment and population classes

class ModelCompartment(Node):
    ''' A class to wrap up data for one compartment within a cascade network. '''
    
    def __init__(self, label = 'default', index = 0, popsize = 0.0):
        Node.__init__(self, label = label, index = index)
        self.popsize = np.array([float(popsize)])   # Number of people in compartment.
        self.tag_dead = False                       # Tag for whether this compartment contains dead people.
        self.junction = False
    

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
    
    def getComp(self, comp_label):
        ''' Allow compartments to be retrieved by label rather than index. '''
        comp_index = self.comp_ids[comp_label]
        return self.comps[comp_index]
        
    def getDep(self, dep_label):
        ''' Allow dependencies to be retrieved by label rather than index. '''
        dep_index = self.dep_ids[dep_label]
        return self.deps[dep_index]
        
    def genCascade(self, settings):
        '''
        Generate a compartmental cascade as defined in a settings object.
        Fill out the compartment and transition lists within the model population object.
        '''
        for k, label in enumerate(settings.node_specs.keys()):
            self.comps.append(ModelCompartment(label = label, index = (self.index,k)))
            if 'tag_dead' in settings.node_specs[label].keys():
                self.comps[-1].tag_dead = True
            if 'junction' in settings.node_specs[label].keys():
                self.comps[-1].junction = True
            self.comp_ids[label] = k
        j = 0
        k = 0
        for label in settings.linkpar_specs.keys():
            if 'tag' in settings.linkpar_specs[label]:
                tag = settings.linkpar_specs[label]['tag']
                for pair in settings.links[tag]:
                    self.links.append(self.getComp(pair[0]).makeLinkTo(self.getComp(pair[1]),link_index=j,link_label=label))
                    if not tag in self.link_ids.keys():
                        self.link_ids[tag] = []
                    self.link_ids[tag].append(j)
                    j += 1
            else:
                self.deps.append(Variable(label=label))
                self.dep_ids[label] = k
                k += 1
        for label in settings.charac_specs.keys():
            if 'par_dependency' in settings.charac_specs[label]:
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
        
        # One more loop to subtract out any included characteristics from each characteristic (e.g. susceptibles = total - infected).
        for charac_label in init_dict.keys():
            entry_point = settings.charac_specs[charac_label]['entry_point']
            for pop_label in parset.pop_labels:
                val = init_dict[charac_label][pop_label]
                flat_list, dep_list = flattenDict(input_dict = settings.charac_specs, base_key = charac_label, sub_keys = ['includes'])
                flat_list.remove(entry_point)
                for include in flat_list:
                    if include in charac_for_entry.keys():
                        val -= init_dict[charac_for_entry[include]][pop_label]
                self.getPop(pop_label).getComp(entry_point).popsize[0] = val
                
                # While looping through characteristics/populations, fill deps list with those that must be calculated per timestep.
                if 'par_dependency' in settings.charac_specs[charac_label]:
                    pass   #NECESSARY?
                    
        
        # Propagating cascade parameter parset values into ModelPops. Handle both 'tagged' links and 'untagged' dependencies.
        for par in parset.pars['cascade']:
            if 'tag' in settings.linkpar_specs[par.label]:
                tag = settings.linkpar_specs[par.label]['tag']                  # Map parameter label to link tag.
                for pop_label in parset.pop_labels:
                    for link_id in self.getPop(pop_label).link_ids[tag]:        # Map link tag to link id in ModelPop.            
                        self.getPop(pop_label).links[link_id].vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_label)
                        self.getPop(pop_label).links[link_id].val_format = par.y_format[pop_label]
            else:
                for pop_label in parset.pop_labels:
                    dep_id = self.getPop(pop_label).dep_ids[par.label]          # Map dependency label to dependency id in ModelPop.
                    self.getPop(pop_label).deps[dep_id].vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_label)
                    self.getPop(pop_label).deps[dep_id].val_format = par.y_format[pop_label]
                
        
        # Propagating transfer parameter parset values into Model object.
        for trans_type in parset.transfers.keys():
            if parset.transfers[trans_type]:
                for pop_source in parset.transfers[trans_type].keys():
                    par = parset.transfers[trans_type][pop_source]
                    for pop_target in par.y:
                        for comp in self.getPop(pop_source).comps:
                            trans_tag = comp.label + '_' + trans_type + '_to_' + pop_target       # NOTE: Perhaps there is a nicer way to set up transfer tagging.
                            if not comp.tag_dead and not comp.junction:
                                num_links = len(self.getPop(pop_source).links)
                                link = comp.makeLinkTo(self.getPop(pop_target).getComp(comp.label),link_index=num_links,link_label=trans_tag)
                                link.vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_target)
                                link.val_format = par.y_format[pop_target]
                                
                                self.getPop(pop_source).links.append(link)
                                self.getPop(pop_source).link_ids[trans_tag] = [num_links]
                                
        self.processJunctions(settings = settings)      # Junctions may be initialised with popsize by this stage. Flow that onwards.
                        
                
    def process(self, settings):
        ''' Run the full model. '''
        
        for t in self.sim_settings['tvec'][1:]:
            self.stepForward(settings = settings, dt = settings.tvec_dt)
            self.processJunctions(settings = settings)
        
        return self.pops, self.sim_settings
        
    
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
                    
                    outputs[cid][pop.label] /= vals
        
        return outputs
        
    def processJunctions(self, settings):
        '''
        For every compartment considered a junction, propagate the contents onwards until all junctions are empty.
        '''
        
        ti = self.t_index
        ti_link = ti - 1
        if ti_link < 0: ti_link = ti    # For the case where junctions are processed immediately after model initialisation.
        final_review = False
        
        while not final_review:
            final_review = True     # Assume that this is the final sweep through junctions to verify their emptiness.
            for pop in self.pops:
                for junction_label in settings.junction_labels:
                    comp = pop.getComp(junction_label)
                    popsize = comp.popsize[ti]
                    if popsize > 0:
                        final_review = False    # Outflows could propagate into other junctions requiring another review.
                        denom_val = sum(pop.links[lid].vals[ti_link] for lid in comp.outlink_ids)
                        if denom_val == 0: raise OptimaException('ERROR: Proportions for junction "%s" outflows sum to zero, resulting in a nonsensical ratio. There may even be (invalidly) no outgoing transitions for this junction.' % junction_label)
                        for lid in comp.outlink_ids:
                            link = pop.links[lid]
                            
                            comp.popsize[ti] -= popsize * link.vals[ti_link] / denom_val
                            pop.getComp(link.label_to).popsize[ti] += popsize * link.vals[ti_link] / denom_val
        

    def stepForward(self, settings, dt = 1.0):
        '''
        Evolve model characteristics by one timestep (defaulting as 1 year).
        Each application of this method writes calculated values to the next position in popsize arrays, regardless of dt.
        Thus the corresponding time vector associated with variable dt steps must be tracked externally.
        '''
        
        ti = self.t_index
        
        # Preallocate a change-in-popsize array. Requires each population to have the same cascade.
        num_pops = len(self.pops)
        num_comps = len(self.pops[0].comps)         # NOTE: Hard-coded reference to 'zeroth' population. Improve later.
        dpopsize = np.zeros(num_pops*num_comps)
        
        # First loop through all pops and comps to calculate value changes.
        k = 0
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
                        
                        # Evaluate custom function for variable if it exists and overwrite any value that is currently stored for transition.
                        if link.label in settings.linkpar_specs:
                            if 'f_stack' in settings.linkpar_specs[link.label]:
                                f_stack = dcp(settings.linkpar_specs[link.label]['f_stack'])
                                deps = dcp(settings.linkpar_specs[link.label]['deps'])
                                for dep_label in deps.keys():
                                    deps[dep_label] = pop.getDep(dep_label).vals[ti]
                                link.vals[ti] = settings.parser.evaluateStack(stack = f_stack, deps = deps)
                        
                        converted_frac = 1 - (1 - link.vals[ti]) ** dt      # A formula for converting from yearly fraction values to the dt equivalent.
                        
                        dpopsize[did_from] -= comp_source.popsize[ti] * converted_frac
                        dpopsize[did_to] += comp_source.popsize[ti] * converted_frac
                
                k += 1

        # Second loop through all pops and comps to apply value changes at next timestep index.
        for pid in xrange(num_pops):
            for cid in xrange(num_comps):
                did = pid * num_comps + cid
                self.pops[pid].comps[cid].popsize[ti+1] += dpopsize[did]
        
        # Update timestep index.
        self.t_index += 1
        for pop in self.pops:
            pop.t_index = self.t_index       # Keeps ModelPopulations synchronised with Model object.
        
    

#%% Model function (simulates epidemic dynamics)

def runModel(settings, parset):
    ''' Processes the TB epidemiological model. '''
    
    m = Model()
    m.build(settings = settings, parset = parset)
    m.process(settings = settings)
    outputs = m.calculateOutputs(settings = settings)
    m_pops = m.pops
    sim_settings = m.sim_settings
    
    
    
    #%% Collect and return raw results    
    
    return m_pops, sim_settings, outputs