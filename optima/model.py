#%% Imports

from utils import odict, OptimaException

import numpy as np
import numpy.random as npr
from copy import deepcopy as dcp



#%% Abstract classes used in model

class Node(object):
    ''' Lightweight abstract class to represent one object within a network. '''
    def __init__(self, label = 'default', index = 0):
        self.label = label          # Reference name for this object.
        self.index = index          # Positional index for container list.
        self.num_outlinks = 0       # Tracks number of nodes this one is linked to as an initial node.
        self.outlink_ids = []       # List of indices corresponding to each outgoing link.
        
    def makeLinkTo(self, other_node, link_index):
        '''
        Link this node to another (i.e. create a transition link).
        Must provide an index that, by intention, represents where this link is positioned in its storage container. 
        '''
        if not isinstance(other_node, Node):
            raise OptimaException('ERROR: Attempting to link compartment to something that is not a compartment.')
        self.num_outlinks += 1
        self.outlink_ids.append(link_index)
        return Link(self, other_node)

class Link(object):
    '''
    Lightweight abstract class to represent unidirectional flow between two objects in a network.
    If used in ModelPop, the Link refers to two cascade compartments within a single population.
    If used in Model, the Link refers to two distinct population groups.
    In the latter case, intended logic should transfer agents between all (non-dead) corresponding compartments.
    '''
    def __init__(self, object_from, object_to, val = 0.0):
        self.index_from = object_from.index
        self.index_to = object_to.index
        self.label_from = object_from.label
        self.label_to = object_to.label
        self.vals = np.array([float(val)])   # An abstract array of values. Usually relating to flow rate.
        self.val_format = 'probability'



#%% Cascade compartment and population classes

class ModelCompartment(Node):
    ''' A class to wrap up data for one compartment within a cascade network. '''
    
    def __init__(self, label = 'default', index = 0, popsize = 0.0):
        Node.__init__(self, label = label, index = index)
        self.popsize = np.array([float(popsize)])   # Number of people in compartment.
        self.tag_dead = False                       # Tag for whether this compartment contains dead people.
    

class ModelPopulation(Node): 
    '''
    A class to wrap up data for one population within model.
    Each model population must contain a set of compartments with equivalent labels.
    '''
    
    def __init__(self, settings, label = 'default', index = 0):
        Node.__init__(self, label = label, index = index)
        self.comps = list()         # List of cascade compartments that this model population subdivides into.
        self.links = list()         # List of intra-population cascade transitions within this model population.
        self.comp_ids = dict()      # Maps label of a compartment to its position index within compartments list.
        self.link_ids = dict()      # Maps cascade transition tag to indices for all relevant transitions within links list.
        self.t_index = 0            # Keeps track of array index for current timepoint data within all compartments.
        
        self.genCascade(settings = settings)    # Convert compartmental cascade into lists of compartment and link objects.
    
    def getComp(self, comp_label):
        ''' Allow compartments to be retrieved by label rather than index. '''
        comp_index = self.comp_ids[comp_label]
        return self.comps[comp_index]
        
    def genCascade(self, settings):
        '''
        Generate a compartmental cascade as defined in a settings object.
        Fill out the compartment and transition lists within the model population object.
        '''
        for k, label in enumerate(settings.node_specs.keys()):
            self.comps.append(ModelCompartment(label = label, index = k))
            if 'tag_dead' in settings.node_specs[label].keys():
                self.comps[-1].tag_dead = True
            self.comp_ids[label] = k
        k = 0
        for tag in settings.links.keys():
            for pair in settings.links[tag]:
                self.links.append(self.getComp(pair[0]).makeLinkTo(self.getComp(pair[1]),link_index=k))
                if not tag in self.link_ids:
                    self.link_ids[tag] = []
                self.link_ids[tag].append(k)
                k += 1
    
    def stepCascadeForward(self, dt = 1.0):
        '''
        Evolve model population characteristics by one timestep (defaulting as 1 year).
        Calculated outputs like popsize will overwrite pre-allocated NaNs.
        In contrast, pre-allocated arrays that are used in calculations (e.g. parameters) must be pre-filled with appropriate values.
        '''
        
        ti = self.t_index
        
        # If not pre-allocated, extend compartment variable arrays. Either way copy current value to the next position in the array.
        for comp in self.comps:
            if not len(comp.popsize) > ti + 1:
                comp.popsize = np.append(comp.popsize, 0.0)
            if not len(comp.popsize) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment %s.' % (comp.label))
            comp.popsize[ti+1] = comp.popsize[ti]
                
        dpopsize = np.zeros(len(self.links))        
        
        # First calculation loop. Extend link variable arrays if not pre-allocated. Calculate value changes for next timestep.
        for k, link in enumerate(self.links):
            if not len(link.vals) > ti + 1:
                link.vals = np.append(link.vals, link.vals[-1])
            if not len(link.vals) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment %s.' % (link.label))
            
            converted_frac = 1 - (1 - link.vals[ti]) ** dt      # A formula for converting from yearly fraction values to the dt equivalent.
            dpopsize[k] = self.getComp(link.label_from).popsize[ti] * converted_frac

        # Second calculation loop. Apply value changes at next timestep.
        for k, link in enumerate(self.links):
            self.getComp(link.label_from).popsize[ti+1] -= dpopsize[k]
            self.getComp(link.label_to).popsize[ti+1] += dpopsize[k]
            
        self.t_index += 1       # Update timestep index.
            
    def printCompVars(self, full = False):
        ''' Loop through all compartments and print out current variable values. '''
        for comp in self.comps:
            if not full:
                print('[Pop: %s][Compartment: %5s][Popsize: %15.4f]' % (self.label, comp.label, comp.popsize[self.t_index]))
            else:
                print('[Pop: %s][Compartment: %5s][Popsize...]' % (self.label, comp.label))
                print(comp.popsize)

    def printLinkVars(self, full = False):
        ''' Loop through all links and print out current variable values. '''
        for link in self.links:
            if not full:
                print('[Pop: %s][%5s --> %-5s][Transit. Frac.: %5.4f]' % (self.label, link.label_from, link.label_to, link.vals[self.t_index]))
            else:
                print('[Pop: %s][%5s --> %-5s][Transit. Fraction...]' % (self.label, link.label_from, link.label_to))
                print(link.vals)

    def makeRandomVars(self, for_comp = True, for_link = True):
        ''' Randomise all compartment and link variables. Method used primarily for debugging. '''
        if for_comp:
            for comp in self.comps:
                comp.popsize[self.t_index] = npr.rand()*1e7
        if for_link:
            for link in self.links:
                link.vals = npr.rand()/self.getComp(link.label_from).num_outlinks     # Scaling makes sure fractions leaving a compartment sum to less than 1.
                
    def preAllocate(self, sim_settings):
        '''
        Pre-allocate variable arrays in compartments and links for faster processing.
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
            
            
            
#%% Model class
            
class Model(object):
    ''' A class to wrap up multiple populations within model and handle cross-population transitions. '''
    
    def __init__(self):
        
        self.pops = list()              # List of population groups that this model subdivides into.
        self.transfers = list()         # List of inter-population transitions (i.e. bulk transfers) within this model.
        self.pop_ids = dict()           # Maps label of a population to its position index within populations list.
        self.transfer_ids = dict()      # Maps label of an inter-population transfer to position index within transfers list.      
        
#        self.pops = odict()
        self.sim_settings = odict()
        
#        self.transfers = odict()
        
        
        
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
            
            self.pops[-1].getComp('sus').popsize[0] = 1000000   # NOTE: Temporary. Initial values inserted here.
            
        # Propagating cascade parameter parset values into ModelPops.
        for par in parset.pars:
            tag = settings.linkpar_specs[par.label]['tag']          # Map parameter label -> link tag.
            for pop_label in parset.pop_labels:
                for link_id in self.getPop(pop_label).link_ids[tag]:           # Map link tag -> link id in ModelPop.            
                    self.getPop(pop_label).links[link_id].vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_label)
                    self.getPop(pop_label).links[link_id].val_format = par.y_format[pop_label]
        
        # Propagating transfer parameter parset values into Model object.
        k = 0
        for trans_type in parset.transfers.keys():
            if parset.transfers[trans_type]:
                for pop_source in parset.transfers[trans_type].keys():
                    par = parset.transfers[trans_type][pop_source]
                    for pop_sink in par.y:
                        trans_tag = trans_type + '_' + pop_source + '_' + pop_sink       # NOTE: Perhaps there is a nicer way to set up transfer tagging.
                        self.transfers.append(self.getPop(pop_source).makeLinkTo(self.getPop(pop_sink),link_index=k))
                        self.transfers[-1].vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_sink)
                        self.transfers[-1].val_format = par.y_format[pop_sink]
                        self.transfer_ids[trans_tag] = k
                        k += 1
                        
                
    def process(self, settings):
        ''' Run the full model. '''
        
        for t in self.sim_settings['tvec'][1:]:
            for pop in self.pops:
                pop.stepCascadeForward(dt = settings.tvec_dt)
                
            # Apply transfers.
            for transfer in self.transfers:
                pop_source = transfer.label_from
                pop_sink = transfer.label_to
                tid_source = self.getPop(pop_source).t_index
                tid_sink = self.getPop(pop_sink).t_index
                for comp_label in self.getPop(pop_source).comp_ids.keys():
                    if not self.getPop(pop_source).getComp(comp_label).tag_dead:
                        
                        # Use 't-1' value of transfer rate on forward-stepped 't' value of popsize.
                        converted_frac = 1 - (1 - transfer.vals[tid_source-1]) ** settings.tvec_dt      # A formula for converting from yearly fraction values to the dt equivalent.
                        num_trans = self.getPop(pop_source).getComp(comp_label).popsize[tid_source] * converted_frac
                        self.getPop(pop_source).getComp(comp_label).popsize[tid_source] -= num_trans
                        self.getPop(pop_sink).getComp(comp_label).popsize[tid_sink] += num_trans
        
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
                        raise OptimaException('ERROR: Compartment or characteristic %s has not been pre-calculated for use in calculating %s.' % (inc_label, cid))
                        
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
                        raise OptimaException('ERROR: Compartment or characteristic %s has not been pre-calculated for use in calculating %s.' % (inc_label, cid))
                    
                    outputs[cid][pop.label] /= vals
        
        return outputs
            
        
    

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