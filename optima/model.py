#%% Imports

from utils import odict, OptimaException

import numpy as np
import numpy.random as npr
from copy import deepcopy as dcp



#%% Model compartment classes

class Node(object):
    ''' Lightweight class to represent one compartment within a population. '''
    def __init__(self, label = 'default', index = 0, popsize = 0.0):
        self.label = label          # Reference name for this object.
        self.index = index          # Positional index for container list.
        self.num_outlinks = 0       # Tracks number of nodes this one is linked to as an initial node.
        self.popsize = np.array([float(popsize)])   # Number of people in compartment.
        self.tag_dead = False       # Tag for whether this compartment contains dead people.
        
    def makeLinkTo(self, other_node):
        ''' Link this node to another (i.e. create a transition link).'''
        if not isinstance(other_node, Node):
            raise OptimaException('ERROR: Attempting to link compartment to something that is not a compartment.')
        self.num_outlinks += 1
        return Link(self, other_node)

class Link(object):
    '''
    Lightweight class to represent unidirectional flow.
    If used in ModelPop, the Link refers to two cascade compartments within a single population.
    If used in Model, the Link refers to two distinct population groups.
    In the latter case, intended logic should transfer agents between all (non-dead) corresponding compartments.
    '''
    def __init__(self, object_from, object_to, transit_frac = 0.0):
        self.index_from = object_from.index
        self.index_to = object_to.index
        self.label_from = object_from.label
        self.label_to = object_to.label
        self.transit_frac = np.array([float(transit_frac)])   # Fraction of compartment 1 to move to compartment 2 per year.

class ModelPop(object): 
    ''' A class to wrap up data for one population within model. '''
    def __init__(self, settings, label = 'default', index = 0):
        self.label = label          # Reference name for this object.
        self.index = index          # Positional index for container list.
        self.nodes = list()
        self.links = list()
        self.node_ids = dict()  # Maps node label to positional index in nodes list.
        self.link_ids = dict()  # Maps link tag to list of positional indices in links list.
        self.t_index = 0        # Keeps track of array index for current data within all nodes.
        
        self.genCascade(settings = settings)
    
    def getNode(self, node_label):
        ''' Allow nodes to be retrieved by label rather than index. '''
        node_index = self.node_ids[node_label]
        return self.nodes[node_index]
        
    # NOTE: Consider generalising this method for future diseases using compartmental model.
    def genCascade(self, settings):
        ''' Generate standard cascade, creating a node for each compartment and linking them appropriately. '''
        for l, label in enumerate(settings.node_specs.keys()):
            self.nodes.append(Node(label = label, index = l))
            if 'tag_dead' in settings.node_specs[label].keys():
                self.nodes[-1].tag_dead = True
            self.node_ids[label] = l
        l = 0
        for tag in settings.links.keys():
            for pair in settings.links[tag]:
#                node_from = self.node_ids[pair[0]]
#                node_to = self.node_ids[pair[1]]
                self.links.append(self.getNode(pair[0]).makeLinkTo(self.getNode(pair[1])))
                if not tag in self.link_ids:
                    self.link_ids[tag] = []
                self.link_ids[tag].append(l)
                l += 1
    
    def stepCascadeForward(self, dt = 1.0):
        '''
        Evolve model population characteristics by one timestep (defaulting as 1 year).
        Calculated outputs like popsize will overwrite pre-allocated NaNs.
        In contrast, pre-allocated arrays that are used in calculations (e.g. parameters) must be pre-filled with appropriate values.
        '''
        
        ti = self.t_index
        
        # If not pre-allocated, extend node variable arrays. Either way copy current value to the next position in the array.
        for node in self.nodes:
            if not len(node.popsize) > ti + 1:
                node.popsize = np.append(node.popsize, 0.0)
            if not len(node.popsize) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in node %s.' % (node.label))
            node.popsize[ti+1] = node.popsize[ti]
                
        dpopsize = np.zeros(len(self.links))        
        
        # First calculation loop. Extend link variable arrays if not pre-allocated. Calculate value changes for next timestep.
        for k, link in enumerate(self.links):
            if not len(link.transit_frac) > ti + 1:
                link.transit_frac = np.append(link.transit_frac, link.transit_frac[-1])
            if not len(link.transit_frac) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in link %s.' % (link.label))
            
#            node_from = self.nodes[link.index_from]
            dpopsize[k] = self.getNode(link.label_from).popsize[ti] * link.transit_frac[ti] * dt    # NOTE: Should this just be times dt...?

        # Second calculation loop. Apply value changes at next timestep.
        for k, link in enumerate(self.links):
            self.getNode(link.label_from).popsize[ti+1] -= dpopsize[k]
            self.getNode(link.label_to).popsize[ti+1] += dpopsize[k]
            
        self.t_index += 1       # Update timestep index.
            
    def printNodeVars(self, full = False):
        ''' Loop through all nodes and print out current variable values. '''
        for node in self.nodes:
            if not full:
                print('[Pop: %s][Node: %5s][Popsize: %15.4f]' % (self.label, node.label, node.popsize[self.t_index]))
            else:
                print('[Pop: %s][Node: %5s][Popsize...]' % (self.label, node.label))
                print(node.popsize)

    def printLinkVars(self, full = False):
        ''' Loop through all links and print out current variable values. '''
        for link in self.links:
            if not full:
                print('[Pop: %s][%5s --> %-5s][Transit. Frac.: %5.4f]' % (self.label, link.label_from, link.label_to, link.transit_frac[self.t_index]))
            else:
                print('[Pop: %s][%5s --> %-5s][Transit. Fraction...]' % (self.label, link.label_from, link.label_to))
                print(link.transit_frac)

    def makeRandomVars(self, for_node = True, for_link = True):
        ''' Randomise all node and link variables. Method used primarily for debugging. '''
        if for_node:
            for node in self.nodes:
                node.popsize[self.t_index] = npr.rand()*1e7
        if for_link:
            for link in self.links:
                link.transit_frac = npr.rand()/self.nodes[link.index_from].num_outlinks     # Scaling makes sure fractions leaving a node sum to less than 1.
                
    def preAllocate(self, sim_settings):
        '''
        Pre-allocate variable arrays in nodes and links for faster processing.
        Array maintains initial value but pre-fills everything else with NaNs.
        Thus errors due to incorrect parset value saturation should be obvious from results.
        '''
        for node in self.nodes:
            init_popsize = node.popsize[0]
            node.popsize = np.ones(len(sim_settings['tvec']))*np.nan
            node.popsize[0] = init_popsize
        for link in self.links:
            init_transit_frac = link.transit_frac[0]
            link.transit_frac = np.ones(len(sim_settings['tvec']))*np.nan
            link.transit_frac[0] = init_transit_frac
            
            
            
#%% Model class
            
class Model(object):
    ''' A class to wrap up multiple populations within model and handle cross-population transitions. '''
    
    def __init__(self):
        
        self.pops = list()
        self.transfers = list()
        self.pop_ids = dict()           # Maps pop label to positional index in pops list.
        self.transfer_ids = dict()      # Maps transfer label to positional index in transfer list.        
        
#        self.pops = odict()
        self.sim_settings = odict()
        
        self.transfers = odict()
        
        
        
    def getPop(self, pop_label):
        ''' Allow model populations to be retrieved by label rather than index. '''
        pop_index = self.pop_ids[pop_label]
        return self.pops[pop_index]
        
        
    def build(self, settings, parset):
        ''' Build the full model. '''
        self.sim_settings['tvec'] = np.arange(settings.tvec_start, settings.tvec_end + settings.tvec_dt/2, settings.tvec_dt)
        
        for k, pop_label in enumerate(parset.pop_labels):
#            pop_label = parset.pop_labels[k]
#            pop_name = parset.pop_names[k]
            self.pops.append(ModelPop(settings = settings, label = pop_label, index = k))
            self.pops[-1].preAllocate(self.sim_settings)     # Memory is allocated, speeding up model. However, values are NaN so as to enforce proper parset value saturation.
            self.pop_ids[pop_label] = k
            
            self.pops[-1].getNode('sus').popsize[0] = 1000000   # NOTE: Temporary. Initial values inserted here.
            
        # Transferring parset values into ModelPops.
        for par in parset.pars:
            tag = settings.linkpar_specs[par.label]['tag']          # Map parameter label -> link tag.
            for pop_label in parset.pop_labels:
#                pop_id = self.pop_ids[pop_label]
                for link_id in self.getPop(pop_label).link_ids[tag]:           # Map link tag -> link id in ModelPop.            
                    self.getPop(pop_label).links[link_id].transit_frac = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_label)
        
        self.transfers = dcp(parset.transfers)
                
    def process(self, settings):
        ''' Run the full model. '''
        
        for t in self.sim_settings['tvec'][1:]:
            for pop in self.pops:
                pop.stepCascadeForward(dt = settings.tvec_dt)
                
            # Apply transfers.
            for trans_type in self.transfers.keys():
                for pop_source in self.transfers[trans_type].keys():
                    for pop_sink in self.transfers[trans_type][pop_source].y:
#                        pid_source = self.pop_ids[pop_source]
#                        pid_sink = self.pop_ids[pop_sink]
                        transit_frac = self.transfers[trans_type][pop_source].y[pop_sink][0]
                        tid_source = self.getPop(pop_source).t_index
                        tid_sink = self.getPop(pop_sink).t_index
                        for node_label in self.getPop(pop_source).node_ids.keys():
#                            nid_source = self.getPop(pop_source).node_ids[node_label]
                            if not self.getPop(pop_source).getNode(node_label).tag_dead:
#                                nid_sink = self.getPop(pop_sink).node_ids[node_label]
                                num_trans = self.getPop(pop_source).getNode(node_label).popsize[tid_source] * transit_frac * settings.tvec_dt
                                self.getPop(pop_source).getNode(node_label).popsize[tid_source] -= num_trans
                                self.getPop(pop_sink).getNode(node_label).popsize[tid_sink] += num_trans
                
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
#                pid = self.pop_ids[pop.label]
                outputs[cid][pop.label] = None
                
                # Sum up all relevant node popsizes (or previously calculated characteristics).
                for inc_label in settings.charac_specs[cid]['includes']:
                    if inc_label in self.getPop(pop.label).node_ids.keys():
#                        nid = self.getPop(pop.label).node_ids[inc_label]
                        vals = self.getPop(pop.label).getNode(inc_label).popsize
                    elif inc_label in outputs.keys()[:-1]:
                        vals = outputs[inc_label][pop.label]
                    else:
                        raise OptimaException('ERROR: Node or characteristic %s has not been pre-calculated for use in calculating %s.' % (inc_label, cid))
                        
                    if outputs[cid][pop.label] is None:
                        outputs[cid][pop.label] = dcp(vals)
                    else:
                        outputs[cid][pop.label] += vals
                
                # Divide by relevant node popsize (or previously calculated characteristic).
                if 'denom' in settings.charac_specs[cid]:
                    den_label = settings.charac_specs[cid]['denom']
                    if den_label in outputs.keys()[:-1]:
                        vals = outputs[den_label][pop.label]
                    elif den_label in self.getPop(pop.label).node_ids.keys():
#                        nid = self.getPop(pop.label).node_ids[den_label]
                        vals = self.getPop(pop.label).getNode(den_label).popsize
                    else:
                        raise OptimaException('ERROR: Node or characteristic %s has not been pre-calculated for use in calculating %s.' % (inc_label, cid))
                    
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