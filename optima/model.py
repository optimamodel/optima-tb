#%% Imports

from utils import odict, OptimaException

import numpy as np
import numpy.random as npr



#%% Model compartment classes

class Node(object):
    ''' Lightweight class to represent one compartment within a population. '''
    def __init__(self, name='default', index=0, popsize=0.0):
        self.name = name
        self.index = index
        self.num_outlinks = 0       # Tracks number of nodes this one is linked to as an initial node.
        self.popsize = np.array([float(popsize)])   # Number of people in compartment.
        
    def makeLinkTo(self, other_node):
        ''' Link this node to another (i.e. create a transition link).'''
        if not isinstance(other_node, Node):
            raise OptimaException('ERROR: Attempting to link compartment to something that is not a compartment.')
        self.num_outlinks += 1
        return Link(self, other_node)

class Link(object):
    ''' Lightweight class to represent unidirectional flow between two compartments within a population. '''
    def __init__(self, node_from, node_to, transit_frac = 0.0):
        self.index_from = node_from.index
        self.index_to = node_to.index
        self.name_from = node_from.name
        self.name_to = node_to.name
        self.transit_frac = np.array([float(transit_frac)])   # Fraction of compartment 1 to move to compartment 2 per year.

class ModelPop(object): 
    ''' A class to wrap up data for one population within model. '''
    def __init__(self, settings, name = 'default'):
        self.name = name      
        self.nodes = list()
        self.links = list()
        self.node_ids = dict()  # Maps node label to positional index in nodes list.
        self.link_ids = dict()  # Maps link tag to list of positional indices in links list.
        self.t_index = 0        # Keeps track of array index for current data within all nodes.
        
        self.genCascade(settings = settings)
    
    def getnode(self, node_name):
        ''' Allow nodes to be retrieved by name rather than index '''
        node_index = self.node_ids[node_name]
        return self.nodes[node_index]
        
    # NOTE: Consider generalising this method for future diseases using compartmental model.
    def genCascade(self, settings):
        ''' Generate standard cascade, creating a node for each compartment and linking them appropriately. '''
        for l, label in enumerate(settings.node_labels):
            self.nodes.append(Node(name = label, index = l))
            self.node_ids[label] = l
        l = 0
        for tag in settings.links.keys():
            for pair in settings.links[tag]:
                node_from = self.node_ids[pair[0]]
                node_to = self.node_ids[pair[1]]
                self.links.append(self.nodes[node_from].makeLinkTo(self.nodes[node_to]))
                if not tag in self.link_ids:
                    self.link_ids[tag] = []
                self.link_ids[tag].append(l)
                l += 1
    
    def stepCascadeForward(self, dt = 1.0):
        '''
        Evolve model population characteristics by one timestep (defaulting as 1 year).
        Calculation outputs like popsize will be propagated across pre-allocated NaNs.
        Other than those, pre-allocated arrays must have appropriate values.
        '''
        
        ti = self.t_index
        
        # If not pre-allocated, extend node variable arrays. Either way copy current value to the next position in the array.
        for node in self.nodes:
            if not len(node.popsize) > ti + 1:
                node.popsize = np.append(node.popsize, 0.0)
            if not len(node.popsize) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in node %s.' % (node.name))
            node.popsize[ti+1] = node.popsize[ti]
                
        dpopsize = np.zeros(len(self.links))        
        
        # First calculation loop. Extend link variable arrays if not pre-allocated. Calculate value changes for next timestep.
        for k,link in enumerate(self.links):
            if not len(link.transit_frac) > ti + 1:
                link.transit_frac = np.append(link.transit_frac, link.transit_frac[-1])
            if not len(link.transit_frac) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in link %s.' % (link.name))
            
            node_from = self.nodes[link.index_from]
            dpopsize[k] = node_from.popsize[ti] * link.transit_frac[ti] * dt    # NOTE: Should this just be times dt...?

        # Second calculation loop. Apply value changes at next timestep.
        for k,link in enumerate(self.links):
            self.nodes[link.index_from].popsize[ti+1] -= dpopsize[k]
            self.nodes[link.index_to].popsize[ti+1] += dpopsize[k]
            
        self.t_index += 1       # Update timestep index.
            
    def printNodeVars(self, full = False):
        ''' Loop through all nodes and print out current variable values. '''
        for node in self.nodes:
            if not full:
                print('[Pop: %s][Node: %5s][Popsize: %15.4f]' % (self.name, node.name, node.popsize[self.t_index]))
            else:
                print('[Pop: %s][Node: %5s][Popsize...]' % (self.name, node.name))
                print(node.popsize)

    def printLinkVars(self, full = False):
        ''' Loop through all links and print out current variable values. '''
        for link in self.links:
            if not full:
                print('[Pop: %s][%5s --> %-5s][Transit. Frac.: %5.4f]' % (self.name, link.name_from, link.name_to, link.transit_frac[self.t_index]))
            else:
                print('[Pop: %s][%5s --> %-5s][Transit. Fraction...]' % (self.name, link.name_from, link.name_to))
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
        self.pops = odict()
        self.sim_settings = odict()
        
        self.transfers = odict()
        
        
    def build(self, settings, parset):
        ''' Build the full model. '''
        self.sim_settings['tvec'] = np.arange(settings.tvec_start, settings.tvec_end + settings.tvec_dt/2, settings.tvec_dt)
        
        for k in xrange(len(parset.pop_labels)):
            pop_label = parset.pop_labels[k]
            pop_name = parset.pop_names[k]
            self.pops[pop_label] = ModelPop(settings = settings, name = pop_name)
            self.pops[pop_label].preAllocate(self.sim_settings)     # Memory is allocated, speeding up model. However, values are NaN so as to enforce proper parset value saturation.
            
            self.pops[pop_label].getnode('sus').popsize[0] = 1000000
            
        # Transferring parset values into ModelPops.
        for par in parset.pars:
            tag = settings.linkpar_specs[par.label]['tag']          # Map parameter label -> link tag.
            for pop_label in parset.pop_labels:
                for link_id in self.pops[pop_label].link_ids[tag]:           # Map link tag -> link id in ModelPop.            
                    self.pops[pop_label].links[link_id].transit_frac = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_label)
        
        for trans_type in parset.transfers:
            for pop_label in parset.transfers[trans_type].keys():
                self.transfers[pop_label] = parset.transfers[trans_type][pop_label]
                
                
            
    def process(self, settings, parset):
        ''' Run the full model. '''
        
        for t in self.sim_settings['tvec'][1:]:
            for pop_label in self.pops:
                self.pops[pop_label].stepCascadeForward(dt = settings.tvec_dt)
                
            # Apply transfers.
            for pop_source in self.transfers.keys():
                pop_sink = self.transfers[pop_source]['target']
                transit_frac = self.transfers[pop_source]['value']
                tid_source = self.pops[pop_source].t_index
                tid_sink = self.pops[pop_sink].t_index
                for node_label in self.pops[pop_source].node_ids.keys():
                    nid_source = self.pops[pop_source].node_ids[node_label]
                    nid_sink = self.pops[pop_sink].node_ids[node_label]
                    num_trans = self.pops[pop_source].nodes[nid_source].popsize[tid_source] * transit_frac * settings.tvec_dt
                    self.pops[pop_source].nodes[nid_source].popsize[tid_source] -= num_trans
                    self.pops[pop_sink].nodes[nid_sink].popsize[tid_sink] += num_trans
                
        return self.pops, self.sim_settings
        
    

#%% Model function (simulates epidemic dynamics)

def runModel(settings, parset):
    ''' Processes the TB epidemiological model. '''
    
    m = Model()
    m.build(settings = settings, parset = parset)
    m_pops, sim_settings = m.process(settings = settings, parset = parset)
    
    
    
    #%% Collect and return raw results    
    
    return m_pops, sim_settings