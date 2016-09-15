#%% Imports

from utils import odict, OptimaException

#from numpy import array, append, zeros, ones, arange
#from numpy.random import rand
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
        self.transit_frac = float(transit_frac)   # Fraction of compartment 1 to move to compartment 2 per year.

class ModelPop(object): 
    ''' A class to wrap up data for one population within model. '''
    def __init__(self, settings, name = 'default'):
        self.name = name        
        self.nodes = list()
        self.links = list()
        self.node_ids = dict()  # Maps node label to positional index in nodes list.
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
            self.nodes.append(Node(name=label, index=l))
            self.node_ids[label] = l
        for pair in settings.links[:]:
            node_from = self.node_ids[pair[0]]
            node_to = self.node_ids[pair[1]]
            self.links.append(self.nodes[node_from].makeLinkTo(self.nodes[node_to]))
    
    def stepForward(self, dt = 1.0):
        ''' Evolve model population characteristics by one timestep (defaulting as 1 year). '''
        
        ti = self.t_index
        
        # If not pre-allocated, extend node variable arrays. Either way copy current value to the next position in the array.
        for node in self.nodes:
            if not len(node.popsize) > ti + 1:
                node.popsize = np.append(node.popsize, 0.0)
            if not len(node.popsize) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in node %s.' % (node.name))
            node.popsize[ti+1] = node.popsize[ti]
        
        dpopsize = np.zeros(len(self.links))        
        
        # First loop. Calculate value changes for next timestep.
        for k,link in enumerate(self.links):
            node_from = self.nodes[link.index_from]
            dpopsize[k] = node_from.popsize[ti] * link.transit_frac * dt    # NOTE: Should this just be times dt...?

        # Second loop. Apply value changes at next timestep.
        for k,link in enumerate(self.links):
            self.nodes[link.index_from].popsize[ti+1] -= dpopsize[k]
            self.nodes[link.index_to].popsize[ti+1] += dpopsize[k]
            
        self.t_index += 1       # Update timestep index.
            
    def printNodeVars(self, full=False):
        ''' Loop through all nodes and print out current variable values. '''
        for node in self.nodes:
            if not full:
                print('[Pop: %s][Node: %5s][Popsize: %15.4f]' % (self.name, node.name, node.popsize[self.t_index]))
            else:
                print('[Pop: %s][Node: %5s][Popsize...]' % (self.name, node.name))
                print(node.popsize)

    def printLinkVars(self):
        ''' Loop through all nodes and print out current variable values. '''
        for link in self.links:
            print('[Pop: %s][%5s --> %-5s][Transit. Frac.: %5.4f]' % (self.name, link.name_from, link.name_to, link.transit_frac))

    def makeRandomVars(self, for_node = True, for_link = True):
        ''' Randomise all node and link variables. Method used primarily for debugging. '''
        if for_node:
            for node in self.nodes:
                node.popsize[self.t_index] = npr.rand()*1e7
        if for_link:
            for link in self.links:
                link.transit_frac = npr.rand()/self.nodes[link.index_from].num_outlinks     # Scaling makes sure fractions leaving a node sum to less than 1.
                
    def preAllocate(self, sim_settings):
        ''' Pre-allocate variable arrays in nodes for faster processing. Array is pre-filled with initial value. '''
        for node in self.nodes:
            init_popsize = node.popsize[0]
            node.popsize = np.ones(len(sim_settings['tvec']))*init_popsize
            

#%% Model function (simulates epidemic dynamics)

def model(settings):
    ''' Processes the TB epidemiological model. '''
    
    #%% Setup
    
    sim_settings = odict()
    dt = 0.25
    sim_settings['tvec'] = np.arange(2000, 2030+dt/2, dt)
    
    m_pops = odict()
    m_pops['kids'] = ModelPop(settings = settings, name = 'kids')
    m_pops['adults'] = ModelPop(settings = settings, name = 'adults')
    
    for oid in m_pops:
        m_pops[oid].makeRandomVars(for_node = False, for_link = True)
        m_pops[oid].getnode('sus').popsize[0] = 1000000
        m_pops[oid].preAllocate(sim_settings)

    #%% Run (i.e. evolve epidemic through time)

    for oid in m_pops:
        for t in sim_settings['tvec'][1:]:
#            print('Time: %.1f' % t)
            m_pops[oid].stepForward(dt = 0.25)
#        m_pops[oid].printLinkVars()
#        m_pops[oid].printNodeVars(full = True)
    
    #%% Collect and return raw results    
    
    return m_pops, sim_settings