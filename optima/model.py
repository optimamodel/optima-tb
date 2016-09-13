#%% Imports

from utils import printv, odict, OptimaException
from uuid import uuid4 as uuid
from numpy import array, append, zeros, arange
from numpy.random import rand



#%% Model compartment classes

# Lightweight class to represent one compartment within a population.
class Node(object):
    def __init__(self, name = 'default', popsize = 0.0):
        self.name = name
        self.num_outlinks = 0       # Tracks number of nodes this one is linked to as an initial node.
#        self.uid = uuid()
        self.popsize = array([float(popsize)])   # Current number of people in compartment.
#        self.pre_popsize = array([])    # Values of popsize for previous timesteps.
        
    # Link this node to another (i.e. create a transition link).
    def makeLinkTo(self, other_node):
        if not isinstance(other_node, Node):
            raise OptimaException('ERROR: Attempting to link compartment to something that is not a compartment.')
        self.num_outlinks += 1
        return Link(self, other_node)

# Lightweight class to represent unidirectional flow between two compartments within a population.
class Link(object):
    def __init__(self, node1, node2, transit_frac = 0.0):
        self.name1 = node1.name
        self.name2 = node2.name
#        self.uid1 = node1.uid
#        self.uid2 = node2.uid
        self.transit_frac = float(transit_frac)   # Fraction of compartment 1 to move to compartment 2 per year.

# A class to wrap up data for one population within model.
class ModelPop(object): 
    def __init__(self, name = 'default'):
        self.name = name        
        self.nodes = odict()
        self.links = odict()
        self.t_index = 0        # Keeps track of array index for current data within all nodes.
        
        self.genCascade()
        
    # Generate standard cascade, creating a node for each compartment and linking them appropriately.
    # NOTE: Consider generalising this method for future diseases using compartmental model.
    def genCascade(self):
        for name in ['sus','vac','lte','lts','ltf','dse','rec','dead']:
            self.nodes[name] = Node(name = name)
        for specs in [('sus','vac'), ('sus','lte'), ('vac','lte'), ('lte','lts'), ('lte','ltf'), ('lts','sus'), ('ltf','sus'), ('lts','dse'), ('ltf','dse'), ('dse','rec'), ('dse','dead'), ('rec','lte')]:
            self.links[specs[0]+'->'+specs[1]] = self.nodes[specs[0]].makeLinkTo(self.nodes[specs[1]])
    
    # Evolve model population characteristics by one timestep (defaulting as 1 year).
    def stepForward(self, dt = 1.0):
        dpopsize = zeros(len(self.links))
        for k in xrange(len(self.links)):
            dpopsize[k] = self.nodes[self.links[k].name1].popsize[self.t_index] * self.links[k].transit_frac * dt
            
            # If not pre-allocated, extend the popsize array. Either way copy current value to the next position in the array.
            if not len(self.nodes[self.links[k].name1].popsize) > self.t_index + 1:
                self.nodes[self.links[k].name1].popsize = append(self.nodes[self.links[k].name1].popsize, 0.0)
            if not len(self.nodes[self.links[k].name2].popsize) > self.t_index + 1:
                self.nodes[self.links[k].name2].popsize = append(self.nodes[self.links[k].name2].popsize, 0.0)
            if len(self.nodes[self.links[k].name1].popsize) <= self.t_index + 1 or len(self.nodes[self.links[k].name2].popsize) <= self.t_index + 1:
                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in a node.')
                
            self.nodes[self.links[k].name1].popsize[self.t_index+1] = self.nodes[self.links[k].name1].popsize[self.t_index]
            self.nodes[self.links[k].name2].popsize[self.t_index+1] = self.nodes[self.links[k].name2].popsize[self.t_index]
        
        # Move to the next timestep and adjust values.
        self.t_index += 1
        for k in xrange(len(self.links)):
            self.nodes[self.links[k].name1].popsize[self.t_index] -= dpopsize[k]
            self.nodes[self.links[k].name2].popsize[self.t_index] += dpopsize[k]
            
    # Loop through all nodes and print out current variable values.
    def printNodeVars(self, full = False):
        for oid in self.nodes:
            if not full:
                print('[Pop: %s][Node: %5s][Popsize: %15.4f]' % (self.name, self.nodes[oid].name, self.nodes[oid].popsize[self.t_index]))
            else:
                print('[Pop: %s][Node: %5s][Popsize...]' % (self.name, self.nodes[oid].name))
                print(self.nodes[oid].popsize)

    # Loop through all nodes and print out current variable values.
    def printLinkVars(self):
        for oid in self.links:
            print('[Pop: %s][%5s --> %-5s][Transit. Frac.: %5.4f]' % (self.name, self.links[oid].name1, self.links[oid].name2, self.links[oid].transit_frac))

    # Randomise all node and link variables. Method used primarily for debugging.
    def makeRandomVars(self, for_node = True, for_link = True):
        if for_node:
            for oid in self.nodes:
                self.nodes[oid].popsize[self.t_index] = rand()*1e7
        if for_link:
            for oid in self.links:
                self.links[oid].transit_frac = rand()/self.nodes[self.links[oid].name1].num_outlinks     # Scaling makes sure fractions leaving a node sum to less than 1.
                
    # Pre-allocate variable arrays in nodes for faster processing.
    def preAllocate(self, sim_settings):
        for oid in self.nodes:
            self.popsize = zeros(len(sim_settings['t_range']))
            

#%% Model function (simulates epidemic dynamics)

def model(verbose = 2):
    ''' Processes the TB epidemiological model. '''
    
    #%% Setup
    
    sim_settings = odict()
    sim_settings['t_range'] = arange(2017,2031)
    
    m_pops = odict()
    m_pops['kids'] = ModelPop(name = 'kids')
    m_pops['adults'] = ModelPop(name = 'adults')
    
    for oid in m_pops:
        m_pops[oid].makeRandomVars(for_node = False, for_link = True)
        m_pops[oid].printLinkVars()
        m_pops[oid].nodes['sus'].popsize[0] = 1000000
        m_pops[oid].preAllocate(sim_settings)

    #%% Run (i.e. evolve epidemic through time)

    for oid in m_pops:
        m_pops[oid].makeRandomVars(for_node = False, for_link = True)
        for t in sim_settings['t_range'][1:]:
            print('Time: %.1f' % t)
            m_pops[oid].printNodeVars()
            m_pops[oid].stepForward()
    
    #%% Collect and return raw results    
    
    return m_pops