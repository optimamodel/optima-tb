#%% Imports

from utils import printv, odict, OptimaException
from uuid import uuid4 as uuid
from numpy import array



#%% Model compartment classes

# Lightweight class to represent one compartment within a population.
class Node(object):
    def __init__(self, name = 'default', popsize = 0.0):
        self.name = name
        self.uid = uuid()
        self.popsize = float(popsize)
        
    # Link this node to another (i.e. create a transition link).
    def makeLinkTo(self, other_node):
        if not isinstance(other_node, Node):
            raise OptimaException('ERROR: Attempting to link compartment to something that is not a compartment.')
        return Link(self, other_node)

# Lightweight class to represent unidirectional flow between two compartments within a population.
class Link(object):
    def __init__(self, node1, node2, frac_transit = 0.0):
        self.name1 = node1.name
        self.name2 = node2.name
        self.uid1 = node1.uid
        self.uid2 = node2.uid
        self.frac_transit = float(frac_transit)   # Fraction of compartment 1 to move to compartment 2 per year.

# A class to wrap up data for one population within model.
class ModelPop(object): 
    def __init__(self, name = 'default'):
        self.name = name        
        self.nodes = odict()
        self.links = odict()
        
        self.genCascade()
        
    # Generate standard cascade, creating a node for each compartment and linking them appropriately.
    # NOTE: Consider generalising this method for future diseases using compartmental model.
    def genCascade(self):
        for name in ['sus','vac','lte','lts','ltf','dse','rec','dead']:
            self.nodes[name] = Node(name = name)
        for couple in [('sus','vac'), ('sus','lte'), ('vac','lte'), ('lte','lts'), ('lte','ltf'),
                     ('lts','dse'), ('ltf','dse'), ('dse','rec'), ('dse','dead')]:
            self.links[couple[0]+'->'+couple[1]] = self.nodes[couple[0]].makeLinkTo(self.nodes[couple[1]])
    
    # Evolve model population characteristics by one timestep (defaulting as 1 year).
    def stepForward(self, dt = 1.0):
        for link in self.links:
            pass
            
    # Loop through all nodes and print out current variable values.
    def printCurrVars(self):
        for oid in self.nodes:
            print('[Pop: %s][Compartment: %s][Popsize: %f]' % (self.name, self.nodes[oid].name, self.nodes[oid].popsize))
        
            

#%% Model function (simulates epidemic dynamics)

def model(verbose = 2):
    ''' Processes the TB epidemiological model. '''
    
    #%% Setup
    
    sim_settings = odict()
    sim_settings['t'] = array([2000])
    m_pops = odict()
    m_pops['kids'] = ModelPop(name = 'kids')

    #%% Run (i.e. evolve epidemic through time)

    for oid in m_pops:
        m_pops[oid].stepForward()
    
    #%% Collect and return raw results    
    
    return m_pops