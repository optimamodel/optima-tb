#%% Imports

from utils import odict, OptimaException
from settings import Settings
from plotting import gridColorMap

from numpy import array, append, zeros, arange
from numpy.random import rand
from pylab import subplots#, fill_between
from copy import deepcopy as dcp



#%% Model compartment classes

# Lightweight class to represent one compartment within a population.
class Node(object):
    def __init__(self, name = 'default', popsize = 0.0):
        self.name = name
        self.num_outlinks = 0       # Tracks number of nodes this one is linked to as an initial node.
        self.popsize = array([float(popsize)])   # Number of people in compartment.
        
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
        self.transit_frac = float(transit_frac)   # Fraction of compartment 1 to move to compartment 2 per year.

# A class to wrap up data for one population within model.
class ModelPop(object): 
    def __init__(self, settings, name = 'default'):
        self.name = name        
        self.nodes = odict()
        self.links = odict()
        self.t_index = 0        # Keeps track of array index for current data within all nodes.
        
        self.genCascade(settings = settings)
        
    # Generate standard cascade, creating a node for each compartment and linking them appropriately.
    # NOTE: Consider generalising this method for future diseases using compartmental model.
    def genCascade(self, settings):
        for label in settings.node_labels:
            self.nodes[label] = Node(name = label)
        for specs in [('sus','vac'), ('sus','lte'), ('vac','lte'), ('lte','ltsu'), ('lte','ltfu'), ('ltsu','sus'), ('ltfu','sus'), ('ltsu','s+e'), ('ltfu','s+e'), ('s+e','rec'), ('s+e','dead'), ('rec','lte')]:
            self.links[specs[0]+'->'+specs[1]] = self.nodes[specs[0]].makeLinkTo(self.nodes[specs[1]])
    
    # Evolve model population characteristics by one timestep (defaulting as 1 year).
    def stepForward(self, dt = 1.0):
        
        ti = self.t_index
        
        # If not pre-allocated, extend node variable arrays. Either way copy current value to the next position in the array.        
        for k in xrange(len(self.nodes)):
            node = self.nodes[k]
            if not len(node.popsize) > ti + 1:
                node.popsize = append(node.popsize, 0.0)
            if not len(node.popsize) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in node %s.' % (node.name))
            node.popsize[ti+1] = node.popsize[ti]
        
        dpopsize = zeros(len(self.links))        
        
        # First loop. Calculate value changes for next timestep.
        for k in xrange(len(self.links)):
            link = self.links[k]
            node1 = self.nodes[link.name1]
            dpopsize[k] = node1.popsize[ti] * link.transit_frac * dt    # NOTE: Should this just be times dt...?

        # Second loop. Apply value changes at next timestep.
        for k in xrange(len(self.links)):
            self.nodes[self.links[k].name1].popsize[ti+1] -= dpopsize[k]
            self.nodes[self.links[k].name2].popsize[ti+1] += dpopsize[k]
            
        self.t_index += 1       # Update timestep index.
            
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
            self.popsize = zeros(len(sim_settings['tvec']))
            

#%% Model function (simulates epidemic dynamics)

def model(settings):
    ''' Processes the TB epidemiological model. '''
    
    #%% Setup
    
    sim_settings = odict()
    sim_settings['tvec'] = arange(2017,2051)
    
    m_pops = odict()
    m_pops['kids'] = ModelPop(settings = settings, name = 'kids')
    m_pops['adults'] = ModelPop(settings = settings, name = 'adults')
    
    for oid in m_pops:
        m_pops[oid].makeRandomVars(for_node = False, for_link = True)
        m_pops[oid].printLinkVars()
        m_pops[oid].nodes['sus'].popsize[0] = 1000000
        m_pops[oid].preAllocate(sim_settings)

    #%% Run (i.e. evolve epidemic through time)

    for oid in m_pops:
        for t in sim_settings['tvec'][1:]:
            print('Time: %.1f' % t)
            m_pops[oid].printNodeVars()
            m_pops[oid].stepForward()
    
    #%% Collect and return raw results    
    
    return m_pops, sim_settings



#%% Test model function

settings = Settings()
test_pops, sim_settings = model(settings = settings)

for pop_oid in test_pops:
    pop = test_pops[pop_oid]
    
    fig, ax = subplots(figsize=(15,10))
    colors = gridColorMap(len(pop.nodes))
    bottom = 0*sim_settings['tvec']
    
    for k in xrange(len(pop.nodes)):
        node = pop.nodes[k]
        top = bottom + node.popsize
        
        ax.fill_between(sim_settings['tvec'], bottom, top, facecolor=colors[k], alpha=1, lw=0)
        ax.plot((0, 0), (0, 0), color=colors[k], linewidth=10)
        bottom = dcp(top)
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])   
    
    legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5)}
    ax.set_title('Cascade - %s' % (pop.name))
    ax.set_xlabel('Year')
    ax.set_ylabel('People')
    ax.set_xlim((sim_settings['tvec'][0], sim_settings['tvec'][-1]))
    ax.set_ylim((0, max(top)))
    cascadenames = [pop.nodes[oid].name for oid in pop.nodes]
    ax.legend(cascadenames, **legendsettings)