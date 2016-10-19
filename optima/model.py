#%% Imports

from utils import flattenDict, odict, OptimaException

import numpy as np
import numpy.random as npr
from copy import deepcopy as dcp

import scipy.optimize as spo



#%% Calculation functions used in model

def transitionRelation(v_test, v_know, v_formats, dt):
    '''
    Solvable function that relates Poisson-based transition event rates with fractions of people moving after time dt.
    Input array v_know can contain any mix of rates and fractions, as long as corresponding list v_formats indicates which are marked by 'rate'.
    A numerical solver would then solve for v_test, a mix of fractions and rates opposite to v_know.
    '''
    
    f = np.zeros(len(v_know))
    F = np.zeros(len(v_know))
    d = np.zeros(len(v_know))
    
    # Rearrange arrays of values provided by solver and user into arrays of rates (f) and fractions (F).
    for k in xrange(len(v_know)):
        if v_formats[k] == 'rate':
            f[k] = v_know[k]
            F[k] = v_test[k]
        else:
            f[k] = v_test[k]
            F[k] = v_know[k]

    # As an array, calculate the discrepancy between calculated and given fraction values.
    # The two correspond if the difference is zero.
    if sum(f) == 0:
        d = -F
    else:
        d = (1-np.exp(-sum(f)*dt))*f/sum(f)-F
    
    return d
            

def convertTransitions(values, value_formats, old_dt, new_dt):
    ''' Function that converts yearly transition values provided in various formats (assuming old_dt is 1) into timestep-relevant fractions. '''
    
    new_vals = np.zeros(len(values))
    
    # Convert any values that are provided as probabilities into Poisson rates (average number of times transition event is encountered per year).
    k = 0
    for val in values:
        if value_formats[k] == 'probability':
            values[k] = -np.log(1-val)/old_dt
            value_formats[k] = 'rate'
        k += 1
    rates = dcp(values)     # This pre-allocated array should currently be a mix of fractions and rates, but will be made all rates later.
    
    # Given a set of rates/fractions, solve what the corresponding fractions/rates are.
    x = spo.fsolve(transitionRelation, np.ones(len(values))/2, args=(values, value_formats, old_dt), full_output = True)
    
    # If the numerical solver does not converge, print the solver output and crash.
    if x[2] != 1: 
        print x        
        raise OptimaException('ERROR: Transitions cannot be reconciled. This may be due to the sum of yearly outflows for a compartment being greater than 100%.')
    
    # If a value passed into the solver was not a rate, then the corresponding output should be. Finish generating a Poisson rates array.
    for k in xrange(len(values)):
        if value_formats[k] != 'rate':
            rates[k] = x[0][k]

    # Convert Poisson rates into fractions of compartment populations that should be moved per timestep new_dt.
    k = 0
    for rate in rates:
        if sum(rates) == 0:
            new_vals[k] = 0.0
        else:
            new_vals[k] = (1-np.exp(-sum(rates)*new_dt))*rates[k]/sum(rates)
        k += 1

    return new_vals



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
            self.comps.append(ModelCompartment(label = label, index = (self.index,k)))
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
    
#    def stepCascadeForward(self, dt = 1.0):
#        '''
#        Evolve model population characteristics by one timestep (defaulting as 1 year).
#        Calculated outputs like popsize will overwrite pre-allocated NaNs.
#        In contrast, pre-allocated arrays that are used in calculations (e.g. parameters) must be pre-filled with appropriate values.
#        This method only works if all links in this population target compartments in the same population. Generally used only for testing/debugging.
#        '''
#        
#        ti = self.t_index
#        
#        # If not pre-allocated, extend compartment variable arrays. Either way copy current value to the next position in the array.
#        for comp in self.comps:
#            if not len(comp.popsize) > ti + 1:
#                comp.popsize = np.append(comp.popsize, 0.0)
#            if not len(comp.popsize) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
#                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment %s.' % (comp.label))
#            comp.popsize[ti+1] = comp.popsize[ti]
#                
#        dpopsize = np.zeros(len(self.links))        
#        
#        # First calculation loop. Extend link variable arrays if not pre-allocated. Calculate value changes for next timestep.
#        for k, link in enumerate(self.links):
#            if not len(link.vals) > ti + 1:
#                link.vals = np.append(link.vals, link.vals[-1])
#            if not len(link.vals) > ti + 1:         # If one extension did not create an index of ti+1, something is seriously wrong...
#                raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment %s.' % (link.label))
#            if link.index_to[0] != self.index:      # If a link in this population targets a compartment in another population, crash this method.
#                raise OptimaException('ERROR: Population %s contains a link from %s to another population. Cannot step cascade independently forward.' % (self.label, link.label_from))
#            
#            converted_frac = 1 - (1 - link.vals[ti]) ** dt      # A formula for converting from yearly fraction values to the dt equivalent.
#            dpopsize[k] = self.getComp(link.label_from).popsize[ti] * converted_frac
#
#        # Second calculation loop. Apply value changes at next timestep.
#        for k, link in enumerate(self.links):
#            self.getComp(link.label_from).popsize[ti+1] -= dpopsize[k]
#            self.getComp(link.label_to).popsize[ti+1] += dpopsize[k]
#            
#        self.t_index += 1       # Update timestep index.
            
    def printCompVars(self, full = False):
        ''' Loop through all compartments and print out current variable values. '''
        for comp in self.comps:
            if not full:
                print('[Pop: %s][Compartment: %5s][Popsize: %15.4f]' % (self.label, comp.label, comp.popsize[self.t_index]))
            else:
                print('[Pop: %s][Compartment: %5s][Popsize...]' % (self.label, comp.label))
                print(comp.popsize)

#    def printLinkVars(self, full = False):
#        ''' Loop through all links and print out current variable values. '''
#        for link in self.links:
#            if not full:
#                print('[Pop: %s][%5s --> %-5s][Transit. Frac.: %5.4f]' % (self.label, link.label_from, link.label_to, link.vals[self.t_index]))
#            else:
#                print('[Pop: %s][%5s --> %-5s][Transit. Fraction...]' % (self.label, link.label_from, link.label_to))
#                print(link.vals)

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
            
#            self.pops[-1].getComp('sus').popsize[0] = 1000000   # NOTE: Temporary. Initial values inserted here.
        
        # Propagating initial characteristic parset values into ModelPops.
        # Critically must be applied in order defined within settings.charac_specs.
        init_dict = {}
        t_init = np.array([self.sim_settings['tvec'][0]])
        for charac_label in settings.charac_specs.keys():
            if settings.charac_specs[charac_label].has_key('entry_point'):
                entry_point = settings.charac_specs[charac_label]['entry_point']
                par = parset.pars['characs'][parset.par_ids['characs'][charac_label]]
                for pop_label in parset.pop_labels:
                    if not init_dict.has_key(pop_label): init_dict[pop_label] = {}
                    val = par.interpolate(tvec = t_init, pop_label = pop_label)
                    for include in flattenDict(input_dict = settings.charac_specs, base_key = charac_label, sub_key = 'includes'):
                        if init_dict[pop_label].has_key(include):
                            val -= init_dict[pop_label][include]
                    self.getPop(pop_label).getComp(entry_point).popsize[0] = val
                    init_dict[pop_label][entry_point] = val
        
        # Propagating cascade parameter parset values into ModelPops.
        for par in parset.pars['cascade']:
            tag = settings.linkpar_specs[par.label]['tag']          # Map parameter label -> link tag.
            for pop_label in parset.pop_labels:
                for link_id in self.getPop(pop_label).link_ids[tag]:           # Map link tag -> link id in ModelPop.            
                    self.getPop(pop_label).links[link_id].vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_label)
                    self.getPop(pop_label).links[link_id].val_format = par.y_format[pop_label]
        
        # Propagating transfer parameter parset values into Model object.
        for trans_type in parset.transfers.keys():
            if parset.transfers[trans_type]:
                for pop_source in parset.transfers[trans_type].keys():
                    par = parset.transfers[trans_type][pop_source]
                    for pop_target in par.y:
                        for comp in self.getPop(pop_source).comps:
                            trans_tag = comp.label + '_' + trans_type + '_to_' + pop_target       # NOTE: Perhaps there is a nicer way to set up transfer tagging.
                            if not comp.tag_dead:
                                num_links = len(self.getPop(pop_source).links)
                                link = comp.makeLinkTo(self.getPop(pop_target).getComp(comp.label),link_index=num_links)
                                link.vals = par.interpolate(tvec = self.sim_settings['tvec'], pop_label = pop_target)
                                link.val_format = par.y_format[pop_target]
                                
                                self.getPop(pop_source).links.append(link)
                                self.getPop(pop_source).link_ids[trans_tag] = [num_links]
                        
                
    def process(self, settings):
        ''' Run the full model. '''
        
        for t in self.sim_settings['tvec'][1:]:
            self.stepForward(dt = settings.tvec_dt)
        
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
        
        
    def stepForward(self, dt = 1.0):
        '''
        Evolve model characteristics by one timestep (defaulting as 1 year).
        Each application of this method writes calculated values to the next position in popsize arrays, regardless of dt.
        Thus the corresponding time vector associated with variable dt steps must be tracked externally.
        '''
        
        ti = self.t_index
        
        # Preallocate a change-in-popsize array. Requires each population to have the same cascade.
        num_pops = len(self.pops)
        num_comps = len(self.pops[0].comps)     # NOTE: Hard-coded check for zeroth population. Improve later.
        dpopsize = np.zeros(num_pops*num_comps)
        
        # First loop through all pops and comps to calculate value changes.
        k = 0
        for pop in self.pops:
            for comp in pop.comps:
                
                # If not pre-allocated, extend compartment variable arrays. Either way copy current value to the next position in the array.
                if not len(comp.popsize) > ti + 1:
                    comp.popsize = np.append(comp.popsize, 0.0)
                if not len(comp.popsize) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                    raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment %s.' % (comp.label))
                comp.popsize[ti+1] = comp.popsize[ti]
                
#                vals = np.zeros(comp.num_outlinks)
#                val_formats = [None]*comp.num_outlinks
#                j = 0
                for lid in comp.outlink_ids:
                    link = pop.links[lid]
                    
                    # If link values are not pre-allocated to the required time-vector length, extend with the same value as the last index.
                    if not len(link.vals) > ti + 1:
                        link.vals = np.append(link.vals, link.vals[-1])
                    if not len(link.vals) > ti + 1:         # If one extension did not create an index of ti+1, something is seriously wrong...
                        raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment %s.' % (link.label))
#                    
#                    # Store values and formats for outlinks relevant to the current compartment.
#                    vals[j] = link.vals[ti]
#                    val_formats[j] = link.val_format
#                    j += 1
#                
#                # If there are transitions to be applied, convert them to movement fractions appropriate to one timestep.
#                if len(vals) > 0:
#                    new_vals = convertTransitions(values = dcp(vals), value_formats = dcp(val_formats), old_dt = 1.0, new_dt = dt)
#                        
#                    j = 0
#                    for lid in comp.outlink_ids:
#                    link = pop.links[lid]
                    
                    did_from = link.index_from[0] * num_comps + link.index_from[1]
                    did_to = link.index_to[0] * num_comps + link.index_to[1]
                    comp_source = self.pops[link.index_from[0]].getComp(link.label_from)
                    
                    converted_frac = 1 - (1 - link.vals[ti]) ** dt      # A formula for converting from yearly fraction values to the dt equivalent.
                    
                    dpopsize[did_from] -= comp_source.popsize[ti] * converted_frac
                    dpopsize[did_to] += comp_source.popsize[ti] * converted_frac
#                    j += 1
                
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