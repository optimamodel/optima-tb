# %% Imports

from optima_tb.utils import flattenDict, odict, OptimaException
from optima_tb.validation import checkTransitionFraction
import optima_tb.settings as project_settings
from optima_tb.results import ResultSet
from optima_tb.parsing import FunctionParser
from optima_tb.ModelPrograms import ModelProgramSet, ModelProgram

import logging
logger = logging.getLogger(__name__)
parser = FunctionParser(debug=False)  # Decomposes and evaluates functions written as strings, in accordance with a grammar defined within the parser object.

import numpy as np
from copy import deepcopy as dcp
import uuid

import matplotlib.pyplot as plt

# np.seterr(all='raise')
# %% Abstract classes used in model

class ModelCompartment(object):
    pass


class Node(object):
    ''' Lightweight abstract class to represent one object within a network. '''
    def __init__(self, label='default'):
        self.uid = uuid.uuid4()
        self.label = label          # Reference name for this object.
        # self.index = index          # Index to denote storage position in Model object. Format is left unspecified.
        #                             # For a ModelPopulation, this is intended as an integer.
        #                             # For a Compartment, this is intended as a tuple with first element denoting index of ModelPopulation.

class Variable(object):
    '''
    Lightweight abstract class to store a variable array of values (presumably corresponding to an external time vector).
    Includes an attribute to describe the format of these values.

    Examples include characteristics and dependent parameters (NB. all non-dependents parameter correspond to links)
    '''
    def __init__(self, label='default'):
        self.uid = uuid.uuid4()
        self.label = label
        self.t = None
        self.vals = None
        self.vals_old = None                # An optional array that stores old values in the case of overwriting.
        self.units = ''

    def preallocate(self,tvec):
        self.t = tvec
        self.vals = np.ones(tvec.shape) * np.nan
        self.vals_old = np.ones(tvec.shape) * np.nan

    def plot(self):
        plt.figure()
        d = self.__dict__
        for label,val in d.items():
            if isinstance(val,np.ndarray) and val is not self.t and val.size == self.t.size:
                plt.plot(self.t,val,label=label)
        plt.legend()
        plt.xlabel('Year')
        plt.ylabel("%s (%s)" % (self.label,self.units))

    def update(self,ti=None):
        # A Variable can have a function to update its value at a given time, which is
        # overloaded differently for Characteristics and Parameters
        return

    def unlink(self):
        # A Variable has an unlink function that replaces all of its internal references 
        # with the corresponding object UUID
        for prop in self.__dict__.keys():
            if isinstance(self.__dict__[prop],list):
                for i,obj in enumerate(self.__dict__[prop]):
                    if hasattr(obj,'uid'):
                        self.__dict__[prop][i] = obj.uid
            else:
                obj = self.__dict__[prop]
                if hasattr(obj,'uid'):
                    self.__dict__[prop] = obj.uid

    def relink(self,objs):
        # Given a dictionary of objects, restore the internal references
        # based on the UUID
        for prop in self.__dict__.keys():
            if prop == 'uid':
                continue
            elif isinstance(self.__dict__[prop],list):
                for i,obj in enumerate(self.__dict__[prop]):
                    if isinstance(obj,uuid.UUID):
                        self.__dict__[prop][i] = objs[obj]
            else:
                obj = self.__dict__[prop]
                if isinstance(obj,uuid.UUID):
                    self.__dict__[prop] = objs[obj]

    def __repr__(self):
        return '%s "%s" (%s)' % (self.__class__.__name__,self.label,self.uid)

class Compartment(Variable):
    ''' A class to wrap up data for one compartment within a cascade network. '''

    def __init__(self, label='default'):
        Variable.__init__(self, label=label)
        self.units = 'people'
        self.tag_birth = False                      # Tag for whether this compartment contains unborn people.
        self.tag_dead = False                       # Tag for whether this compartment contains dead people.
        self.is_junction = False

        self.outlinks = []
        self.inlinks = []

    def __repr__(self, *args, **kwargs):
        return "Compartment %s: %g" % (self.label, self.vals[0])

    def getValue(self, ti):
        """ Get value of population at timestep ti """
        return self.vals[ti]

class Characteristic(Variable):
    ''' A characteristic represents a grouping of compartments 
    '''
    def __init__(self, label='default'):
        # includes is a list of Compartments, whose values are summed
        # the denominator is another Characteristic that normalizes this one
        # All passed by reference so minimal performance impact
        Variable.__init__(self, label=label)
        self.units = 'people'
        self.includes = []
        self.denominator = None
        self.dependency = False # This flag indicates whether another variable depends on this one, indicating the value needs to be computed during integration

    def add_include(self,x):
        assert isinstance(x,Compartment) or isinstance(x,Characteristic)
        self.includes.append(x)
        if isinstance(x,Characteristic):
            x.dependency = True

    def add_denom(self,x):
        assert isinstance(x,Compartment) or isinstance(x,Characteristic)
        self.denominator = x
        if isinstance(x,Characteristic):
            x.dependency = True
        self.units = 'proportion'


    def update(self,ti=None):
        # Read popsizes at time ti from includes, and update the value of this characteristic
        # If ti is none, then use all vals
        if ti is None:
            ti = np.arange(0,self.vals.size) # This corresponds to every time point
        else:
            ti = np.array(ti).ravel()

        self.vals[ti] = 0

        for comp in self.includes:
            self.vals[ti] += comp.vals[ti]

        if self.denominator is not None:
            denom  = self.denominator.vals[[ti]]

            ## Todo - do this better!
            for i in xrange(0,ti.size):
                if denom[i] > 0:
                    self.vals[ti[i]] /= denom[i]
                elif self.vals[ti[i]] < project_settings.TOLERANCE:
                    self.vals[ti[i]] = 0  # Given a zero/zero case, make the answer zero.
                else:
                    self.vals[ti[i]] = np.inf  # Given a non-zero/zero case, keep the answer infinite.
    
class Parameter(Variable):
    # A parameter is a Variable that can have a value computed via an f_stack and a list of 
    # dependent Variables. This class may need to be renamed to avoid confusion with
    # the class in Parameter.py which provides a means of computing the Parameter that is used by the model.
    # This is a Parameter in the cascade.xlsx sense - there is one Parameter object for every item in the 
    # Parameters sheet. A parameter that maps to multiple transitions (e.g. doth_rate) will have one parameter
    # and multiple Link instances that depend on the same Parameter instance
    def __init__(self, label='default',f_stack = None, deps = None, limits = None):
        Variable.__init__(self, label=label)
        if deps is not None:
            for dep in deps:
                if hasattr(dep,'dependency'):
                    # Mark all required objects as dependent if they have such a property
                    # This will happen for Characteristics and other Parameters, but not for Compartments
                    dep.dependency = True 

        self.deps = deps
        self.f_stack = f_stack
        self.limits = limits
        self.dependency = False 
        self.scale_factor = 1.0
        self.links = [] # References to links that derive from this parameter
        self.source_popsize_cache_time = None
        self.source_popsize_cache_val = None

    def constrain(self,ti):
        # NB. Must be an array, so ti must must not be supplied
        if self.limits is not None:
            self.vals[ti] = max(self.limits[0],self.vals[ti])
            self.vals[ti] = min(self.limits[1],self.vals[ti])

    def update(self,ti=None):
        # Update the value of this program at time index ti 
        # by evaluating its f_stack function using the 
        # current values of all dependent variables at time index ti
        
        if ti is None:
            ti = np.arange(0,self.vals.size) # This corresponds to every time point
        else:
            ti = np.array(ti)

        if self.f_stack is None:
            return

        dep_vals = {}
        for dep in self.deps:
            dep_vals[dep.label] = dep.vals[[ti]]
        self.vals[ti] = parser.evaluateStack(stack=self.f_stack[0:], deps=dep_vals)   # self.f_stack[0:] makes a copy

    def source_popsize(self,ti):
        # Get the total number of people covered by this program
        # i.e. the sum of the source compartments of all links that
        # derive from this program
        # If impact_labels is specified, it must be a list of link labels
        # Then only links whose label is in that list will be included
        if ti == self.source_popsize_cache_time:
            return self.source_popsize_cache_val
        else:
            n = 0
            for link in self.links:
                n += link.source.vals[ti]
            self.source_popsize_cache_time = ti
            self.source_popsize_cache_val = n
            return n

class Link(Variable):
    '''
    A Link is a Variable that maps to a transition between compartments. As
    such, it contains an inflow and outflow compartment.  A Link's value is
    drawn from a Parameter, and there may be multiple links that draw values
    from the same parameter. The values stored in this extended version of
    Variable refer to flow rates. If used in ModelPop, the Link references two
    cascade compartments within a single population.
    '''
    def __init__(self, parameter, object_from, object_to, is_transfer=False):
        Variable.__init__(self, label=parameter.label) # A link should be labelled with the Parameter's label so it can be associated with that parameter later
        self.units = 'people/timestep'

        self.parameter = parameter # Source parameter where the unscaled link value is drawn from (a single parameter may have multiple links)
        self.parameter.dependency = True # A transition parameter must be updated during integration

        self.source = object_from # Compartment to remove people from
        self.dest = object_to # Compartment to add people to

        # Wire up references to this object
        self.parameter.links.append(self)
        self.source.outlinks.append(self)
        self.dest.inlinks.append(self)

        self.is_transfer = is_transfer # A transfer connections compartments across populations
        
        # Link vals stores the number of people actually transferred
        # The target flow also stores for each time point, the number of people proposed to move
        # The original parameter value is available from the Link's bound parameter
        self.target_flow = None # For each time point, store the number of people that were proposed to move (in units of number of people)
        
    def __repr__(self, *args, **kwargs):
        return "Link %s - %s to %s" % (self.label, self.source.label, self.dest.label)

    def preallocate(self,tvec):
        Variable.preallocate(self, tvec)
        self.target_flow = np.ones(tvec.shape) * np.nan

    def plot(self):
        Variable.plot(self)
        plt.title('Link %s to %s' % (self.source.label,self.dest.label))



# %% Cascade compartment and population classes
class ModelPopulation(Node):
    '''
    A class to wrap up data for one population within model.
    Each model population must contain a set of compartments with equivalent labels.
    '''

    def __init__(self, settings, label='default'):
        Node.__init__(self, label=label)
        self.comps = list()         # List of cascade compartments that this model population subdivides into.
        self.characs = list()       # List of output characteristics and parameters (dependencies computed during integration, pure outputs added after)
        self.links = list()         # List of intra-population cascade transitions within this model population.
        self.pars = list()
        
        self.comp_ids = dict()      # Maps label of a compartment to its position index within compartments list.
        self.charac_ids = dict()      # Maps cascade transition tag to indices for all relevant transitions within links list.
        self.par_ids = dict()      # Maps cascade transition tag to indices for all relevant transitions within links list.

        self.genCascade(settings=settings)    # Convert compartmental cascade into lists of compartment and link objects.

    def __repr__(self):
        return '%s "%s" (%s)' % (self.__class__.__name__,self.label,self.uid)

    def popsize(self,ti):
        # A population's popsize is the sum of all of the people in its compartments, excluding
        # birth and death compartments
        n = 0
        for comp in self.comps:
            if not comp.tag_birth and not comp.tag_dead:
                n += comp.vals[ti]
        return n

    def getModelState(self, ti):
        states = [c.getValue(ti) for c in self.comps]
        return states

    def getComp(self, comp_label):
        ''' Allow compartments to be retrieved by label rather than index. Returns a Compartment. '''
        comp_index = self.comp_ids[comp_label]
        return self.comps[comp_index]

    def getLinks(self, link_tag):
        ''' Retrieve Links associated with a Parameter tag '''
        par = self.getPar(link_tag)
        return par.links

    def getCharac(self, charac_label):
        ''' Allow dependencies to be retrieved by label rather than index. Returns a Variable. '''
        index = self.charac_ids[charac_label]
        return self.characs[index]

    def getPar(self, par_label):
        ''' Allow dependencies to be retrieved by label rather than index. Returns a Variable. '''
        index = self.par_ids[par_label]
        return self.pars[index]

    def genCascade(self, settings):
        '''
        Generate a compartmental cascade as defined in a settings object.
        Fill out the compartment, transition and dependency lists within the model population object.
        Maintaining order as defined in a cascade workbook is crucial due to cross-referencing.
        '''

        # First, make a Compartment for every compartment in the cascade
        tvec = settings.tvec

        for comp_id, label in enumerate(settings.node_specs.keys()):
            self.comps.append(Compartment(label=label))
            if 'tag_birth' in settings.node_specs[label]:
                self.comps[-1].tag_birth = True
            if 'tag_dead' in settings.node_specs[label]:
                self.comps[-1].tag_dead = True
            if 'junction' in settings.node_specs[label]:
                self.comps[-1].is_junction = True
            self.comp_ids[label] = comp_id

        # First pass, instantiate objects
        for charac_id,label in enumerate(settings.charac_specs.keys()):
            self.characs.append(Characteristic(label=label))
            self.charac_ids[label] = charac_id

        # Second pass, add includes and denominator
        for charac_id,label in enumerate(settings.charac_specs.keys()):
            charac = self.characs[charac_id]

            inc_labels = settings.charac_specs[label]['includes']
            den_label = settings.charac_specs[label]['denom'] if 'denom' in settings.charac_specs[label] else None

            includes = []
            for inc_label in inc_labels:
                if inc_label in self.comp_ids:
                    charac.add_include(self.comps[self.comp_ids[inc_label]])
                elif inc_label in self.charac_ids:
                    charac.add_include(self.characs[self.charac_ids[inc_label]])

            if den_label in self.comp_ids: 
                charac.add_denom(self.comps[self.comp_ids[den_label]])
            elif den_label in self.charac_ids:
                charac.add_denom(self.characs[self.charac_ids[den_label]])

        # Todo - Instantiation will work fine out of order now, but evaluation will not
        # Could safely reorder the list now to resolve dependencies and allow users to
        # define Characteristics in any order as long as a valid order of execution exists

        # Todo - perform similar adaptation for parameters
        # Next, create parameters
        for par_id,label in enumerate(settings.linkpar_specs.keys()):
            spec = settings.linkpar_specs[label]

            if 'f_stack' in spec:
                f_stack = dcp(spec['f_stack'])
                deps = []
                # Todo - clean up finding the dependent parameters
                for dep_label in spec['deps']:
                    if dep_label in self.comp_ids:
                        deps.append(self.comps[self.comp_ids[dep_label]])
                    elif dep_label in self.charac_ids:
                        deps.append(self.characs[self.charac_ids[dep_label]])
                    elif dep_label in self.par_ids:
                        deps.append(self.pars[self.par_ids[dep_label]])
                    else:
                        raise OptimaException('Could not find dependency %s for Parameter %s' % (dep_label,label))
            else:
                f_stack = None
                deps = None

            # Add limits
            if 'min' in spec or 'max' in spec:
                limits = [-np.inf,np.inf]
                if 'min' in spec:
                    limits[0] = spec['min']
                if 'max' in spec:
                    limits[1] = spec['max']
            else:
                limits = None

            self.pars.append(Parameter(label=label,deps=deps,f_stack=f_stack,limits=limits))
            self.par_ids[label] = par_id

        # Finally, create links between compartments
        for label,spec in settings.linkpar_specs.items():
            par = self.pars[self.par_ids[label]]
            if 'tag' in settings.linkpar_specs[label]:
                tag = settings.linkpar_specs[label]['tag']
                for pair in settings.links[tag]:
                    src = self.getComp(pair[0])
                    dst = self.getComp(pair[1])
                    new_link = Link(par,src,dst) # The link needs to be labelled with the Parameter it derives from so that Results can find it later
                    self.links.append(new_link)

    def preAllocate(self, sim_settings):
        '''
        Pre-allocate variable arrays in compartments, links and dependent variables for faster processing.
        Array maintains initial value but pre-fills everything else with NaNs.
        Thus errors due to incorrect parset value saturation should be obvious from results.
        '''
        tvec = sim_settings['tvec']
        for obj in self.comps + self.characs + self.links + self.pars:
            obj.preallocate(tvec)

# %% Model class
class Model(object):
    ''' A class to wrap up multiple populations within model and handle cross-population transitions. '''

    def __init__(self):

        self.pops = list()              # List of population groups that this model subdivides into.
        self.pop_ids = dict()           # Maps label of a population to its position index within populations list.
        self.contacts = dict()          # Maps interactions 'from' (i.e. a->b for [a][b]) and 'into' (i.e. a<-b for [a][b]) ModelPopulations, marking them with a weight.
        self.sim_settings = odict()
        self.t_index = 0                # Keeps track of array index for current timepoint data within all compartments.
        self.programs_active = None     # True or False depending on whether Programs will be used or not
        self.pset = None                # Instance of ModelProgramSet

    def unlink(self):
        # Break cycles when deepcopying or pickling by swapping them for UIDs
        # Primary storage is in the comps, links, and outputs properties
        for pop in self.pops:
             for obj in pop.comps + pop.characs + pop.pars + pop.links:
                 obj.unlink()
        if self.pset is not None:
            self.pset.unlink()
        self.pars_by_pop = None

    def relink(self):
        objs = {}
        self.pars_by_pop = dict()
        for pop in self.pops:
            for obj in pop.comps + pop.characs + pop.pars + pop.links:
                objs[obj.uid] = obj
                if isinstance(obj,Parameter):
                    if obj.label in self.pars_by_pop:
                        self.pars_by_pop[obj.label].append(obj)
                    else:
                        self.pars_by_pop[obj.label] = [obj]

        for pop in self.pops:
            for obj in pop.comps + pop.characs + pop.pars + pop.links:
                obj.relink(objs)

        if self.pset is not None:
            self.pset.relink(objs)

    def __getstate__(self):
        self.unlink()
        d = dcp(self.__dict__) # Deepcopy so that relinking doesn't modify the originals
        self.relink() # Relink, otherwise the original object gets unlinked
        return d

    def __setstate__(self, d):
        self.__dict__ = d
        self.relink()

    def getPop(self, pop_label):
        ''' Allow model populations to be retrieved by label rather than index. '''
        pop_index = self.pop_ids[pop_label]
        return self.pops[pop_index]

    def build(self, settings, parset, progset=None, options=None):
        ''' Build the full model. '''

        if options is None: options = dict()

        self.sim_settings['tvec'] = settings.tvec # NB. returning a mutable variable in a class @property method returns a new object each time
        self.sim_settings['tvec_dt'] = settings.tvec_dt

        self.sim_settings['impact_pars_not_func'] = []      # Program impact parameters that are not functions of other parameters and thus already marked for dynamic updating.
                                                            # This is only non-empty if a progset is being used in the model.
        parser.debug = settings.parser_debug

        for k, pop_label in enumerate(parset.pop_labels):
            self.pops.append(ModelPopulation(settings=settings, label=pop_label))
            self.pops[-1].preAllocate(self.sim_settings)     # Memory is allocated, speeding up model. However, values are NaN so as to enforce proper parset value saturation.
            self.pop_ids[pop_label] = k

        self.contacts = dcp(parset.contacts)    # Simple propagation of interaction details from parset to model.

        # Propagating initial characteristic parset values into ModelPops.
        # NOTE: Extremely involved process, so might be worth extracting the next few paragraphs as a separate method.

        t_init = np.array([self.sim_settings['tvec'][0]])

        seed_dict = _calculate_compartment_initialization(parset,settings,t_init)

        # Now initialise all model compartments with these calculated values.
        for seed_label in seed_dict:
            for pop_label in parset.pop_labels:
                val = seed_dict[seed_label][pop_label]
                if abs(val) < project_settings.TOLERANCE:
                    val = 0
                elif val < 0.:
                    raise OptimaException('ERROR: Initial value calculated for compartment "%s" in population "%s" is %f. Review and make sure each characteristic has at least as many people as the sum of all included compartments.' % (seed_label, pop_label, val))
                self.getPop(pop_label).getComp(seed_label).vals[0] = val

        # Propagating cascade parameter parset values into ModelPops. Handle both 'tagged' links and 'untagged' dependencies.
        for cascade_par in parset.pars['cascade']:
            for pop_label in parset.pop_labels:
                pop = self.getPop(pop_label)
                par = pop.pars[pop.par_ids[cascade_par.label]] # Find the parameter with the requested label

                par.vals = cascade_par.interpolate(tvec=self.sim_settings['tvec'], pop_label=pop_label)
                par.units = cascade_par.y_format[pop_label]
                par.scale_factor = cascade_par.y_factor[pop_label]

                # Now that we constrain everything, can leave it out here
                # But consider putting it back if the code is too slow
                # # Apply min/max restrictions on all parameters that are not functions.
                # # Functional parameters will be calculated and constrained during a run, hence they can be np.nan at this stage.
                # if not cascade_par.label in settings.par_funcs.keys():
                #     if 'min' in settings.linkpar_specs[cascade_par.label]:
                #         for pop_label in parset.pop_labels:
                #             for link_id in self.getPop(pop_label).link_ids[tag]:
                #                 vals = self.getPop(pop_label).links[link_id].vals
                #                 self.getPop(pop_label).links[link_id].vals[vals < settings.linkpar_specs[cascade_par.label]['min']] = settings.linkpar_specs[cascade_par.label]['min']
                #     if 'max' in settings.linkpar_specs[cascade_par.label]:
                #         for pop_label in parset.pop_labels:
                #             for link_id in self.getPop(pop_label).link_ids[tag]:
                #                 vals = self.getPop(pop_label).links[link_id].vals
                #                 self.getPop(pop_label).links[link_id].vals[vals > settings.linkpar_specs[cascade_par.label]['max']] = settings.linkpar_specs[cascade_par.label]['max']

        # Propagating transfer parameter parset values into Model object.
        # For each population pair, instantiate a Parameter with the values from the databook
        # For each compartment, instantiate a set of Links that all derive from that Parameter
        # NB. If a Program somehow the transfer parameter, those values will automatically
        # propagate to the links
        for trans_type in parset.transfers:
            if parset.transfers[trans_type]:
                for pop_source in parset.transfers[trans_type]:

                    transfer_parameter = parset.transfers[trans_type][pop_source] # This contains the data for all of the destrination pops

                    pop = self.getPop(pop_source)


                    for pop_target in transfer_parameter.y:

                        # Create the parameter object for this link (shared across all compartments)
                        par_label = trans_type + '_' + pop_source + '_to_' + pop_target # e.g. 'aging_0-4_to_15-64'
                        par = Parameter(label=par_label,f_stack = None, deps = None, limits = None)
                        par.preallocate(self.sim_settings['tvec'])
                        val = transfer_parameter.interpolate(tvec=self.sim_settings['tvec'], pop_label=pop_target)
                        par.vals = val
                        par.scale_factor = transfer_parameter.y_factor[pop_target]
                        par.units = transfer_parameter.y_format[pop_target]
                        par_id = len(pop.pars)
                        pop.pars.append(par)
                        pop.par_ids[par_label] = par_id

                        target_pop_obj = self.getPop(pop_target)

                        for source in pop.comps:
                            if not (source.tag_birth or source.tag_dead or source.is_junction):
                                # Instantiate a link between corresponding compartments
                                dest = target_pop_obj.getComp(source.label) # Get the corresponding compartment
                                link_tag = par_label + source.label # e.g. 'aging_0-4_to_15-64_sus' - this is probably never used?
                                link = Link(par, source, dest, is_transfer=True)
                                link.preallocate(self.sim_settings['tvec'])
                                link_id = len(pop.links)
                                pop.links.append(link)

        # # Make a lookup dict for programs
        self.pars_by_pop = {}
        for pop in self.pops:
            for par in pop.pars:
                if par.label in self.pars_by_pop:
                    self.pars_by_pop[par.label].append(par)
                else:
                    self.pars_by_pop[par.label] = [par]

        # Finally, prepare ModelProgramSet helper if programs are going to be used
        if 'progs_start' in options:
            if progset is not None:
                self.programs_active = True
                self.sim_settings['progs_start'] = options['progs_start']

                if 'progs_end' in options:
                    self.sim_settings['progs_end'] = options['progs_end']
                else:
                    self.sim_settings['progs_end'] = np.inf # Neverending programs
                if 'init_alloc' in options:
                    self.sim_settings['init_alloc'] = options['init_alloc']
                else: self.sim_settings['init_alloc'] = {}
                if 'constraints' in options:
                    self.sim_settings['constraints'] = options['constraints']
                if 'alloc_is_coverage' in options:
                    self.sim_settings['alloc_is_coverage'] = options['alloc_is_coverage']
                else: self.sim_settings['alloc_is_coverage'] = False
                if 'saturate_with_default_budgets' in options:
                    self.sim_settings['saturate_with_default_budgets'] = options['saturate_with_default_budgets']
                for impact_label in progset.impacts:
                    if impact_label not in settings.par_funcs:
                        self.sim_settings['impact_pars_not_func'].append(impact_label)

                self.pset = ModelProgramSet(progset,self.pops) # Make a ModelProgramSet wrapper
                self.pset.load_constraints(self.sim_settings['constraints'])
                alloc = self.pset.get_alloc(self.sim_settings)[0]
                self.pset.update_cache(alloc,self.sim_settings['tvec'],self.sim_settings['tvec_dt']) # Perform precomputations

            else:
                raise OptimaException('ERROR: A model run was initiated with instructions to activate programs, but no program set was passed to the model.')
        else:
            self.programs_active = False

        # Make sure initially-filled junctions are processed and initial dependencies are calculated.
        self.updateValues(settings=settings, do_special=False)     # Done first just in case junctions are dependent on characteristics.
                                                                                                # No special rules are applied at this stage, otherwise calculations would be iterated twice before the first step forward.
                                                                                                # NOTE: If junction outflows were to be tagged by special rules, initial calculations may be off. Return to this later and consider logic rigorously.
        self.processJunctions(settings=settings)
        self.updateValues(settings=settings)


        # set up sim_settings for later use wrt population tags
        self.sim_settings['tag_birth'] = []
        self.sim_settings['tag_death'] = []
        self.sim_settings['tag_no_plot'] = []
        for node_label in settings.node_specs:
            for tag in ['tag_birth', 'tag_death', 'tag_no_plot']:
                if tag in settings.node_specs[node_label]:
                    self.sim_settings[tag] = node_label


    def process(self, settings, progset=None):
        ''' 
        Run the full model.
        '''

        for t in self.sim_settings['tvec'][1:]:
            self.stepForward(settings=settings, dt=settings.tvec_dt)
            self.processJunctions(settings=settings)
            self.updateValues(settings=settings)

        for pop in self.pops:
            [par.update() for par in pop.pars if not par.dependency]
            [charac.update() for charac in pop.characs if not charac.dependency]

        return self.pops, self.sim_settings

    def stepForward(self, settings, dt=1.0):
        '''
        Evolve model characteristics by one timestep (defaulting as 1 year).
        Each application of this method writes calculated values to the next position in popsize arrays, regardless of dt.
        Thus the corresponding time vector associated with variable dt steps must be tracked externally.
        '''

        ti = self.t_index

        for pop in self.pops:
            for comp in pop.comps:
                comp.vals[ti+1] = comp.vals[ti]

        for pop in self.pops:

            for comp_source in pop.comps:

                if not comp_source.is_junction:  # Junctions collect inflows during this step. They do not process outflows here.

                    outlinks = comp_source.outlinks  # List of outgoing links
                    outflow = np.zeros(len(comp_source.outlinks))  # Outflow for each link # TODO - store in the link objects?

                    for i, link in enumerate(outlinks):

                        # Compute the number of people that are going out of each link
                        converted_amt = 0
                        transition = link.parameter.vals[ti]

                        if link.parameter.scale_factor is not None and link.parameter.scale_factor != project_settings.DO_NOT_SCALE:  # scale factor should be available to be used
                            transition *= link.parameter.scale_factor

                        if link.parameter.units == 'fraction':
                            # check if there are any violations, and if so, deal with them
                            if transition > 1.:
                                transition = checkTransitionFraction(transition, settings.validation)
                            converted_frac = 1 - (1 - transition) ** dt  # A formula for converting from yearly fraction values to the dt equivalent.
                            converted_amt = comp_source.vals[ti] * converted_frac
                        elif link.parameter.units == 'proportion':
                            converted_amt = comp_source.vals[ti] * transition # Use the value directly - NB. in theory this branch should only run for Junction links, which aren't here (i.e. it should never be used at the moment)
                        elif link.parameter.units == 'number':
                            converted_amt = transition * dt
                            if link.is_transfer:
                                transfer_rescale = comp_source.vals[ti] / pop.getCharac(settings.charac_pop_count).vals[ti]
                                converted_amt *= transfer_rescale
                        else:
                            raise OptimaException('Unknown parameter units!')

                        outflow[i] = converted_amt
                        link.target_flow[ti] = converted_amt

                    # Prevent negative population by proportionately downscaling the outflow
                    # if there are insufficient people _currently_ in the compartment
                    # Rescaling is performed if the validation setting is 'avert', otherwise
                    # either a warning will be displayed or an error will be printed
                    if not comp_source.tag_birth and np.sum(outflow) > comp_source.vals[ti]:
                        validation_level = settings.validation['negative_population']

                        if validation_level == project_settings.VALIDATION_AVERT or validation_level == project_settings.VALIDATION_WARN:
                            outflow = outflow / np.sum(outflow) * comp_source.vals[ti]
                        else:
                            warning = "Negative value encountered for: (%s - %s) at ti=%g : popsize = %g, outflow = %g" % (pop.label,comp_source.label,ti,comp_source.vals[ti],sum(outflow))
                            if validation_level == project_settings.VALIDATION_ERROR:
                                raise OptimaException(warning)
                            elif validation_level == project_settings.VALIDATION_WARN:
                                logger.warn(warning)

                    # Apply the flows to the compartments
                    for i, link in enumerate(outlinks):
                        link.dest.vals[ti+1] += outflow[i]
                        link.vals[ti] = outflow[i]

                    comp_source.vals[ti+1] -= np.sum(outflow)


        # Guard against populations becoming negative due to numerical artifacts
        for pop in self.pops:
            for comp in pop.comps:
                comp.vals[ti+1] = max(0,comp.vals[ti+1])

        # Update timestep index.
        self.t_index += 1

    def processJunctions(self, settings):
        '''
        For every compartment considered a junction, propagate the contents onwards until all junctions are empty.
        '''

        ti = self.t_index
        ti_link = ti - 1
        if ti_link < 0: 
            ti_link = ti    # For the case where junctions are processed immediately after model initialisation.

        review_required = True
        review_count = 0
        while review_required:
            review_count += 1
            review_required = False # Don't re-run unless a junction has refilled

            if review_count > settings.recursion_limit:
                raise OptimaException('ERROR: Processing junctions (i.e. propagating contents onwards) for timestep %i is taking far too long. Infinite loop suspected.' % ti_link)

            for pop in self.pops:
                junctions = [comp for comp in pop.comps if comp.is_junction]

                for junc in junctions:

                    # Stores junction popsize values before emptying - so vals_old stores the total number of people that transitioned out of this junction while vals = 0
                    if review_count == 1:
                        junc.vals_old[ti] = junc.vals[ti] # Back up the old value (should have been preallocated)
                    else:
                        # If a junction is being reviewed again, it means that it received inflow before emptying.
                        # Add this inflow to the stored popsize.
                        junc.vals_old[ti] += junc.vals[ti]

                    # If the compartment is numerically empty, make it empty
                    if junc.vals[ti] <= project_settings.TOLERANCE:   # Includes negative values.
                        junc.vals[ti] = 0
                    else:
                        current_size = junc.vals[ti]
                        denom_val = sum(link.parameter.vals[ti_link] for link in junc.outlinks) # This is the total number of people in the outflow compartments, at the previous timestep, used for splitting the outputs
                        if denom_val == 0:
                            raise OptimaException('ERROR: Proportions for junction "%s" outflows sum to zero, resulting in a nonsensical ratio. There may even be (invalidly) no outgoing transitions for this junction.' % junction_label)
                        for link in junc.outlinks:
                            if review_count == 1:
                                link.vals[ti] = 0
                                link.target_flow[ti] = 0
                            flow = current_size * link.parameter.vals[ti_link] / denom_val
                            link.source.vals[ti] -= flow
                            link.dest.vals[ti]   += flow
                            link.vals[ti] += flow
                            link.target_flow[ti] += flow
                            if link.dest.is_junction:
                                review_required = True # Need to review if a junction received an inflow at this step

    def updateValues(self, settings, do_special=True):
        '''
        Run through all parameters and characteristics flagged as dependencies for custom-function parameters and evaluate them for the current timestep.
        These dependencies must be calculated in the same order as defined in settings, characteristics before parameters, otherwise references may break.
        Also, parameters in the dependency list do not need to be calculated unless explicitly depending on another parameter.
        Parameters that have special rules are usually dependent on other population values, so are included here.
        '''

        ti = self.t_index

        # The output list maintains the order in which characteristics and parameters appear in the
        # settings, and also puts characteristics before parameters. So iterating through that list
        # will automatically compute them in the correct order

        # First, compute dependent characteristics, as parameters might depend on them
        for pop in self.pops:
            for charac in pop.characs:
                if charac.dependency:
                    charac.update(ti)
        
        # 1st:  Update parameters if they have an f_stack variable
        # 2nd:  Any parameter that is overwritten by a program-based cost-coverage-impact transformation.
        # 3rd:  Any parameter that is overwritten by a special rule, e.g. averaging across populations.
        # 4th:  Any parameter that is restricted within a range of values, i.e. by min/max values.
        # Looping through populations must be internal so that all values are calculated before special inter-population rules are applied.
        # We resolve one parameter at a time, in dependency order
        do_program_overwrite = self.programs_active and self.sim_settings['tvec'][ti] >= self.sim_settings['progs_start'] and self.sim_settings['tvec'][ti] <= self.sim_settings['progs_end']
        if do_program_overwrite:
            prog_vals = self.pset.compute_pars(ti)[0]

        for par_label in (settings.par_funcs.keys() + self.sim_settings['impact_pars_not_func']):
            pars = self.pars_by_pop[par_label] # All of the parameters with this label, across populations. There should be one for each population (these are Parameters, not Links)

            # First - update parameters that are dependencies, evaluating f_stack if required
            for par in pars:
                if par.dependency:
                    par.update(ti)

            # Then overwrite with program values
            if do_program_overwrite:
                for par in pars:
                    if par.uid in prog_vals:
                        par.vals[ti] = prog_vals[par.uid]

            # Handle parameters tagged with special rules. Overwrite vals if necessary.
            if do_special and 'rules' in settings.linkpar_specs[par_label]:
                pars = self.pars_by_pop[par_label]  # All of the parameters with this label, across populations. There should be one for each population (these are Parameters, not Links)

                for par in pars:
                    par.vals_old[ti] = par.vals[ti]

                rule = settings.linkpar_specs[par_label]['rules']
                for pop in self.pops:
                    if rule == 'avg_contacts_in':
                        from_list = self.contacts['into'][pop.label].keys()

                        # If interactions with a pop are initiated by the same pop, no need to proceed with special calculations. Else, carry on.
                        if not ((len(from_list) == 1 and from_list[0] == pop.label)):
                            old_vals = np.ones(len(from_list)) * np.nan
                            weights = np.ones(len(from_list)) * np.nan
                            pop_counts = np.ones(len(from_list)) * np.nan
                            if len(from_list) == 0:
                                new_val = 0.0
                            else:
                                k = 0
                                for from_pop in from_list:
                                    # All transition links with the same par_label are identically valued. For calculations, only one is needed for reference.
                                    par = self.getPop(from_pop).getPar(par_label)
                                    old_vals[k] = par.vals_old[ti]
                                    weights[k] = self.contacts['into'][pop.label][from_pop]
                                    pop_counts[k] = self.getPop(from_pop).getCharac(settings.charac_pop_count).vals[ti]
                                    k += 1
                                wpc = np.multiply(weights, pop_counts)          # Population counts weighted by contact rate.
                                wpc_sum = sum(wpc)                              # Normalisation factor for weighted population counts.
                                if abs(wpc_sum) > project_settings.TOLERANCE:
                                    new_val = np.dot(old_vals, wpc / wpc_sum)         # Do a weighted average of the parameter values pertaining to contact-initiating pop groups.
                                else:
                                    new_val = 0.0   # Only valid because if the weighted sum is zero, all pop_counts must be zero, meaning that the numerator is zero.

                            # Update the parameter's value in this population - will propagate to links in next stage
                            pop.getPar(par_label).vals[ti] = new_val

            # Restrict the parameter's value if a limiting range was defined
            for par in pars:
                par.constrain(ti)

    def calculateOutputs(self, settings):
        '''
        Calculate outputs (called cascade characteristics in settings).
        These outputs must be calculated in the same order as defined in settings, otherwise references may break.
        Include any parameters marked as an output in the cascade sheet.
        Return a dictionary with all of the outputs

        As outputs are all Variables or Characteristics they use the vals property for their value
        '''

        outputs = odict()

        for pop in self.pops:
            for output in pop.pars + pop.characs + pop.links:
                if isinstance(output,Parameter) and len(output.links) > 0:
                    continue # Load in links instead of this parameter
                if output.label not in outputs:
                    outputs[output.label] = odict()
                if pop.label in outputs[output.label]:
                    assert isinstance(output,Link), 'There is a duplicate output label that is NOT a link - this is not supposed to happen'
                    outputs[output.label][pop.label] += output.vals
                else:
                    outputs[output.label][pop.label] = output.vals
                
        return outputs

    def printModelState(self, ti):
        for pop in self.pops:
            logging.info("Population: %s" % pop.label)
            logging.info(pop.getModelState(ti))


def runModel(settings, parset, progset=None, options=None):
    '''
    Processes the TB epidemiological model.
    Parset-based overwrites are generally done externally, so the parset is only used for model-building.
    Progset-based overwrites take place internally and must be part of the processing step.
    The options dictionary is usually passed in with progset to specify when the overwrites take place.
    '''

    m = Model()
    m.build(settings=settings, parset=parset, progset=progset, options=options)
    m.process(settings=settings, progset=progset)

    results = ResultSet(m, parset, settings, progset, options)    # NOTE: Progset may need to be passed to results. Depends on what results object stores.

    return results


def _calculate_compartment_initialization(parset,settings,t_init):
    seed_dict = {}              # Compartment values to initialise with.
    include_dict = odict()      # Lowest-level nodes included for each entry-point characteristic, keyed by entry-point.
                                # They are often inserted in definitional order, so useful to keep as an odict. Speeds up calculation process.
    calc_done = {}

    # All compartments are either characteristic entry-points or contain zero people.
    # First, generate a dictionary of prospective values to seed compartments with, all assumed to be zero.
    # Assume the values for all compartments have been calculated.
    for node_label in settings.node_specs:
        seed_dict[node_label] = odict()
        calc_done[node_label] = True    # Values are technically already settled for nodes if a node is not an entry-point.
                                        # Alternatively, the node can be an entry-point of a characteristic including only the entry-point.
        for pop_label in parset.pop_labels:
            seed_dict[node_label][pop_label] = 0.

    # Now update assumptions by looping through all characteristics containing entry-points.
    # We note that initial values for entry-points are derived from the value of the characteristic minus the values of all other 'included' compartments.
    # For inclusions of only the entry-point and nothing else, its seeding value is simply updated with the interpolated value of the characteristic.
    # For inclusions of more compartments (once flattened out), calculations are more difficult.
    # The characteristic seeding value is still updated, but the entry-point must be removed from the dictionary that tracks calculated compartments.
    for charac_label in settings.charac_specs:
        if 'entry_point' in settings.charac_specs[charac_label]:
            ep_label = settings.charac_specs[charac_label]['entry_point']
            flat_list, dep_list = flattenDict(input_dict=settings.charac_specs, base_key=charac_label, sub_keys=['includes'])
            flat_list.remove(ep_label)
            if len(flat_list) > 0:
                del calc_done[ep_label]
                include_dict[ep_label] = dcp(flat_list)
            par = parset.pars['characs'][parset.par_ids['characs'][charac_label]]
            for pop_label in parset.pop_labels:
                val = par.interpolate(tvec=t_init, pop_label=pop_label)[0]
                seed_dict[ep_label][pop_label] = val

    # Loop through and handle prevalences (i.e. characteristics with denominators).
    # The denominator value will be drawn from the parset, not the seed dictionary, so beware a prevalence as a denominator.
    for charac_label in settings.charac_specs:
        if 'entry_point' in settings.charac_specs[charac_label]:
            if 'denom' in settings.charac_specs[charac_label]:
                ep_label = settings.charac_specs[charac_label]['entry_point']
                denom_label = settings.charac_specs[charac_label]['denom']
                par = parset.pars['characs'][parset.par_ids['characs'][denom_label]]
                for pop_label in parset.pop_labels:
                    val = par.interpolate(tvec=t_init, pop_label=pop_label)[0]
                    seed_dict[ep_label][pop_label] *= val

    # Now loop through all 'uncalculated entry-points'.
    # If any of their remaining inclusions have been calculated in previous loops, stop tracking the included compartment and subtract its value from the entry-point.
    # Eventually, an entry-point should be equal to its characteristic seeding value minus the correction of all other included compartments.
    # With no more inclusions to keep track of, this entry-point is considered fully calculated and can be subtracted from other 'uncalculated entry-points'.
    # Eventually there will be no inclusions left for any entry-point, meaning all initial values are calculated.
    review_count = 0
    while len(include_dict) > 0:
        for entry_point in dcp(include_dict.keys()):
            for include in dcp(include_dict[entry_point]):

                # Subtract the values of any included already-calculated nodes from the value of an entry-point.
                if include in calc_done:
                    for pop_label in parset.pop_labels:
                        val = seed_dict[entry_point][pop_label] - seed_dict[include][pop_label]
                        if val < 0 and abs(val) > project_settings.TOLERANCE:
                            logger.error('Negative value encountered for Entry point: %s, Pop_label: %s, Compartment: %s    Entry point size: %f, compartment size: %f' % (entry_point, pop_label, include, seed_dict[entry_point][pop_label], seed_dict[include][pop_label]))
                        seed_dict[entry_point][pop_label] -= seed_dict[include][pop_label]
                    include_dict[entry_point].remove(include)

                # If all included nodes have been calculated and subtracted from an entry-point, then the entry-point is now a fully calculated node.
                # It can be used in subtractions for other nodes.
                if len(include_dict[entry_point]) == 0:
                    calc_done[entry_point] = True
                    del include_dict[entry_point]
        review_count += 1
        if review_count > len(settings.node_specs):
            raise OptimaException('ERROR: Calculation phase for initial compartment values has looped more times than the number of compartments. Something is likely wrong with characteristic definitions.')

    return seed_dict