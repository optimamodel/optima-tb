# %% Imports

from optima_tb.utils import flattenDict, odict, OptimaException
import optima_tb.settings as project_settings
from optima_tb.results import ResultSet
from optima_tb.parsing import parse_function
from optima_tb.ModelPrograms import ModelProgramSet, ModelProgram
from collections import defaultdict
import cPickle as pickle
import logging
logger = logging.getLogger(__name__)

import numpy as np
import scipy.optimize
from copy import deepcopy as dcp
import uuid

import matplotlib.pyplot as plt



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
        self.dt = None
        if 'vals' not in dir(self): # characteristics already have a vals method
            self.vals = None
        self.units = 'unknown'              # 'unknown' units are distinct to dimensionless units, that have value ''

    def preallocate(self,tvec,dt):
        self.t = tvec
        self.dt = dt
        self.vals = np.ones(tvec.shape) * np.nan

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

    def set_dependent(self):
        # Make the variable a dependency. For Compartments and Links, this does nothing. For Characteristics and
        # Parameters, it will set the dependent flag, but in addition, any validation constraints e.g. a Parameter
        # that depends on Links cannot itself be a dependency, will be enforced
        return

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
        self.vals = None

        self.outlinks = []
        self.inlinks = []

    def unlink(self):
        self.outlinks = [x.uid for x in self.outlinks]
        self.inlinks = [x.uid for x in self.inlinks]

    def relink(self,objs):
        self.outlinks = [objs[x] for x in self.outlinks]
        self.inlinks = [objs[x] for x in self.inlinks]

    @property
    def outflow(self):
        # Return the outflow at each timestep - for a junction, this is equal to the number
        # of people that were in the junction
        x = np.zeros(self.t.shape)
        if self.outlinks:
            for link in self.outlinks:
                x += link.vals
        return x

    def expected_duration(self,ti=None):
        # Returns the expected number of years that an individual is expected to remain
        # in this compartment for, if the outgoing flow rates are maintained
        if ti is None:
            ti = np.arange(0,len(self.t))

        if self.tag_birth or self.tag_dead or self.is_junction:
            return 0

        outflow_probability = 0
        for link in self.outlinks:
            if link.parameter.units == 'fraction':
                outflow_probability += 1 - (1 - link.parameter.vals[ti]) ** self.dt  # A formula for converting from yearly fraction values to the dt equivalent.
            elif link.parameter.units == 'number':
                outflow_probability += link.parameter.vals[ti]*self.dt/self.vals[ti]
            else:
                raise OptimaException('Unknown parameter units')

        remain_probability = 1-outflow_probability
        dur = np.zeros(outflow_probability.shape)
        # From OSL/HMM-MAR:
        # Given a transition probability, what is the expected lifetime in units of steps?
        # This can be determined using a symbolic integration, below
        # syms k_step p;
        # assume(p>0);
        # f =(k_step)*p^k_step*(1-p); # (probability that the state lasts k *more* steps, multiplied by lifetime which is k)
        # fa = 1+int(f,k_step,0,inf); # Add 1 to the lifetime because all states last at least 1 sample
        # f = @(x) double(subs(fa,p,x));
        #
        # However, the result of the symbolic integration contains a limit which is
        # zero unless p=1, but p=1 is not needed because we know it means the compartment immediately empties.
        # So instead of the function above, can instead drop the limit term and write the rest of the
        # expression out which gives identical results from p=0.00001 to p=0.99 (note that if dirichletdiag>=1)
        # then the upper bound on the transition probability is p=0.5 anyway for K=2
        dur[dur<1] = (1 - (1. / np.log(remain_probability[dur<1])**2)*(remain_probability[dur<1] - 1))*self.dt
        return dur

    def expected_outflow(self,ti):
        # After 1 year, where are people expected to be? If people would leave in less than a year,
        # then the numbers correspond to when the compartment is empty
        popsize = self.vals[ti] # This is the number of people who are leaving

        outflow = {link.dest.label: popsize * link.vals[ti]/self.dt for link in self.outlinks}
        rescale = popsize/sum([y for _,y in outflow.items()])
        for x in outflow.keys():
            outflow[x] *= rescale
        outflow[self.label] = popsize - sum([y for _,y in outflow.items()])

        return outflow

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
        self.internal_vals = None

    def preallocate(self,tvec,dt):
        self.t = tvec
        self.dt = dt
        self.internal_vals = np.ones(tvec.shape) * np.nan

    def get_included_comps(self):
        includes = []
        for inc in self.includes:
            if isinstance(inc,Characteristic):
                includes += inc.get_included_comps()
            else:
                includes.append(inc)
        return includes

    @property
    def vals(self):
        if self.internal_vals is None:
            vals = np.zeros(self.t.shape)

            for comp in self.includes:
                vals += comp.vals

            if self.denominator is not None:
                denom = self.denominator.vals
                vals_zero = vals < project_settings.TOLERANCE
                vals[denom > 0] /= denom[denom > 0]
                vals[vals_zero] = 0.0
                vals[(denom<=0) & (~vals_zero)] = np.inf

            return vals
        else:
            return self.internal_vals


    def set_dependent(self):
        self.dependency = True

        for inc in self.includes:
            inc.set_dependent()

        if self.denominator is not None:
            self.denominator.set_dependent()

    def unlink(self):
        self.includes = [x.uid for x in self.includes]
        self.denominator = self.denominator.uid if self.denominator is not None else None

    def relink(self,objs):
        # Given a dictionary of objects, restore the internal references
        # based on the UUID
        self.includes = [objs[x] for x in self.includes]
        self.denominator = objs[self.denominator] if self.denominator is not None else None

    def add_include(self,x):
        assert isinstance(x,Compartment) or isinstance(x,Characteristic)
        self.includes.append(x)

    def add_denom(self,x):
        assert isinstance(x,Compartment) or isinstance(x,Characteristic)
        self.denominator = x
        self.units = ''

    def update(self,ti):
        self.internal_vals[ti] = 0
        for comp in self.includes:
            self.internal_vals[ti] += comp.vals[ti]
        if self.denominator is not None:
            denom  = self.denominator.vals[ti]
            if denom > 0:
                self.internal_vals[ti] /= denom
            elif self.internal_vals[ti] < project_settings.TOLERANCE:
                self.internal_vals[ti] = 0  # Given a zero/zero case, make the answer zero.
            else:
                self.internal_vals[ti] = np.inf  # Given a non-zero/zero case, keep the answer infinite.

class Parameter(Variable):
    # A parameter is a Variable that can have a value computed via an f_stack and a list of 
    # dependent Variables. This class may need to be renamed to avoid confusion with
    # the class in Parameter.py which provides a means of computing the Parameter that is used by the model.
    # This is a Parameter in the cascade.xlsx sense - there is one Parameter object for every item in the 
    # Parameters sheet. A parameter that maps to multiple transitions (e.g. doth_rate) will have one parameter
    # and multiple Link instances that depend on the same Parameter instance
    #
    #  *** Parameter values are always annualized ***
    def __init__(self, label='default'):
        Variable.__init__(self, label=label)
        self.vals = None
        self.deps = None
        self.fcn_str = None
        self._fcn = None

        self.limits = None # Can be a two element vector [min,max]

        self.dependency = False 
        self.scale_factor = 1.0
        self.links = [] # References to links that derive from this parameter
        self.source_popsize_cache_time = None
        self.source_popsize_cache_val = None

    def set_fcn(self, fcn_str, pop):
        self.fcn_str = fcn_str
        self._fcn, dep_list = parse_function(fcn_str)
        deps = {}
        for dep_name in dep_list:
            deps[dep_name] = pop.getVariable(dep_name)
        self.deps = deps

        # If this Parameter has links, it must be marked as dependent for evaluation during integration
        if self.links:
            self.set_dependent()

    def set_dependent(self):
        self.dependency = True
        if self.deps is not None: # Make all dependencies dependent too, this will propagate through dependent parameters
            for deps in self.deps.values():
                for dep in deps:
                    if isinstance(dep,Link):
                        raise OptimaException('A Parameter that depends on transition flow rates cannot be a dependency, it must be output only')
                    dep.set_dependent()

    def unlink(self):
        self.links = [x.uid for x in self.links]
        if self.deps is not None:
            for dep_name in self.deps:
                self.deps[dep_name] = [x.uid for x in self.deps[dep_name]]
        if self._fcn is not None:
            self._fcn = None

    def relink(self,objs):
        # Given a dictionary of objects, restore the internal references
        # based on the UUID
        self.links = [objs[x] for x in self.links]

        if self.deps is not None:
            for dep_name in self.deps:
                self.deps[dep_name] = [objs[x] for x in self.deps[dep_name]]

        if self.fcn_str:
            self._fcn, dep_list = parse_function(self.fcn_str)


    def constrain(self,ti):
        # NB. Must be an array, so ti must must not be supplied
        if self.limits is not None:
            self.vals[ti] = max(self.limits[0],self.vals[ti])
            self.vals[ti] = min(self.limits[1],self.vals[ti])

    def update(self,ti=None):
        # Update the value of this Parameter at time index ti
        # by evaluating its f_stack function using the 
        # current values of all dependent variables at time index ti
        if not self._fcn:
            return

        if ti is None:
            ti = np.arange(0,self.vals.size) # This corresponds to every time point
        else:
            ti = np.array(ti)

        dep_vals = dict.fromkeys(self.deps,0.0)
        for dep_name,deps in self.deps.items():
            for dep in deps:
                if isinstance(dep, Link):
                    dep_vals[dep_name] += dep.vals[[ti]] / dep.dt
                else:
                    dep_vals[dep_name] += dep.vals[[ti]]

        self.vals[ti] = self.scale_factor*self._fcn(**dep_vals)

    def source_popsize(self,ti):
        # Get the total number of people covered by this program
        # i.e. the sum of the source compartments of all links that
        # derive from this program
        # If impact_labels is specified, it must be a list of link labels
        # Then only links whose label is in that list will be included
        if ti == self.source_popsize_cache_time:
            return self.source_popsize_cache_val
        else:
            n = 0.0
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
    # *** Link values are always dt-based ***
    def __init__(self, parameter, object_from, object_to,tag):
        # Note that the Link's label is the transition tag
        Variable.__init__(self, label=tag)
        self.vals = None
        self.tag = tag
        self.units = 'people'

        self.parameter = parameter # Source parameter where the unscaled link value is drawn from (a single parameter may have multiple links)

        self.source = object_from # Compartment to remove people from
        self.dest = object_to # Compartment to add people to

        # Wire up references to this object
        self.parameter.links.append(self)
        self.source.outlinks.append(self)
        self.dest.inlinks.append(self)

    def unlink(self):
        self.parameter = self.parameter.uid
        self.source = self.source.uid
        self.dest = self.dest.uid

    def relink(self,objs):
        # Given a dictionary of objects, restore the internal references
        # based on the UUID
        self.parameter = objs[self.parameter]
        self.source = objs[self.source]
        self.dest = objs[self.dest]

    def __repr__(self, *args, **kwargs):
        return "Link %s (parameter %s) - %s to %s" % (self.label, self.parameter.label, self.source.label, self.dest.label)

    def plot(self):
        Variable.plot(self)
        plt.title('Link %s to %s' % (self.source.label,self.dest.label))

# %% Cascade compartment and population classes
class ModelPopulation(object):
    '''
    A class to wrap up data for one population within model.
    Each model population must contain a set of compartments with equivalent labels.
    '''

    def __init__(self, settings, label='default'):
        self.uid = uuid.uuid4()
        self.label = label          # Reference name for this object.

        self.comps = list()         # List of cascade compartments that this model population subdivides into.
        self.characs = list()       # List of output characteristics and parameters (dependencies computed during integration, pure outputs added after)
        self.links = list()         # List of intra-population cascade transitions within this model population.
        self.pars = list()
        
        self.comp_lookup = dict()      # Maps label of a compartment to a compartment
        self.charac_lookup = dict()
        self.par_lookup= dict()
        self.link_lookup = dict() # Map label of Link to a list of Links with that label

        self.genCascade(settings=settings)    # Convert compartmental cascade into lists of compartment and link objects.

        self.popsize_cache_time = None
        self.popsize_cache_val = None

    def __repr__(self):
        return '%s "%s" (%s)' % (self.__class__.__name__,self.label,self.uid)

    def unlink(self):
        for obj in self.comps + self.characs + self.pars + self.links:
            obj.unlink()
        self.comp_lookup = None
        self.charac_lookup = None
        self.par_lookup = None
        self.link_lookup = None

    def relink(self,objs):
        for obj in self.comps + self.characs + self.pars + self.links:
            obj.relink(objs)
        self.comp_lookup = {comp.label:comp for comp in self.comps}
        self.charac_lookup = {charac.label:charac for charac in self.characs}
        self.par_lookup = {par.label:par for par in self.pars}
        link_labels = set([link.label for link in self.links])
        self.link_lookup = {label:[link for link in self.links if link.label == label] for label in link_labels}

    def popsize(self,ti=None):
        # A population's popsize is the sum of all of the people in its compartments, excluding
        # birth and death compartments
        if ti is None:
            return np.sum([comp.vals for comp in self.comps if (not comp.tag_birth and not comp.tag_dead)],axis=0)

        if ti == self.popsize_cache_time:
            return self.popsize_cache_val
        else:
            n = 0
            for comp in self.comps:
                if not comp.tag_birth and not comp.tag_dead:
                    n += comp.vals[ti]

            return n

    def getVariable(self,label):
        # Returns a list of variables whose label matches the requested label
        # At the moment, labels are unique across object types and within object
        # types except for links, but if that logic changes, simple modifications can
        # be made here

        label = label.replace('___',':')

        if label in self.comp_lookup:
            return [self.comp_lookup[label]]
        elif label in self.charac_lookup:
            return [self.charac_lookup[label]]
        elif label in self.par_lookup:
            return [self.par_lookup[label]]
        elif label in self.link_lookup:
            return self.link_lookup[label]
        elif ':' in label:
            # Link labels in Atomica end in 'par_name:flow' so they are handled above
            # This branch of the if statement is exclusively for compartment lookup
            # Allowed syntax is
            # 'source_label:' - All links going out from source
            # ':dest_label' - All links going into destination
            # 'source_label:dest_label' - All links from Source to Dest
            # 'source_label:dest_label:par_label' - All links from Source to Dest belonging to a given Parameter
            # ':dest:par_label'
            # 'source::par_label' - As per above
            # Note - because compartment labels are resolved within compartments, this function currently cannot
            # be used to look up transfer links
            label_tokens = label.split(':')
            if len(label_tokens) == 2:
                label_tokens.append('')
            src,dest,par = label_tokens

            if src and dest:
                links = [l for l in self.getComp(src).outlinks if l.dest.label == dest]
            elif src:
                links = self.getComp(src).outlinks
            elif dest:
                links = self.getComp(dest).inlinks
            else:
                links = self.links

            if par:
                links = [l for l in links if l.parameter.label == par]

            return links

        else:
            raise OptimaException('Object not found: (%s)' % (label))

    def getComp(self, comp_label):
        ''' Allow compartments to be retrieved by label rather than index. Returns a Compartment. '''
        return self.comp_lookup[comp_label]

    def getLinks(self, label):
        ''' Retrieve Links'''
        # Links can be looked up by parameter label or by link label, unlike getVariable. This is because
        # getLinks() is guaranteed to return a list of Link objects
        # As opposed to getVariable which would retrieve the Parameter for 'doth rate' and the Links for 'z'
        if label in self.par_lookup:
            return self.par_lookup[label].links
        elif label in self.link_lookup:
            return self.link_lookup[label]
        else:
            raise OptimaException('Object %s not found' % (label))

    def getCharac(self, charac_label):
        ''' Allow dependencies to be retrieved by label rather than index. Returns a Variable. '''
        return self.charac_lookup[charac_label]

    def getPar(self, par_label):
        ''' Allow dependencies to be retrieved by label rather than index. Returns a Variable. '''
        return self.par_lookup[par_label]

    def genCascade(self, settings):
        '''
        Generate a compartmental cascade as defined in a settings object.
        Fill out the compartment, transition and dependency lists within the model population object.
        Maintaining order as defined in a cascade workbook is crucial due to cross-referencing.
        '''

        # Instantiate Compartments
        for label,spec in settings.node_specs.items():
            self.comps.append(Compartment(label=label))
            if 'tag_birth' in spec:
                self.comps[-1].tag_birth = True
            if 'tag_dead' in spec:
                self.comps[-1].tag_dead = True
            if 'junction' in spec:
                self.comps[-1].is_junction = True
        self.comp_lookup = {comp.label:comp for comp in self.comps}

        # Characteristics first pass, instantiate objects
        for label in settings.charac_specs:
            self.characs.append(Characteristic(label=label))
        self.charac_lookup = {charac.label:charac for charac in self.characs}

        # Characteristics second pass, add includes and denominator
        for label,spec in settings.charac_specs.items():
            charac = self.getCharac(label)
            for inc_label in spec['includes']:
                charac.add_include(self.getVariable(inc_label)[0]) # nb. We expect to only get one match for the label, so use index 0
            if 'denom' in spec:
                charac.add_denom(self.getVariable(spec['denom'])[0])

        # Finally, make sure the 'alive' characteristic is kept up to date
        self.getCharac(settings.charac_pop_count).set_dependent()


        # Parameters first pass, create parameter objects and links
        for label,spec in settings.linkpar_specs.items():
            par = Parameter(label=label)
            self.pars.append(par)
            if 'tag' in settings.linkpar_specs[label]:
                tag = settings.linkpar_specs[label]['tag']
                for pair in settings.links[tag]:
                    src = self.getComp(pair[0])
                    dst = self.getComp(pair[1])
                    new_link = Link(par,src,dst,tag) # The link needs to be labelled with the Parameter it derives from so that Results can find it later
                    if tag not in self.link_lookup:
                        self.link_lookup[tag] = [new_link]
                    else:
                        self.link_lookup[tag].append(new_link)
                    self.links.append(new_link)
        self.par_lookup = {par.label:par for par in self.pars}

        # Parameters second pass, process f_stacks, deps, and limits
        for label,spec in settings.linkpar_specs.items():
            par = self.getPar(label)

            if ('min' in spec) or ('max' in spec):
                par.limits = [-np.inf, np.inf]
                if 'min' in spec:
                    par.limits[0] = spec['min']
                if 'max' in spec:
                    par.limits[1] = spec['max']

            if 'fcn_str' in spec:
                par.set_fcn(spec['fcn_str'],self)

    def preallocate(self, tvec, dt):
        '''
        Pre-allocate variable arrays in compartments, links and dependent variables for faster processing.
        Array maintains initial value but pre-fills everything else with NaNs.
        Thus errors due to incorrect parset value saturation should be obvious from results.
        '''
        for obj in self.comps + self.characs + self.links + self.pars:
            obj.preallocate(tvec,dt)

    def initialize_compartments(self,parset,settings,t_init):
        # Given a set of characteristics and their initial values, compute the initial
        # values for the compartments by solving the set of characteristics simultaneously

        characs = [c for c in self.characs if settings.charac_specs[c.label]['databook_order']!=-1]
        comps = [c for c in self.comps if not (c.tag_birth or c.tag_dead)]
        charac_indices = {c.label:i for i,c in enumerate(characs)} # Make lookup dict for characteristic indices
        comp_indices = {c.label:i for i,c in enumerate(comps)} # Make lookup dict for compartment indices

        b = np.zeros((len(characs),1))
        A = np.zeros((len(characs),len(comps)))

        # Construct the characteristic value vector (b) and the includes matrix (A)
        for i,c in enumerate(characs):
            # Look up the characteristic value
            par = parset.pars['characs'][parset.par_ids['characs'][c.label]]
            b[i] = par.interpolate(tvec=np.array([t_init]), pop_label=self.label)[0]
            if c.denominator is not None:
                denom_par = parset.pars['characs'][parset.par_ids['characs'][c.denominator.label]]
                b[i] *= denom_par.interpolate(tvec=np.array([t_init]), pop_label=self.label)[0]
            for inc in c.get_included_comps():
                A[i,comp_indices[inc.label]] = 1.0

        # Solve the linear system - allowing negative popsizes to be computed
        x, residual, rank, _ = np.linalg.lstsq(A,b,rcond=-1)

        # Solve the linear system - do not allow negative popsizes
        # x, residual = scipy.optimize.nnls(A,b.ravel())
        # rank = np.linalg.matrix_rank(A)

        # Halt if the solution is not unique (could relax this check later)
        if rank < A.shape[1]:
            raise OptimaException('Characteristics are not full rank, cannot determine a unique initialization')

        # Print warning for characteristics that are not well matched by the compartment size solution
        proposed = np.matmul(A,x)
        for i in xrange(0,len(characs)):
            if abs(proposed[i]-b[i]) > project_settings.TOLERANCE:
                logger.warn('Mismatched characteristic %s %s - Requested %f, Calculated %f' % (self.label,characs[i].label,b[i],proposed[i]))
        
        # Print diagnostic output for compartments that were assigned a negative value
        def report_characteristic(charac,n_indent=0):
            if charac.label in charac_indices:
                logger.warn(n_indent * '\t' + 'Characteristic %s: Requested = %f, Calculated %f' % (charac.label,b[charac_indices[charac.label]],proposed[charac_indices[charac.label]]))
            else:
                logger.warn(n_indent * '\t' + 'Characteristic %s not in databook: Requested = N/A' % (charac.label))

            n_indent += 1
            for inc in charac.includes:
                if isinstance(inc,Characteristic):
                    report_characteristic(inc,n_indent)
                else:
                    logger.warn(n_indent * '\t' + 'Compartment %s: Computed value = %f' % (inc.label,x[comp_indices[inc.label]]))

        # Halt for any negative popsizes - print diagnostic for negative compartment
        if np.any(x < -project_settings.TOLERANCE):
            logger.warn('NEGATIVE POPSIZE DETECTED - Dumping output')
            for i in xrange(0, len(comps)):
                if x[i] < -project_settings.TOLERANCE:
                    logger.warn('Compartment %s %s - Calculated %f' % (self.label, comps[i].label, x[i]))
                    for charac in characs:
                        if comps[i] in charac.get_included_comps():
                            report_characteristic(charac)
            raise OptimaException('Negative initial popsizes')

        # Halt for an unsatisfactory overall solution (could relax this check later)
        # Print diagnostic for all characteristics
        if residual > project_settings.TOLERANCE:
            logger.warn('LARGE RESIDUAL - Dumping output')
            for charac in characs:
                report_characteristic(charac)
            raise OptimaException('Residual was %f which is unacceptably large (should be < %f) - this points to a probable inconsistency in the initial values' % (residual,project_settings.TOLERANCE))



        # Otherwise, insert the values
        for i,c in enumerate(comps):
            c.vals[0] = max(0.0,x[i])

# %% Model class
class Model(object):
    ''' A class to wrap up multiple populations within model and handle cross-population transitions. '''

    def __init__(self,settings, parset, progset=None, options=None,name=None):

        self.pops = list()              # List of population groups that this model subdivides into.
        self.pop_ids = dict()           # Maps label of a population to its position index within populations list.
        self.contacts = dict()          # Maps interactions 'from' (i.e. a->b for [a][b]) and 'into' (i.e. a<-b for [a][b]) ModelPopulations, marking them with a weight.
        self.sim_settings = odict()
        self.t_index = 0                # Keeps track of array index for current timepoint data within all compartments.
        self.programs_active = None     # True or False depending on whether Programs will be used or not
        self.pset = None                # Instance of ModelProgramSet
        self.t = None
        self.dt = None
        self.uid = uuid.uuid4()

        self.build(settings, parset, progset, options)

    def unlink(self):
        # Break cycles when deepcopying or pickling by swapping them for UIDs
        # Primary storage is in the comps, links, and outputs properties

        # If we are already unlinked, do nothing
        if self.pars_by_pop is None:
            return

        for pop in self.pops:
            pop.unlink()
        if self.pset is not None:
            self.pset.unlink()
        self.pars_by_pop = None

    def relink(self):
        # Need to enumerate objects at Model level because transitions link across pops

        # If we are already linked, do nothing
        if self.pars_by_pop is not None:
            return

        objs = {}
        for pop in self.pops:
            for obj in pop.comps + pop.characs + pop.pars + pop.links:
                objs[obj.uid] = obj

        for pop in self.pops:
            pop.relink(objs)

        self.pars_by_pop = dict()
        for pop in self.pops:
            for par in pop.pars:
                if par.label in self.pars_by_pop:
                    self.pars_by_pop[par.label].append(par)
                else:
                    self.pars_by_pop[par.label] = [par]

        if self.pset is not None:
            self.pset.relink(objs)

    def __getstate__(self):
        self.unlink()
        d = pickle.dumps(self.__dict__, protocol=-1) # Pickling to string results in a copy
        self.relink() # Relink, otherwise the original object gets unlinked
        return d

    def __setstate__(self, d):
        self.__dict__ = pickle.loads(d)
        self.relink()

    def getPop(self, pop_label):
        ''' Allow model populations to be retrieved by label rather than index. '''
        pop_index = self.pop_ids[pop_label]
        return self.pops[pop_index]

    def build(self, settings, parset, progset=None, options=None):
        ''' Build the full model. '''

        if options is None: options = dict()

        self.t = settings.tvec # NB. returning a mutable variable in a class @property method returns a new object each time
        self.dt = settings.tvec_dt

        self.sim_settings['impact_pars_not_func'] = []      # Program impact parameters that are not functions of other parameters and thus already marked for dynamic updating.
                                                            # This is only non-empty if a progset is being used in the model.

        for k, pop_label in enumerate(parset.pop_labels):
            self.pops.append(ModelPopulation(settings=settings, label=pop_label))
            self.pops[-1].preallocate(self.t,self.dt)     # Memory is allocated, speeding up model. However, values are NaN so as to enforce proper parset value saturation.
            self.pop_ids[pop_label] = k
            self.pops[-1].initialize_compartments(parset,settings,self.t[0])

        self.contacts = dcp(parset.contacts)    # Simple propagation of interaction details from parset to model.

        # Propagating cascade parameter parset values into ModelPops. Handle both 'tagged' links and 'untagged' dependencies.
        for cascade_par in parset.pars['cascade']:
            for pop_label in parset.pop_labels:
                pop = self.getPop(pop_label)
                par = pop.getPar(cascade_par.label) # Find the parameter with the requested label
                par.vals = cascade_par.interpolate(tvec=self.t, pop_label=pop_label)
                par.scale_factor = cascade_par.y_factor[pop_label]
                if par.links:
                    par.units = cascade_par.y_format[pop_label]

        # Propagating transfer parameter parset values into Model object.
        # For each population pair, instantiate a Parameter with the values from the databook
        # For each compartment, instantiate a set of Links that all derive from that Parameter
        # NB. If a Program somehow targets the transfer parameter, those values will automatically
        for trans_type in parset.transfers:
            if parset.transfers[trans_type]:
                for pop_source in parset.transfers[trans_type]:

                    transfer_parameter = parset.transfers[trans_type][pop_source] # This contains the data for all of the destrination pops

                    pop = self.getPop(pop_source)


                    for pop_target in transfer_parameter.y:

                        # Create the parameter object for this link (shared across all compartments)
                        par_label = trans_type + '_' + pop_source + '_to_' + pop_target # e.g. 'aging_0-4_to_15-64'
                        par = Parameter(label=par_label)
                        par.preallocate(self.t,self.dt)
                        val = transfer_parameter.interpolate(tvec=self.t, pop_label=pop_target)
                        par.vals = val
                        par.scale_factor = transfer_parameter.y_factor[pop_target]
                        par.units = transfer_parameter.y_format[pop_target]
                        pop.pars.append(par)
                        pop.par_lookup[par_label] = par

                        target_pop_obj = self.getPop(pop_target)

                        for source in pop.comps:
                            if not (source.tag_birth or source.tag_dead or source.is_junction):
                                # Instantiate a link between corresponding compartments
                                dest = target_pop_obj.getComp(source.label) # Get the corresponding compartment
                                link_tag = par_label + '_' + source.label # e.g. 'aging_0-4_to_15-64_sus'
                                link = Link(par, source, dest, link_tag)
                                link.preallocate(self.t,self.dt)
                                pop.links.append(link)
                                if link.label in pop.link_lookup:
                                    pop.link_lookup[link.label].append(link)
                                else:
                                    pop.link_lookup[link.label] = [link]

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
                else:
                    self.sim_settings['constraints'] = None
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
                alloc = self.pset.get_alloc(self.t,self.dt,self.sim_settings)
                self.pset.update_cache(alloc,self.t,self.dt) # Perform precomputations

            else:
                raise OptimaException('ERROR: A model run was initiated with instructions to activate programs, but no program set was passed to the model.')
        else:
            self.programs_active = False

        # set up sim_settings for later use wrt population tags
        self.sim_settings['tag_birth'] = []
        self.sim_settings['tag_death'] = []
        self.sim_settings['tag_no_plot'] = []
        for node_label in settings.node_specs:
            for tag in ['tag_birth', 'tag_death', 'tag_no_plot']:
                if tag in settings.node_specs[node_label]:
                    self.sim_settings[tag] = node_label

    def process(self, settings, progset,full_output):
        ''' 
        Run the full model.
        '''

        assert self.t_index == 0  # Only makes sense to process a simulation once, starting at ti=0

        self.update_pars(settings=settings, do_special=False)
        self.update_junctions(settings=settings)
        self.update_pars(settings=settings)
        self.update_links(settings=settings)

        for t in self.t[1:]:
            self.update_comps()
            self.t_index += 1  # Step the simulation forward
            self.update_junctions(settings=settings)
            self.update_pars(settings=settings)
            self.update_links(settings=settings)

        for pop in self.pops:

            # First, clear characteristic internal vals so non-dependent characteristics correctly appear in the parameter update
            for charac in pop.characs:
                charac.internal_vals = None # Wipe out characteristic vals to save space

            [par.update() for par in pop.pars if not par.dependency] # Update any remaining parameters

            if not full_output:
                for par in pop.pars:
                    if (not par.label in settings.linkpar_specs) or (not 'output' in settings.linkpar_specs[par.label]) or (settings.linkpar_specs[par.label] != 'y'):
                        par.vals = None
                        for link in par.links:
                            link.vals = None

        return self.pops, self.sim_settings

    def update_links(self, settings):
        '''
        Evolve model characteristics by one timestep (defaulting as 1 year).
        Each application of this method writes calculated values to the next position in popsize arrays, regardless of dt.
        Thus the corresponding time vector associated with variable dt steps must be tracked externally.
        '''

        ti = self.t_index

        for pop in self.pops:

            # First, populate all of the link values without any outflow constraints
            for par in pop.pars:
                if par.links:
                    transition = par.vals[ti]

                    if not transition:
                        for link in par.links:
                            link.vals[ti] = 0.0
                        continue
                    elif par.units == 'fraction':
                        if transition > 1.0:
                            validation_level = settings.validation['transition_fraction']
                            if validation_level == project_settings.VALIDATION_ERROR:
                                raise OptimaException(warning)
                            elif validation_level == project_settings.VALIDATION_WARN:
                                warning = "(t=%.2f) Link %s-%s (%s) has transition value = %.3f (>1)" % (self.t[ti], link.source.label, link.dest.label, pop.label, transition)
                                logger.warn(warning)
                            transition = 1.0
                        converted_frac = 1 - (1 - transition) ** self.dt  # A formula for converting from yearly fraction values to the dt equivalent.

                        for link in par.links:
                            if link.source.tag_birth:
                                n_alive = 0
                                for p in self.pops:
                                    n_alive += p.popsize(ti)
                                link.vals[ti] = n_alive * converted_frac
                            else:
                                link.vals[ti] = link.source.vals[ti] * converted_frac
                    elif par.units == 'number':
                        converted_amt = transition * self.dt
                        if len(par.links) > 1:
                            for link in par.links:
                                link.vals[ti] = converted_amt * link.source.vals[ti] / par.source_popsize(ti)
                        else:
                            par.links[0].vals[ti] = converted_amt

            # Then, adjust outflows to prevent negative popsizes
            for comp_source in pop.comps:
                if not (comp_source.is_junction or comp_source.tag_birth):
                    outflow = 0.0
                    for link in comp_source.outlinks:
                        outflow += link.vals[ti]

                    if outflow > comp_source.vals[ti]:
                        validation_level = settings.validation['negative_population']
                        warning = "(t=%.2f) Negative value encountered for: (%s - %s) at ti=%g : popsize = %g, outflow = %g" % (self.t[ti],pop.label,comp_source.label,ti,comp_source.vals[ti],outflow)
                        if validation_level == project_settings.VALIDATION_ERROR:
                            raise OptimaException(warning)
                        elif validation_level == project_settings.VALIDATION_WARN:
                            logger.warn(warning)

                        rescale = comp_source.vals[ti]/outflow
                        for link in comp_source.outlinks:
                            link.vals[ti] *= rescale


    def update_comps(self):
        """
        Set the compartment values at self.t_index+1 based on the current values at self.t_index
        and the link values at self.t_index. Values are updated by iterating over all outgoing links

        """

        ti = self.t_index

        # Pre-populate the current value - need to iterate over pops here because transfers
        # will cross population boundaries
        for pop in self.pops:
            for comp in pop.comps:
                comp.vals[ti + 1] = comp.vals[ti]

        for pop in self.pops:
            for comp in pop.comps:
                if not comp.is_junction:
                    for link in comp.outlinks:
                        if link.vals[ti]:
                            link.source.vals[ti + 1] -= link.vals[ti]
                            link.dest.vals[ti + 1] += link.vals[ti]

        # Guard against populations becoming negative due to numerical artifacts
        for pop in self.pops:
            for comp in pop.comps:
                comp.vals[ti + 1] = max(0, comp.vals[ti + 1])

    def update_junctions(self, settings):
        '''
        For every compartment considered a junction, propagate the contents onwards until all junctions are empty.
        '''
        # This is run once per timestep...all compartments are resolved in update_comps() so for the links to behave
        # consistently, we consider the junctions to be updated BEFORE stepping t_index forward
        # At the very first timestep, the junctions MIGHT already have people in them
        # In that case, those people should be flushed to

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

                    if review_count == 1:
                        for link in junc.outlinks:
                            link.vals[ti_link] = 0

                    # If the compartment is numerically empty, make it empty
                    if junc.vals[ti] <= project_settings.TOLERANCE:   # Includes negative values.
                        junc.vals[ti] = 0
                    else:
                        current_size = junc.vals[ti]
                        denom_val = sum(link.parameter.vals[ti_link] for link in junc.outlinks) # This is the total fraction of people requested to leave, so outflows are scaled to the entire compartment size
                        if denom_val == 0:
                            raise OptimaException('ERROR: Proportions for junction "%s" outflows sum to zero, resulting in a nonsensical ratio. There may even be (invalidly) no outgoing transitions for this junction.' % junc.label)
                        for link in junc.outlinks:
                            flow = current_size * link.parameter.vals[ti_link] / denom_val
                            link.source.vals[ti] -= flow
                            link.dest.vals[ti]   += flow
                            link.vals[ti_link] += flow
                            if link.dest.is_junction:
                                review_required = True # Need to review if a junction received an inflow at this step

    def update_pars(self, settings, do_special=True):
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
        do_program_overwrite = self.programs_active and self.t[ti] >= self.sim_settings['progs_start'] and self.t[ti] <= self.sim_settings['progs_end']
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

                old_vals = {par.uid: par.vals[ti] for par in self.pars_by_pop[par_label]}

                rule = settings.linkpar_specs[par_label]['rules']
                for pop in self.pops:
                    if rule == 'avg_contacts_in':

                        
                        # par_vals = [par.vals[ti] for par in self.pars_by_pop[par_label]]
                        # contacts = np.eye(6) * 9 + np.ones((6, 6))
                        # popsizes = [pop.getCharac(settings.charac_pop_count).internal_vals[ti] for pop in self.pops]
                        # par_vals = np.array(par_vals).reshape(6, 1)
                        # popsizes = np.array(popsizes).reshape(6, 1)
                        # weights = contacts * popsizes
                        # weights = weights / np.sum(weights, axis=0)
                        # par_vals = np.matmul(weights.T, par_vals)
                        # for par, val in zip(self.pars_by_pop[par_label], par_vals):
                        #     par.vals[ti] = val



                        from_list = self.contacts['into'][pop.label].keys()

                        # If interactions with a pop are initiated by the same pop, no need to proceed with special calculations. Else, carry on.
                        if not ((len(from_list) == 1 and from_list[0] == pop.label)):

                            if len(from_list) == 0:
                                new_val = 0.0
                            else:
                                val_sum = 0.0
                                weights = 0.0

                                for k,from_pop in enumerate(from_list):
                                    # All transition links with the same par_label are identically valued. For calculations, only one is needed for reference.
                                    par = self.getPop(from_pop).getPar(par_label)
                                    weight = self.contacts['into'][pop.label][from_pop]*self.getPop(from_pop).getCharac(settings.charac_pop_count).internal_vals[ti]
                                    val_sum += old_vals[par.uid]*weight
                                    weights += weight

                                if abs(val_sum) > project_settings.TOLERANCE:
                                    new_val = val_sum / weights
                                else:
                                    new_val = 0.0   # Only valid because if the weighted sum is zero, all pop_counts must be zero, meaning that the numerator is zero.

                            # Update the parameter's value in this population - will propagate to links in next stage
                            pop.getPar(par_label).vals[ti] = new_val

            # Restrict the parameter's value if a limiting range was defined
            for par in pars:
                par.constrain(ti)

def runModel(settings, parset, progset=None, options=None,full_output=True,name=None):
    '''
    Processes the TB epidemiological model.
    Parset-based overwrites are generally done externally, so the parset is only used for model-building.
    Progset-based overwrites take place internally and must be part of the processing step.
    The options dictionary is usually passed in with progset to specify when the overwrites take place.
    - If full_output = False, non-output Parameters (and corresponding links) will be set to None
    '''

    m = Model(settings, parset, progset, options)
    m.process(settings, progset,full_output)

    results = ResultSet(m, parset, settings, progset, options, name)

    return results
