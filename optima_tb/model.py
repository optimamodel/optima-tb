# %% Imports

from optima_tb.utils import flattenDict, odict, OptimaException
from optima_tb.validation import checkNegativePopulation, checkTransitionFraction
import optima_tb.settings as project_settings
from optima_tb.results import ResultSet
from optima_tb.parsing import FunctionParser


import logging
logger = logging.getLogger(__name__)

import numpy as np
from copy import deepcopy as dcp



# %% Abstract classes used in model

class Node(object):
    ''' Lightweight abstract class to represent one object within a network. '''
    def __init__(self, label='default', index=0):
        self.label = label          # Reference name for this object.
        self.index = index          # Index to denote storage position in Model object. Format is left unspecified.
                                    # For a ModelPopulation, this is intended as an integer.
                                    # For a ModelCompartment, this is intended as a tuple with first element denoting index of ModelPopulation.
        self.num_outlinks = 0       # Tracks number of nodes this one is linked to as an initial node.
        self.outlink_ids = []       # List of tupled indices corresponding to each outgoing link.
        self.inlink_ids = []        # List of tupled indices corresponding to each ingoing link.
                                    # Note that the second element of each tuple is the storage id of the link.
                                    # The first element of each tuple is the storage id of the superset container (e.g. ModelPopulation)

    def makeLinkTo(self, other_node, link_index, link_label, is_transfer=False):
        '''
        Link this node to another (i.e. create a transition link).
        Must provide an index that, by intention, uniquely represents where this link is positioned in its storage container.
        Must also provide the variable label for the link, in case Settings.linkpar_specs[link_label] need to be referred to.
        '''
        if not isinstance(other_node, Node):
            raise OptimaException('ERROR: Attempting to link compartment to something that is not a compartment.')
        self.num_outlinks += 1
        self.outlink_ids.append(link_index)
        other_node.inlink_ids.append(link_index)
        return Link(object_from=self, object_to=other_node, label=link_label, is_transfer=is_transfer)

class Variable(object):
    '''
    Lightweight abstract class to store a variable array of values (presumably corresponding to an external time vector).
    Includes an attribute to describe the format of these values.
    '''
    def __init__(self, label='default', val=0.0):
        self.label = label
        self.vals = np.array([float(val)])  # An abstract array of values.
        self.vals_old = None                # An optional array that stores old values in the case of overwriting.
        if val > 1:
            self.val_format = 'number'
        else:
            self.val_format = 'fraction'


class Link(Variable):
    '''
    A more involved version of Variable, representing unidirectional flow between two objects in a network.
    The values stored in this extended version of Variable refer to flow rates.
    If used in ModelPop, the Link references two cascade compartments within a single population.
    '''
    def __init__(self, object_from, object_to, label='default', val=0.0, scale_factor=1.0, is_transfer=False):
        Variable.__init__(self, label=label, val=val)
        self.index_from = object_from.index
        self.index_to = object_to.index
        self.label_from = object_from.label
        self.label_to = object_to.label
        self.scale_factor = scale_factor
        self.is_transfer = is_transfer

    def __repr__(self, *args, **kwargs):
        print "self.scale_factor = ", self.scale_factor, type(self.scale_factor)
        return "Link %s: from (%s, %s) to (%s, %s) with scale_factor=%g" % (self.label, self.index_from, self.label_from, self.index_to, self.label_to, self.scale_factor)

# %% Cascade compartment and population classes

class ModelCompartment(Node):
    ''' A class to wrap up data for one compartment within a cascade network. '''

    def __init__(self, label='default', index=0, popsize=0.0):
        Node.__init__(self, label=label, index=index)
        self.popsize = np.array([float(popsize)])   # Number of people in compartment.
        self.popsize_old = None                     # An array that stores old popsize values in the case of junctions (and perhaps other overwriting).
        self.tag_birth = False                      # Tag for whether this compartment contains unborn people.
        self.tag_dead = False                       # Tag for whether this compartment contains dead people.
        self.junction = False

    def __repr__(self, *args, **kwargs):
        return "%s: %g" % (self.label, self.popsize[0])

    def getValue(self, ti):
        """ Get value of population at timestep ti """
        return self.popsize[ti]

class ModelPopulation(Node):
    '''
    A class to wrap up data for one population within model.
    Each model population must contain a set of compartments with equivalent labels.
    '''

    def __init__(self, settings, label='default', index=0):
        Node.__init__(self, label=label, index=index)
        self.comps = list()         # List of cascade compartments that this model population subdivides into.
        self.links = list()         # List of intra-population cascade transitions within this model population.
        self.deps = list()          # List of value arrays for link-dependencies (both characteristic and untagged parameters) required by model.
        self.comp_ids = dict()      # Maps label of a compartment to its position index within compartments list.
        self.link_ids = dict()      # Maps cascade transition tag to indices for all relevant transitions within links list.
        self.dep_ids = dict()       # Maps label of a relevant characteristic/parameter to its position index within dependencies list.
        self.t_index = 0            # Keeps track of array index for current timepoint data within all compartments.
        self.id = index             # A numerical id of the Population. Ideally corresponds to its storage index within Model container.

        self.genCascade(settings=settings)    # Convert compartmental cascade into lists of compartment and link objects.

    def __repr__(self, *args, **kwargs):
        return "".join("%s" % self.comps)

    def getModelState(self, ti):
        states = [c.getValue(ti) for c in self.comps]
        return states

    def getComp(self, comp_label):
        ''' Allow compartments to be retrieved by label rather than index. Returns a ModelCompartment. '''
        comp_index = self.comp_ids[comp_label]
        return self.comps[comp_index]

    def getLinks(self, link_tag):
        ''' Allow links to be retrieved by tag rather than index. Returns a list of Links. '''
        link_index_list = self.link_ids[link_tag]
        link_list = []
        for link_index in link_index_list:
            link_list.append(self.links[link_index])
        return link_list

    def getDep(self, dep_label):
        ''' Allow dependencies to be retrieved by label rather than index. Returns a Variable. '''
        dep_index = self.dep_ids[dep_label]
        return self.deps[dep_index]

    def genCascade(self, settings):
        '''
        Generate a compartmental cascade as defined in a settings object.
        Fill out the compartment, transition and dependency lists within the model population object.
        Maintaining order as defined in a cascade workbook is crucial due to cross-referencing.
        '''
        for k, label in enumerate(settings.node_specs.keys()):
            self.comps.append(ModelCompartment(label=label, index=(self.index, k)))
            if 'tag_birth' in settings.node_specs[label]:
                self.comps[-1].tag_birth = True
            if 'tag_dead' in settings.node_specs[label]:
                self.comps[-1].tag_dead = True
            if 'junction' in settings.node_specs[label]:
                self.comps[-1].junction = True
            self.comp_ids[label] = k

        k = 0
        # Create actual link objects for parameters with tags, a.k.a. transitions.
        for label in settings.linkpar_specs:
            if 'tag' in settings.linkpar_specs[label]:
                tag = settings.linkpar_specs[label]['tag']
                for pair in settings.links[tag]:
                    self.links.append(self.getComp(pair[0]).makeLinkTo(self.getComp(pair[1]), link_index=(self.id, k), link_label=label))
                    if not tag in self.link_ids:
                        self.link_ids[tag] = []
                    self.link_ids[tag].append(k)
                    k += 1

        k = 0
        # Now create variable objects for dependencies, starting with characteristics.
        for label in settings.charac_specs:
            if 'par_dependency' in settings.charac_specs[label]:
                self.deps.append(Variable(label=label))
                self.dep_ids[label] = k
                k += 1

        # Finish dependencies with untagged parameters, i.e. ones that are not considered transition flow rates.
        for label in settings.par_deps:
            self.deps.append(Variable(label=label))
            self.dep_ids[label] = k
            k += 1

    def preAllocate(self, sim_settings):
        '''
        Pre-allocate variable arrays in compartments, links and dependent variables for faster processing.
        Array maintains initial value but pre-fills everything else with NaNs.
        Thus errors due to incorrect parset value saturation should be obvious from results.
        '''
        for comp in self.comps:
            init_popsize = comp.popsize[0]
            comp.popsize = np.ones(len(sim_settings['tvec'])) * np.nan
            comp.popsize[0] = init_popsize
        for link in self.links:
            init_val = link.vals[0]
            link.vals = np.ones(len(sim_settings['tvec'])) * np.nan
            link.vals[0] = init_val
        for dep in self.deps:
            init_val = dep.vals[0]
            dep.vals = np.ones(len(sim_settings['tvec'])) * np.nan
            dep.vals[0] = init_val



# %% Model class

class Model(object):
    ''' A class to wrap up multiple populations within model and handle cross-population transitions. '''

    def __init__(self):

        self.pops = list()              # List of population groups that this model subdivides into.
        self.pop_ids = dict()           # Maps label of a population to its position index within populations list.

        self.contacts = dict()          # Maps interactions 'from' (i.e. a->b for [a][b]) and 'into' (i.e. a<-b for [a][b]) ModelPopulations, marking them with a weight.

        self.sim_settings = odict()

        self.t_index = 0                # Keeps track of array index for current timepoint data within all compartments.

        self.parser = FunctionParser(debug=False)  # Decomposes and evaluates functions written as strings, in accordance with a grammar defined within the parser object.

        self.prog_vals = dict()         # Stores coverage and impact values for programs, given budget info passed into model.
                                        # Intended to avoid constant getImpact() calls during parameter value overwrites.


    def getPop(self, pop_label):
        ''' Allow model populations to be retrieved by label rather than index. '''
        pop_index = self.pop_ids[pop_label]
        return self.pops[pop_index]


    def preCalculateProgsetVals(self, settings, progset):
        ''' Work out program coverages and impacts ahead of the model run. '''

        try: start_year = self.sim_settings['progs_start']
        except: raise OptimaException('ERROR: Pre-calculation of program set values has been initiated without specifying a start year.')

        try: init_alloc = self.sim_settings['init_alloc']
        except: raise OptimaException('ERROR: Pre-calculation of program set values has been initiated without specifying a starting allocation, empty or otherwise.')

        try: alloc_is_coverage = self.sim_settings['alloc_is_coverage']
        except: raise OptimaException('ERROR: Pre-calculation of program set values has been initiated without specifying whether starting allocation is in money or coverage.')

        for prog in progset.progs:

            # Store budgets/coverages for programs that are initially allocated.
            if prog.label in init_alloc:
                alloc = init_alloc[prog.label]

                # If ramp constraints are active, stored cost and coverage needs to be a fully time-dependent array corresponding to timevec.
                if 'constraints' in self.sim_settings and 'max_yearly_change' in self.sim_settings['constraints'] and prog.label in self.sim_settings['constraints']['max_yearly_change']:
                    if alloc_is_coverage:
                        default = prog.getCoverage(budget=prog.getDefaultBudget(year=start_year))
                    else:
                        default = prog.getDefaultBudget(year=start_year)
                    default = prog.getDefaultBudget(year=start_year)
                    if np.abs(alloc - default) > project_settings.TOLERANCE:
#                        print 'Start it up...'
                        alloc_def = self.sim_settings['tvec'] * 0.0 + default
                        alloc_new = self.sim_settings['tvec'] * 0.0 + alloc
                        try: eps = self.sim_settings['constraints']['max_yearly_change'][prog.label]['val']
                        except: raise OptimaException('ERROR: A maximum yearly change constraint was passed to the model for "%s" but had no value associated with it.' % prog.label)
                        if 'rel' in self.sim_settings['constraints']['max_yearly_change'][prog.label] and self.sim_settings['constraints']['max_yearly_change'][prog.label]['rel'] is True:
                            eps *= default
                        if np.isnan(eps): eps = np.inf  # Eps likely becomes a nan if it was infinity multiplied by zero.
                        if np.abs(eps * settings.tvec_dt) < np.abs(alloc - default):
                            if np.abs(eps) < project_settings.TOLERANCE:
                                raise OptimaException('ERROR: The change in budget for ramp-constrained "%s" is effectively zero. Model will not continue running; change in program funding would be negligible.' % prog.label)
                            alloc_ramp = default + (self.sim_settings['tvec'] - start_year) * eps * np.sign(alloc - default)
                            if alloc >= default: alloc_ramp = np.minimum(alloc_ramp, alloc_new)
                            else: alloc_ramp = np.maximum(alloc_ramp, alloc_new)
                            alloc = alloc_def * (self.sim_settings['tvec'] < start_year) + alloc_ramp * (self.sim_settings['tvec'] >= start_year)

                if alloc_is_coverage:
                    self.prog_vals[prog.label] = {'cost':prog.getBudget(coverage=alloc), 'cov':alloc, 'impact':{}}
                else:
                    self.prog_vals[prog.label] = {'cost':alloc, 'cov':prog.getCoverage(budget=alloc), 'impact':{}}

            # Store default budgets/coverages for all other programs if saturation is selected.
            elif 'saturate_with_default_budgets' in self.sim_settings and self.sim_settings['saturate_with_default_budgets'] is True:
                self.prog_vals[prog.label] = {'cost':prog.getDefaultBudget(year=start_year), 'cov':prog.getCoverage(budget=prog.getDefaultBudget(year=start_year)), 'impact':{}}

            else:
                logger.warn("Program '%s' was not contained in init_alloc and not saturated, therefore was not created." % prog.label)

            # Convert coverage into impact for programs.
            if prog.label in self.prog_vals:
                if 'cov' in self.prog_vals[prog.label]:
                    cov = self.prog_vals[prog.label]['cov']

                    # Check if program attributes have multiple distinct values across time.
                    do_full_tvec_check = {}
                    for att_label in prog.attributes:
                        att_vals = prog.attributes[att_label]
                        if len(set(att_vals[~np.isnan(att_vals)])) <= 1:
                            do_full_tvec_check[att_label] = False
                        else:
                            do_full_tvec_check[att_label] = True

                    # If attributes change over time and coverage values are greater than zero, impact functions based on them need to be interpolated across all time.
                    # Otherwise interpolating for the very last timestep alone should be sufficient.
                    # This means impact values can be stored as full arrays, single-element arrays or scalars in the case that an impact function has no attributes.
                    for par_label in prog.target_pars:
                        do_full_tvec = False
                        if 'attribs' in prog.target_pars[par_label] and np.sum(cov) > project_settings.TOLERANCE:
                            for att_label in prog.target_pars[par_label]['attribs']:
                                if att_label in do_full_tvec_check and do_full_tvec_check[att_label] is True:
                                    do_full_tvec = True
                        if do_full_tvec is True:
                            years = self.sim_settings['tvec']
                        else:
                            years = [self.sim_settings['tvec'][-1]]
                        self.prog_vals[prog.label]['impact'][par_label] = prog.getImpact(cov, impact_label=par_label, parser=self.parser, years=years, budget_is_coverage=True)


    def build(self, settings, parset, progset=None, options=None):
        ''' Build the full model. '''

        if options is None: options = dict()

        self.sim_settings['tvec'] = np.arange(settings.tvec_start, settings.tvec_end + settings.tvec_dt / 2, settings.tvec_dt)
        self.sim_settings['impact_pars_not_func'] = []      # Program impact parameters that are not functions of other parameters and thus already marked for dynamic updating.
                                                            # This is only non-empty if a progset is being used in the model.
        if 'progs_start' in options:
            if progset is not None:
                self.sim_settings['progs_start'] = options['progs_start']

                if 'progs_end' in options:
                    self.sim_settings['progs_end'] = options['progs_end']
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

                self.preCalculateProgsetVals(settings=settings, progset=progset)   # For performance.
            else:
                raise OptimaException('ERROR: A model run was initiated with instructions to activate programs, but no program set was passed to the model.')

        self.parser.debug = settings.parser_debug

        for k, pop_label in enumerate(parset.pop_labels):
            self.pops.append(ModelPopulation(settings=settings, label=pop_label, index=k))
            self.pops[-1].preAllocate(self.sim_settings)     # Memory is allocated, speeding up model. However, values are NaN so as to enforce proper parset value saturation.
            self.pop_ids[pop_label] = k

        self.contacts = dcp(parset.contacts)    # Simple propagation of interaction details from parset to model.

        # Propagating initial characteristic parset values into ModelPops.
        # NOTE: Extremely involved process, so might be worth extracting the next few paragraphs as a separate method.

        t_init = np.array([self.sim_settings['tvec'][0]])
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

        # Now initialise all model compartments with these calculated values.
        for seed_label in seed_dict:
            for pop_label in parset.pop_labels:
                val = seed_dict[seed_label][pop_label]
                if abs(val) < project_settings.TOLERANCE:
                    val = 0
                elif val < 0.:
                    raise OptimaException('ERROR: Initial value calculated for compartment "%s" in population "%s" is %f. Review and make sure each characteristic has at least as many people as the sum of all included compartments.' % (seed_label, pop_label, val))
                self.getPop(pop_label).getComp(seed_label).popsize[0] = val



        # Propagating cascade parameter parset values into ModelPops. Handle both 'tagged' links and 'untagged' dependencies.
        for par in parset.pars['cascade']:

            if 'tag' in settings.linkpar_specs[par.label]:
                tag = settings.linkpar_specs[par.label]['tag']                  # Map parameter label to link tag.
                for pop_label in parset.pop_labels:
                    for link_id in self.getPop(pop_label).link_ids[tag]:        # Map link tag to link id in ModelPop.
                        self.getPop(pop_label).links[link_id].vals = par.interpolate(tvec=self.sim_settings['tvec'], pop_label=pop_label)
                        self.getPop(pop_label).links[link_id].val_format = par.y_format[pop_label]
                        self.getPop(pop_label).links[link_id].scale_factor = par.y_factor[pop_label]
                        
                # Apply min/max restrictions on all parameters that are not functions.
                # Functional parameters will be calculated and constrained during a run, hence they can be np.nan at this stage.
                if not par.label in settings.par_funcs.keys():
                    if 'min' in settings.linkpar_specs[par.label]:
                        for pop_label in parset.pop_labels:
                            for link_id in self.getPop(pop_label).link_ids[tag]:
                                vals = self.getPop(pop_label).links[link_id].vals
                                self.getPop(pop_label).links[link_id].vals[vals < settings.linkpar_specs[par.label]['min']] = settings.linkpar_specs[par.label]['min']
                    if 'max' in settings.linkpar_specs[par.label]:
                        for pop_label in parset.pop_labels:
                            for link_id in self.getPop(pop_label).link_ids[tag]:
                                vals = self.getPop(pop_label).links[link_id].vals
                                self.getPop(pop_label).links[link_id].vals[vals > settings.linkpar_specs[par.label]['max']] = settings.linkpar_specs[par.label]['max']

            else:
                for pop_label in parset.pop_labels:
                    dep_id = self.getPop(pop_label).dep_ids[par.label]          # Map dependency label to dependency id in ModelPop.
                    self.getPop(pop_label).deps[dep_id].vals = par.interpolate(tvec=self.sim_settings['tvec'], pop_label=pop_label)
                    self.getPop(pop_label).deps[dep_id].val_format = par.y_format[pop_label]
                    self.getPop(pop_label).deps[dep_id].scale_factor = par.y_factor[pop_label]
                
                # Apply min/max restrictions on all parameters that are not functions.
                # Functional parameters will be calculated and constrained during a run, hence they can be np.nan at this stage.
                if not par.label in settings.par_funcs.keys():
                    if 'min' in settings.linkpar_specs[par.label]:
                        for pop_label in parset.pop_labels:
                            dep_id = self.getPop(pop_label).dep_ids[par.label]
                            vals = self.getPop(pop_label).deps[dep_id].vals
                            self.getPop(pop_label).deps[dep_id].vals[vals < settings.linkpar_specs[par.label]['min']] = settings.linkpar_specs[par.label]['min']
                    if 'max' in settings.linkpar_specs[par.label]:
                        for pop_label in parset.pop_labels:
                            dep_id = self.getPop(pop_label).dep_ids[par.label]
                            vals = self.getPop(pop_label).deps[dep_id].vals
                            self.getPop(pop_label).deps[dep_id].vals[vals > settings.linkpar_specs[par.label]['max']] = settings.linkpar_specs[par.label]['max']


        # Propagating transfer parameter parset values into Model object.
        for trans_type in parset.transfers:
            if parset.transfers[trans_type]:
                for pop_source in parset.transfers[trans_type]:
                    par = parset.transfers[trans_type][pop_source]

                    for pop_target in par.y:
                        for comp in self.getPop(pop_source).comps:
                            trans_tag = comp.label + '_' + trans_type + '_to_' + pop_target       # NOTE: Perhaps there is a nicer way to set up transfer tagging.
                            if not comp.tag_birth and not comp.tag_dead and not comp.junction:

                                num_links = len(self.getPop(pop_source).links)
                                link = comp.makeLinkTo(self.getPop(pop_target).getComp(comp.label), link_index=(self.getPop(pop_source).id, num_links), link_label=trans_tag, is_transfer=True)
                                link.vals = par.interpolate(tvec=self.sim_settings['tvec'], pop_label=pop_target)
                                link.val_format = par.y_format[pop_target]
                                link.scale_factor = par.y_factor[pop_target]
#                                if link.val_format == 'number': link.vals /= settings.num_transfer_nodes

                                self.getPop(pop_source).links.append(link)
                                self.getPop(pop_source).link_ids[trans_tag] = [num_links]

        # Make sure initially-filled junctions are processed and initial dependencies are calculated.
        self.updateValues(settings=settings, progset=progset, do_special=False)     # Done first just in case junctions are dependent on characteristics.
                                                                                                # No special rules are applied at this stage, otherwise calculations would be iterated twice before the first step forward.
                                                                                                # NOTE: If junction outflows were to be tagged by special rules, initial calculations may be off. Return to this later and consider logic rigorously.
        self.processJunctions(settings=settings)
        self.updateValues(settings=settings, progset=progset)


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
#            self.printModelState(self.t_index)
            self.stepForward(settings=settings, dt=settings.tvec_dt)
            self.processJunctions(settings=settings)
            self.updateValues(settings=settings, progset=progset)

        return self.pops, self.sim_settings

    def _calculateDPops(self, settings, ti, dt, modified_compartment={}, empty_compartment=[], reset_compartment=[]):
        """
        
        This function gets called from model.stepForward() but also from within a validation check. 
        
        Params:
            settings
            ti
            dt
            modified_compartment optional: Dictionary of (did,reduction_ratio) pairs that indicate by what percentage an dpop amount should be reduced.  
            empty_compartment    optional: List of compartment did that are empty
            reset_compartment    optional: List of compartment did that already have a negative value
        
        """

        # Preallocate a change-in-popsize array. Requires each population to have the same cascade.
        num_pops = len(self.pops)
        num_comps = len(self.pops[0].comps)         # NOTE: Hard-coded reference to 'zeroth' population. Improve later.
        dpopsize = np.zeros(num_pops * num_comps)
        # Variables for tracking in and out amounts to a compartment, separately
        # dpop_in = np.zeros(num_pops*num_comps)
        dpop_out = np.zeros(num_pops * num_comps)

        for pop in self.pops:
            for comp in pop.comps:

                # If not pre-allocated, extend compartment variable arrays. Either way copy current value to the next position in the array.
                if not len(comp.popsize) > ti + 1:
                    comp.popsize = np.append(comp.popsize, 0.0)
                if not len(comp.popsize) > ti + 1:      # If one extension did not create an index of ti+1, something is seriously wrong...
                    raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment "%s".' % (comp.label))
                comp.popsize[ti + 1] = comp.popsize[ti]

                if comp.label not in settings.junction_labels:      # Junctions collect inflows during this step. They do not process outflows here.
                    for lid_tuple in comp.outlink_ids:
                        lid = lid_tuple[-1]
                        link = pop.links[lid]

                        # If link values are not pre-allocated to the required time-vector length, extend with the same value as the last index.
                        if not len(link.vals) > ti + 1:
                            link.vals = np.append(link.vals, link.vals[-1])
                        if not len(link.vals) > ti + 1:         # If one extension did not create an index of ti+1, something is seriously wrong...
                            raise OptimaException('ERROR: Current timepoint in simulation does not mesh with array length in compartment "%s".' % (link.label))

                        did_from = link.index_from[0] * num_comps + link.index_from[1]
                        did_to = link.index_to[0] * num_comps + link.index_to[1]
                        comp_source = self.pops[link.index_from[0]].getComp(link.label_from)

                        converted_amt = 0
                        
#                        # TODO: Make this range check more efficient by pre-validating. Could use extra validation for user error.
#                        if link.label in settings.linkpar_specs:
#                            if 'min' in settings.linkpar_specs[link.label] and settings.linkpar_specs[link.label]['min'] > link.vals[ti]:
#                                link.vals[ti] = settings.linkpar_specs[link.label]['min']
#                            if 'max' in settings.linkpar_specs[link.label] and settings.linkpar_specs[link.label]['max'] < link.vals[ti]:
#                                link.vals[ti] = settings.linkpar_specs[link.label]['max']

                        transition = link.vals[ti]

                        if link.scale_factor is not None and link.scale_factor != project_settings.DO_NOT_SCALE : # scale factor should be available to be used
                            transition *= link.scale_factor


                        if link.val_format == 'fraction':
                            # check if there are any violations, and if so, deal with them
                            if transition > 1.:
                                transition = checkTransitionFraction(transition, settings.validation)

                            converted_frac = 1 - (1 - transition) ** dt      # A formula for converting from yearly fraction values to the dt equivalent.
                            converted_amt = comp_source.popsize[ti] * converted_frac
                        elif link.val_format == 'number':
                            converted_amt = transition * dt
                            if link.is_transfer:
                                transfer_rescale = comp_source.popsize[ti] / pop.getDep(settings.charac_pop_count).vals[ti]
                                converted_amt *= transfer_rescale

                        else:
                            raise OptimaException("Unknown link type: %s in model\nObserved for population %s, compartment %s" % (link.val_format, pop.label, comp.label))


                        # If we've observed that this link contributes to the creation of a negative population size, we modify the amounts:
                        if did_from in empty_compartment:
                            # If the compartment is an empty one, we shouldn't update
                            logging.debug("Empty compartment: no change made for link from=%g,to=%g" % (did_from, did_to))
                            converted_amt = 0

                        elif did_from in modified_compartment:
                            # during validation, we encountered an outflow from a compartment that was more than the population size
                            logging.debug("Modifying flow amount for link (from=%g,to=%g) by %g: prev amount = %g, new amount = %g" % (did_from, did_to, modified_compartment[did_from], converted_amt, modified_compartment[did_from] * converted_amt))
                            converted_amt *= modified_compartment[did_from]

                        elif did_from in reset_compartment:
                            # set converted amount in the other direction. Hmm still not perfect
                            converted_amt = 0


                        dpopsize[did_from] -= converted_amt
                        dpopsize[did_to] += converted_amt

                        # Tracking in and out for a compartment
                        # dpop_in[did_to] += converted_amt
                        dpop_out[did_from] += converted_amt

        return dpopsize, dpop_out



    def stepForward(self, settings, dt=1.0):
        '''
        Evolve model characteristics by one timestep (defaulting as 1 year).
        Each application of this method writes calculated values to the next position in popsize arrays, regardless of dt.
        Thus the corresponding time vector associated with variable dt steps must be tracked externally.
        '''

        ti = self.t_index

        # Requires each population to have the same cascade.
        num_pops = len(self.pops)
        num_comps = len(self.pops[0].comps)         # NOTE: Hard-coded reference to 'zeroth' population. Improve later.


        # First loop through all pops and compartments to calculate value changes.
        dpopsize, dpop_out = self._calculateDPops(settings, ti, dt)

        # Check calculated rates for proposed dpop value. Note that this should be made an iterative process for
        # when we avert an incorrect dpop change, as the avert may result in new errors.
        # TODO: create validation as an iterative process (currently found that it was possible that there was no good solution found
        """
        # SUGGESTION: 
        count = 0 
        while count < some_limit_count: 
            if isDPopValid(self,settings,dpopsize,ti,validation):
                break
            dpopsize,dpop_out,reset_ids = checkNegativePopulation(self,settings,dpopsize,dpop_out,ti,validation)
        """
        dpopsize, dpop_out, reset_ids = checkNegativePopulation(self, settings, dpopsize, dpop_out, ti, dt, settings.validation)

        # Second loop through all pops and comps to apply value changes at next timestep index.
        for pid in xrange(num_pops):
            for cid in xrange(num_comps):
                did = pid * num_comps + cid
                self.pops[pid].comps[cid].popsize[ti + 1] += dpopsize[did]

        # Update timestep index.
        self.t_index += 1
        for pop in self.pops:
            pop.t_index = self.t_index       # Keeps ModelPopulations synchronised with Model object.


    def processJunctions(self, settings):
        '''
        For every compartment considered a junction, propagate the contents onwards until all junctions are empty.
        '''

        ti = self.t_index
        ti_link = ti - 1
        if ti_link < 0: ti_link = ti    # For the case where junctions are processed immediately after model initialisation.
        final_review = False

        review_count = 0
        while not final_review:
            if review_count > settings.recursion_limit: raise OptimaException('ERROR: Processing junctions (i.e. propagating contents onwards) for timestep %i is taking far too long. Infinite loop suspected.' % ti_link)
            final_review = True     # Assume that this is the final sweep through junctions to verify their emptiness.
            for pop in self.pops:
                for junction_label in settings.junction_labels:
                    comp = pop.getComp(junction_label)

                    # Stores junction popsize values before emptying.
                    if review_count == 0:
                        if comp.popsize_old is None:
                            comp.popsize_old = dcp(comp.popsize)
                        else:
                            comp.popsize_old[ti] = comp.popsize[ti]
                    # If a junction is being reviewed again, it means that it received inflow before emptying.
                    # Add this inflow to the stored popsize.
                    else:
                        comp.popsize_old[ti] += comp.popsize[ti]
#                        print 'huzzah'
#                        print comp.popsize
#                        print comp.popsize_old

                    popsize = comp.popsize[ti]

                    if popsize <= project_settings.TOLERANCE:   # Includes negative values.
                        comp.popsize[ti] = 0
                        popsize = 0

                    elif popsize > project_settings.TOLERANCE:
                        final_review = False    # Outflows could propagate into other junctions requiring another review.
                        denom_val = sum(pop.links[lid_tuple[-1]].vals[ti_link] for lid_tuple in comp.outlink_ids)
                        if denom_val == 0: raise OptimaException('ERROR: Proportions for junction "%s" outflows sum to zero, resulting in a nonsensical ratio. There may even be (invalidly) no outgoing transitions for this junction.' % junction_label)
                        for lid_tuple in comp.outlink_ids:
                            lid = lid_tuple[-1]
                            link = pop.links[lid]

                            comp.popsize[ti] -= popsize * link.vals[ti_link] / denom_val
                            pop.getComp(link.label_to).popsize[ti] += popsize * link.vals[ti_link] / denom_val



            review_count += 1



    def updateValues(self, settings, progset=None, do_special=True):
        '''
        Run through all parameters and characteristics flagged as dependencies for custom-function parameters and evaluate them for the current timestep.
        These dependencies must be calculated in the same order as defined in settings, characteristics before parameters, otherwise references may break.
        Also, parameters in the dependency list do not need to be calculated unless explicitly depending on another parameter.
        Parameters that have special rules are usually dependent on other population values, so are included here.
        '''

        ti = self.t_index

        for pop in self.pops:

            # Characteristics that are dependencies first...
            for dep in pop.deps:
                if dep.label in settings.charac_deps:
                    dep.vals[ti] = 0

                    # Sum up all relevant compartment popsizes (or previously calculated characteristics).
                    for inc_label in settings.charac_specs[dep.label]['includes']:
                        if inc_label in pop.comp_ids:
                            val = pop.getComp(inc_label).popsize[ti]
                        elif inc_label in pop.dep_ids:    # NOTE: This should not select a parameter-type dependency due to settings validation, but can validate here if desired.
                            val = pop.getDep(inc_label).vals[ti]
                        else:
#                            print inc_label
#                            print settings.charac_specs[dep.label]['includes']
#                            print pop.comp_ids
#                            print pop.dep_ids
                            raise OptimaException('ERROR: Compartment or characteristic "%s" has not been pre-calculated for use in calculating "%s".' % (inc_label, dep.label))

                        dep.vals[ti] += val

                    # Divide by relevant compartment popsize (or previously calculated characteristic).
                    if 'denom' in settings.charac_specs[dep.label]:
                        den_label = settings.charac_specs[dep.label]['denom']
                        if den_label in pop.dep_ids:  # NOTE: See above note for avoiding parameter-type dependencies.
                            val = pop.getDep(den_label).vals[ti]
                        elif den_label in pop.comp_ids:
                            val = pop.getComp(den_label).popsize[ti]
                        else:
                            raise OptimaException('ERROR: Compartment or characteristic "%s" has not been pre-calculated for use in calculating "%s".' % (inc_label, dep.label))

                        if val == 0:
                            if abs(dep.vals[ti]) < project_settings.TOLERANCE:
                                dep.vals[ti] = 0        # Given a zero/zero case, make the answer zero.
                            else:
                                dep.vals[ti] = np.inf   # Given a non-zero/zero case, keep the answer infinite.
                        else:
                            dep.vals[ti] /= val

        # Handle parameters now that require calculations between timestep-based updates.
        # 1st:  Any parameter that is a function of others (i.e. of previously calculated dependencies).
        # 2nd:  Any parameter that is overwritten by a program-based cost-coverage-impact transformation.
        # 3rd:  Any parameter that is overwritten by a special rule, e.g. averaging across populations.
        # 4th:  Any parameter that is restricted within a range of values, i.e. by min/max values.
        # Looping through populations must be internal so that all values are calculated before special inter-population rules are applied.
        for par_label in (settings.par_funcs.keys() + self.sim_settings['impact_pars_not_func']):
            for pop in self.pops:
                pars = []
                if par_label in settings.par_deps:
                    pars.append(pop.getDep(par_label))
                else:
                    pars = pop.getLinks(settings.linkpar_specs[par_label]['tag'])
                specs = settings.linkpar_specs[par_label]

                # Calculate the value of a parameter from its function if it exists, otherwise maintain the current value.
                if 'f_stack' in specs:
                    f_stack = dcp(specs['f_stack'])
                    deps = dcp(specs['deps'])
                    for dep_label in deps:
                        if dep_label in settings.par_deps or dep_label in settings.charac_deps:
                            val = pop.getDep(dep_label).vals[ti]
                        else:
                            val = pop.getLinks(settings.linkpar_specs[dep_label]['tag'])[0].vals[ti]    # As links are duplicated for the same tag, can pull values from the zeroth one.
                        deps[dep_label] = val
                    new_val = self.parser.evaluateStack(stack=f_stack, deps=deps)
                else:
                    new_val = pars[0].vals[ti]      # As links are duplicated for the same tag, can pull values from the zeroth one.

                # WARNING: CURRENTLY IS NOT RELIABLE FOR IMPACT PARAMETERS THAT ARE DUPLICATE LINKS.
                if 'progs_start' in self.sim_settings and self.sim_settings['tvec'][ti] >= self.sim_settings['progs_start']:
                    if not ('progs_end' in self.sim_settings and self.sim_settings['tvec'][ti] >= self.sim_settings['progs_end']):
                        if par_label in progset.impacts:
                            first_prog = True   # True if program in prog_label loop is the first one in the impact dict list.
                            impact_list = []    # Notes for each program what its impact would be before coverage limitations.
                            overflow_list = []  # Notes for each program how much greater its funded coverage is than people available to be covered.
                            for prog_label in progset.impacts[par_label]:
                                prog = progset.getProg(prog_label)
                                prog_type = prog.prog_type

                                # Make sure the population in the loop is a target of this program.
                                if pop.label not in prog.target_pops:
                                    continue
                                if not ('init_alloc' in self.sim_settings and prog_label in self.sim_settings['init_alloc']):
                                    continue

                                # If a program target is a transition parameter, its coverage must be distributed and format-converted.
                                if 'tag' in settings.linkpar_specs[par_label]:
                                    # Coverage is assumed to be across a compartment over a set of populations, not a single element, so scaling is required.
                                    source_element_size = self.pops[pars[0].index_from[0]].comps[pars[0].index_from[1]].popsize[ti]
                                    source_set_size = 0
                                    for from_pop in prog.target_pops:
                                        source_set_size += self.getPop(from_pop).comps[pars[0].index_from[1]].popsize[ti]
    
                                    # Coverage is also split across the source compartments of grouped impact parameters, as specified in the cascade sheet.
                                    # NOTE: This might be a place to improve performance.
                                    if 'group' in settings.progtype_specs[prog_type]['impact_pars'][par_label]:
                                        group_label = settings.progtype_specs[prog_type]['impact_pars'][par_label]['group']
                                        for alt_par_label in settings.progtype_specs[prog_type]['impact_par_groups'][group_label]:
                                            if not alt_par_label == par_label:
                                                alt_pars = pop.getLinks(settings.linkpar_specs[alt_par_label]['tag'])
                                                for from_pop in prog.target_pops:
                                                    source_set_size += self.getPop(from_pop).comps[alt_pars[0].index_from[1]].popsize[ti]
    
                                    # Impact functions can be parsed/calculated as time-dependent arrays or time-independent scalars.
                                    # Makes sure that the right value is selected.
                                    try: net_impact = self.prog_vals[prog_label]['impact'][par_label][ti]
                                    except: net_impact = self.prog_vals[prog_label]['impact'][par_label]
    
                                    # Make sure each program impact is in the format of the parameter it affects.
                                    # NOTE: Can be compressed. Currently left explicit.
                                    if pars[0].val_format == 'fraction':
                                        if prog.cov_format == 'fraction':
                                            impact = net_impact
                                        elif prog.cov_format == 'number':
                                            if source_element_size <= project_settings.TOLERANCE:
                                                impact = 0.0
                                            else:
                                                impact = net_impact / source_set_size
        #                                    if impact > 1.0: impact = 1.0   # Maximum fraction allowable due to timestep conversion.
                                    elif pars[0].val_format == 'number':
                                        if prog.cov_format == 'fraction':
                                            impact = net_impact * source_element_size
                                        elif prog.cov_format == 'number':
                                            if source_element_size <= project_settings.TOLERANCE:
                                                impact = 0.0
                                            else:
                                                impact = net_impact * source_element_size / source_set_size
    
                                    # Calculate how excessive program coverage is for the provided budgets.
                                    # If coverage is a fraction, excess is compared to unity.
                                    # If coverage is a number, excess is compared to the total number of people available for coverage.
                                    overflow_factor = 0
                                    try: net_cov = self.prog_vals[prog_label]['cov'][ti]
                                    except: net_cov = self.prog_vals[prog_label]['cov']
                                    if prog.cov_format == 'fraction':
                                        overflow_factor = net_cov
                                    elif prog.cov_format == 'number':
                                        if float(source_set_size) <= project_settings.TOLERANCE:
                                            overflow_factor = np.inf
                                        else:
                                            overflow_factor = net_cov / float(source_set_size)
                                    overflow_list.append(overflow_factor)
                                    
                                # If a program target is any other parameter, the parameter value is directly overwritten by coverage.
                                # TODO: Decide how to handle coverage distribution.
                                else:
                                    try: impact = self.prog_vals[prog_label]['impact'][par_label][ti]
                                    except: impact = self.prog_vals[prog_label]['impact'][par_label]
                                    

#                                year_check = 2015   # Hard-coded check.
#                                par_check = ['spdno_rate','sndno_rate']#['spdsuc_rate','spdno_rate']
#                                if par_label in par_check:
#                                    if self.sim_settings['tvec'][ti] >= year_check and self.sim_settings['tvec'][ti] < year_check + 0.5*settings.tvec_dt:
#                                        print('Year: %s' % self.sim_settings['tvec'][ti])
#                                        print('Program Name: %s' % prog.name)
#                                        print('Program %s: %f' % ('Coverage' if self.sim_settings['alloc_is_coverage'] else 'Budget', self.prog_vals[prog_label]['cost']))
#                                        print('Target Population: %s' % pop.label)
#                                        print('Target Parameter: %s' % par_label)
#                                        print('Unit Cost: %f' % prog.func_specs['pars']['unit_cost'])
#                                        print('Program Coverage: %f' % self.prog_vals[prog_label]['cov'])
#                                        print('Program Impact: %f' % net_impact)
#                                        print('Program Impact Format: %s' % prog.cov_format)
#                                        print('Source Compartment Size (Target Pop): %f' % source_element_size)
#                                        print('Source Compartment Size (Aggregated Over Target Pops): %f' % source_set_size)
#                                        print('Converted Impact: %f' % impact)
#                                        print('Converted Impact Format: %s' % pars[0].val_format)
#                                        print

                                if first_prog:
                                    new_val = 0  # Zero out the new impact parameter for the first program that targets it within an update, just to make sure the overwrite works.
                                    first_prog = False
#                                new_val += impact
                                impact_list.append(impact)

#                            print overflow_list

                            # Checks to make sure that the net coverage of all programs targeting a parameters is capped by those that are available to be covered.
                            # Otherwise renormalises impacts.
                            if len(overflow_list) > 0 and sum(overflow_list) > 1:
                                impact_list = np.multiply(impact_list, 1 / sum(overflow_list))


                            new_val += np.sum(impact_list)

                            # Handle impact constraints.
                            # Note: This applies to any parameter that is impacted by the progset, not just for programs that are in the allocation.
                            if 'constraints' in self.sim_settings and 'impacts' in self.sim_settings['constraints'] and par_label in self.sim_settings['constraints']['impacts']:
                                try: vals = self.sim_settings['constraints']['impacts'][par_label]['vals']
                                except: raise OptimaException('ERROR: An impact constraint was passed to the model for "%s" but had no values associated with it.' % par_label)
                                if not len(vals) == 2: raise OptimaException('ERROR: Constraints for impact "%s" must be provided as a list or tuple of two values, i.e. a lower and an upper constraint.' % par_label)
                                if new_val < vals[0]: new_val = vals[0]
                                if new_val > vals[1]: new_val = vals[1]

                for par in pars:

#                    year_check = 2015   # Hard-coded check.
#                    par_check = ['spdno_rate','sndno_rate']#['spdsuc_rate','spdno_rate']
#                    if par_label in par_check:
#                        if self.sim_settings['tvec'][ti] >= year_check and self.sim_settings['tvec'][ti] < year_check + 0.5*settings.tvec_dt:
#                            print('Year: %s' % self.sim_settings['tvec'][ti])
#                            print('Target Population: %s' % pop.label)
#                            print('Target Parameter: %s' % par_label)
#                            print('Final Impact: %f' % new_val)
#                            print

                    par.vals[ti] = new_val

                    # Backup the values of parameters that are tagged with special rules.
                    if 'rules' in settings.linkpar_specs[par_label]:
                        if par.vals_old is None: par.vals_old = dcp(par.vals)
                        par.vals_old[ti] = new_val

            # Handle parameters tagged with special rules. Overwrite vals if necessary.
            if do_special and 'rules' in settings.linkpar_specs[par_label]:
                rule = settings.linkpar_specs[par_label]['rules']
                for pop in self.pops:
                    if rule == 'avg_contacts_in':
                        from_list = self.contacts['into'][pop.label]

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
                                    if par_label in settings.par_deps:
                                        par = self.getPop(from_pop).getDep(par_label)
                                    else:
                                        par = self.getPop(from_pop).getLinks(settings.linkpar_specs[par_label]['tag'])[0]
                                    old_vals[k] = par.vals_old[ti]
                                    weights[k] = self.contacts['into'][pop.label][from_pop]
                                    pop_counts[k] = self.getPop(from_pop).getDep(settings.charac_pop_count).vals[ti]
                                    k += 1
                                wpc = np.multiply(weights, pop_counts)          # Population counts weighted by contact rate.
                                wpc_sum = sum(wpc)                              # Normalisation factor for weighted population counts.
                                if abs(wpc_sum) > project_settings.TOLERANCE:
                                    new_val = np.dot(old_vals, wpc / wpc_sum)         # Do a weighted average of the parameter values pertaining to contact-initiating pop groups.
                                else:
                                    new_val = 0.0   # Only valid because if the weighted sum is zero, all pop_counts must be zero, meaning that the numerator is zero.

#                            if ti == 0:
#                                print
#                                print('Timestep: %s' % ti)
#                                print('Weighted population-contact averaging will affect "%s" for "%s".' % (par_label, pop.label))
#                                print('Populations initiating contact with "%s"...' % pop.label)
#                                print from_list
#                                print('These populations have "%s" values of...' % par_label)
#                                print old_vals
#                                print('Their pop counts are...')
#                                print pop_counts
#                                print('These are further weighted by...')
#                                print weights
#                                print('The new weighted population-contact average value of "%s" for "%s" is: %f' % (par_label, pop.label, new_val))

                            # Need to update all untagged/tagged links with the new value, hence the list of links.
                            pars = []
                            if par_label in settings.par_deps:
                                pars.append(pop.getDep(par_label))
                            else:
                                pars = pop.getLinks(settings.linkpar_specs[par_label]['tag'])
                            for par in pars:
                                par.vals[ti] = new_val
            
            # Do a final range check on a parameter, if applicable, and restrict its value.
            # TODO: Explore efficiency gains.                
            if 'min' in settings.linkpar_specs[par_label]:
                for pop in self.pops:
                    pars = []
                    if par_label in settings.par_deps:
                        pars.append(pop.getDep(par_label))
                    else:
                        pars = pop.getLinks(settings.linkpar_specs[par_label]['tag'])
                    for par in pars:
                        if settings.linkpar_specs[par_label]['min'] > par.vals[ti]:
                            par.vals[ti] = settings.linkpar_specs[par_label]['min']
            if 'max' in settings.linkpar_specs[par_label]:
                for pop in self.pops:
                    pars = []
                    if par_label in settings.par_deps:
                        pars.append(pop.getDep(par_label))
                    else:
                        pars = pop.getLinks(settings.linkpar_specs[par_label]['tag'])
                    for par in pars:
                        if settings.linkpar_specs[par_label]['max'] < par.vals[ti]:
                            par.vals[ti] = settings.linkpar_specs[par_label]['max']



    def calculateOutputs(self, settings):
        '''
        Calculate outputs (called cascade characteristics in settings).
        These outputs must be calculated in the same order as defined in settings, otherwise references may break.
        Include any parameters marked as an output in the cascade sheet.
        '''

        outputs = odict()

        # For characteristics...
        for cid in settings.charac_specs:
            outputs[cid] = odict()
            for pop in self.pops:
                outputs[cid][pop.label] = None

                # Sum up all relevant compartment popsizes (or previously calculated characteristics).
                # TODO: Don't recalculate if characteristic was calculated as a dependency.
                for inc_label in settings.charac_specs[cid]['includes']:
                    if inc_label in self.getPop(pop.label).comp_ids:
                        vals = self.getPop(pop.label).getComp(inc_label).popsize
                    elif inc_label in outputs.keys()[:-1]:
                        vals = outputs[inc_label][pop.label]
                    else:
                        raise OptimaException('ERROR: Compartment or characteristic "%s" has not been pre-calculated for use in calculating "%s".' % (inc_label, cid))

                    if outputs[cid][pop.label] is None:
                        outputs[cid][pop.label] = dcp(vals)
                    else:
                        outputs[cid][pop.label] += vals

                # Divide by relevant compartment popsize (or previously calculated characteristic).
                if 'denom' in settings.charac_specs[cid]:
                    den_label = settings.charac_specs[cid]['denom']
                    if den_label in outputs.keys()[:-1]:
                        vals = outputs[den_label][pop.label]
                    elif den_label in self.getPop(pop.label).comp_ids:
                        vals = self.getPop(pop.label).getComp(den_label).popsize
                    else:
                        raise OptimaException('ERROR: Compartment or characteristic "%s" has not been pre-calculated for use in calculating "%s".' % (inc_label, cid))

                    outputs[cid][pop.label] /= (vals + project_settings.TOLERANCE)

        # For parameters...
        for cid in settings.linkpar_outputs:
            outputs[cid] = odict()
            for pop in self.pops:
                if 'tag' in settings.linkpar_specs[cid]:    # If a parameter is a transition, grab values from its (first, if duplicated,) link object.
                    tag = settings.linkpar_specs[cid]['tag']
                    outputs[cid][pop.label] = pop.getLinks(tag)[0].vals
                elif cid in pop.dep_ids:             # If a parameter is a dependency, grab values from the dependency object.
                    outputs[cid][pop.label] = pop.getDep(cid).vals
                else:
                    raise OptimaException('ERROR: Output parameter "%s" was not calculated as a "dependency" during the model run.' % cid)

        return outputs

    def printModelState(self, ti):
        for pop in self.pops:
            logging.info("Population: %s" % pop.label)
            logging.info(pop.getModelState(ti))


    def isBirthCompartment(self, comp_id, settings):
        """
        Returns True if comp_label is label for a compartment tagged as births. 
        """
        if 'tag_birth' in settings.node_specs[comp_id]:
            return True
        return False

    def isDeadCompartment(self, comp_id, settings):
        """
        Returns True if comp_label is label for a compartment tagged as dead. 
    
        """
        if 'tag_dead' in settings.node_specs[comp_id]:
            return True
        return False


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


