import logging
logger = logging.getLogger(__name__)

from optima_tb.utils import OptimaException, odict, defaultrepr, objrepr
import optima_tb.settings as project_settings

import numpy as np
from math import ceil, floor
from uuid import uuid4 as uuid
from copy import deepcopy as dcp
from collections import defaultdict

# %% Resultset class that contains one set of results
class ResultSet(object):
    """
    Class to hold one set of Results
    Fields: 
            name               name of result
            parset_name        name (not index) of parset used to create this result set
            parset_id          uuid of relevant parset (in case of duplication or missing parset_name)
            sim_settings       settings for simulation
            pop_labels         population tags for which siumulation was run (e.g. [])
            char_label         characteristics tags for which simulation was run
            ind_observed_data  indices for t that correspond to times at which observed data was collected
                               so that the observed data points can be compared against simulated data
                               using data_observed - data_sim[indices_observed_data]
            t_observed_data    t values for indices_observed_data points
                               popdata_sim        simulated compartment data points (format: ?)
            t_step             t_steps in real time 
            dt                 dt used to create this resultset
            model              simulated model object
            
        Optional: ------------------------------------------------------
            data_observed      datapoints supplied (format: dataset)
            calibration_fit    calibration fit score
            calibration_metric calibration metric used
    
        Example
            if dt = 0.25, and observed data was taken annually
            then 
                indices_observed_data = [0,4,8,..]. This can be determined directly from dt
                t_observed_data       = [2000.,2001.,2002., ...]
    """

    def __init__(self, model, parset, settings, progset=None, budget_options=None, name=None):

        self.uuid = uuid()

        if name is None:
            self.name = parset.name + ('-%s' % (progset.name) if progset is not None else '')
        else:
            self.name = name

        self.parset_name = parset.name
        self.parset_id = parset.uid

        self.dt = settings.tvec_dt
        self.t_step = model.t
        self.indices_observed_data = np.where(self.t_step % 1.0 == 0)
        self.t_observed_data = self.t_step[self.indices_observed_data]

        # Set up for future use
        self.calibration_fit = None
        self.calibration_score = None

        """
        # remaining fields to be set:
        self.popdata_sim
        self.chardata_sim
        self.data_observed
        """

        # work-in-progress: in time, these sections should be removed and only the data
        # points we are interested should be listed
        self.model = model # Reference to original model - each model is only wrapped in one result, so no need to copy it

        self.pop_label_index = {}
        for i, pop in enumerate(self.model.pops):
            self.pop_label_index[pop.label] = i

        self.sim_settings = model.sim_settings
        self.pop_labels = [pop.label for pop in model.pops]
        self.comp_labels = settings.node_names
        self.comp_specs = settings.node_specs
        self.comp_label_names = self.__generateLabelNames(self.comp_specs.keys(), self.comp_labels)
        self.char_labels = [c.label for c in self.model.pops[0].characs] # Double check if output parameters are included here too?
        self.link_labels = []
        for pop in model.pops:
            for link in pop.links:
                if link.label not in self.link_labels:
                    self.link_labels.append(link.label)

        self.budgets = {} # placeholders
        self.coverages = {}
        if progset is not None and budget_options is not None:
            budget_options = dcp(budget_options)
#            print [p.label for p in progset.progs]
            # we have results for progsets and budget_options
            if budget_options.has_key('alloc_is_coverage') and budget_options['alloc_is_coverage']: # TODO update post-Belarus
                self.coverages = budget_options['init_alloc']
                self.budgets = progset.getBudgets(self.coverages)
            else:
                self.budgets = budget_options['init_alloc']
                self.coverages = progset.getCoverages(self.budgets)

        # /work-in-progress

    def __repr__(self):
        ''' Print out useful information when called'''
        output = '====================================================================================\n'
        output += objrepr(self)
        return output


#     def extractSimulationData(self):
#         """
#
#         Currently, only extract the simulated data for characteristics, as we should only use
#         this in fitting for calibration.
#         """
#         pass

    def __generateLabelNames(self, labels, names):
        """
        
        Note: this could potentially go in utils or another location, as it's not really 
        specific to results.py. 
        """
        label_names = odict()
        for (i, label) in enumerate(labels):
            label_names[label] = names[i]
        return label_names

    def getValuesAt(self, label, year_init, year_end=None, pop_labels=None, integrated=False):
        """
        Derives transition flow rates, characteristic or compartment values for results, according to label provided.
       
        The range of values correspond to timepoints from year_init (inclusive) to year_end (exclusive).
        In the absence of a year_end, the value corresponding to only the year_init timepoint is returned.
        Values are summed across all population groups unless a subset is specified with pop_labels.
        Values can also be optionally integrated to return a scalar.
        
        Outputs values as array or integrated scalar, as well as the the corresponding timepoints.
        """

        dt = self.dt

        if pop_labels is None:
            pop_labels = self.pop_labels

        # Get time array ids for values between initial (inclusive) and end year (exclusive).
        tvals = np.array(self.model.t)
        if year_end is None: year_end = year_init + dt
        idx = (tvals >= year_init - dt / 2) * (tvals <= year_end - dt / 2)

        # Initialise output as appropriate array of zeros.
        output = np.zeros(len(tvals[idx]))
        units = ""

        # Find values in results and add them to the output array per relevant population group.
        if label in self.model.pops[0].link_lookup:
            values = self.getFlow(label, pop_labels=pop_labels)[0]
            for pop in values.keys():
                popvalues = values[pop]
                output += popvalues[idx]
        else:
            for pop in pop_labels:
                var = self.model.getPop(pop).getVariable(label)[0]
                output += var.vals[idx]

        # Do a simple integration process if specified by user.
        if integrated:
            output = output.sum() * dt

        return output, tvals[idx], units

    def getCompartmentSizes(self, pop_labels=None, comp_label=None, use_observed_times=False):
        """
        #PopDatapoints
        
        Returns the data points for the simulation, for each compartment, 
        that should correspond to times of the observed data. This method is intended for
        use with calibration. 
        
        [Comment to DK: will we need compartments at all during calibration, scenarios or optimization?
        Possibly for post-hoc analysis, but perhaps it would be worthwhile to assume no, and only
        for some cases should a user - perhaps via a flag in settings - save all data.]
        
        [Intended]
        If pop_id or comp_id are specified, then the method returns the simulated data points for
        just those populations or compartments.
        Otherwise, the simulated data points for all populations and compartments are returned.

        @param pop_id: population id, either as single id or as list
        @param comp_id: compartment id, either as single id or as list
                
        """
        if pop_labels is not None:
            if isinstance(pop_labels, list):
                pops = pop_labels
            else:
                pops = [pop_labels]
        else:
            pops = self.pop_labels

        if comp_label is not None:
            if isinstance(comp_label, list):
                comps = comp_label
            elif isinstance(comp_label, odict):
                comps = comp_label.keys()
            else:
                comps = [comp_label]
        else:
            comps = self.comp_specs.keys()

        # this currently uses odict for the output ...
        datapoints = odict()
        comp_labels = []
        for pi in pops:

            datapoints[pi] = odict()
            p_index = self.pop_label_index[pi]

            for ci, c_label in enumerate(comps):

                comp = self.model.pops[p_index].getComp(c_label)
                if use_observed_times:
                    datapoints[pi][c_label] = comp[self.indices_observed_data]
                else:
                    datapoints[pi][c_label] = comp
                comp_labels.append(c_label)

        units = 'people'

        return datapoints, pops, comps, units


    def getCharacteristicDatapoints(self, pop_label=None, char_label=None, use_observed_times=False):
        """
        Returns the data points for the simulation, for each characteristics, 
        that should correspond to times of the observed data. This method is intended for
        use with calibration. 
        
        [Intended]
        If pop_id or char_id are specified, then the method returns the simulated data points for
        just those populations or compartments.
        Otherwise, the simulated data points for all popuations and compartments are returned.
        
        Param:
             pop_label                          population labels, either as single id or as list
             char_label                         characteristic id, either as single label or as list of labels
             use_observed_times                 boolean flag. If True, returns only datapoints for which there was a corresponding data observation; 
                                                else, returns all characteristic data points
        
        """
        if pop_label is not None:
            if isinstance(pop_label, list):
                pop_labels = pop_label
            else:
                pop_labels = [pop_label]
        else:
            pop_labels = self.pop_labels

        if char_label is not None:
            if isinstance(char_label, list):
                char_labels = char_label
            else:
                char_labels = [char_label]
        else:
            char_labels = self.char_labels

        datapoints = defaultdict(dict)

        for pop_label in pop_labels:
            for char_label in char_labels:
                if use_observed_times:
                    datapoints[char_label][pop_label] = self.model.getPop(pop_label).getVariable(char_label)[0].vals[self.indices_observed_data]
                else:
                    datapoints[char_label][pop_label] = self.model.getPop(pop_label).getVariable(char_label)[0].vals

        units = ''

        return datapoints, char_label, pop_label, units


    def getFlow(self, link_tag, pop_labels=None,annualize=True,as_fraction=False):
        """
        Return the flow at each time point in the simulation for a single parameter

        INPUTS
        - par_label : A string specifying a single Parameter to retrieve flow rates for
        - pop_label : A list or list of strings of population labels. If None, use all populations
        - annualize : Boolean which specifies if the number of moved people should be annualized or not. If True, an annual average is computed; if False, the number of people per time step is computed
        - as_fraction : Boolean which specifies if the flow rate should be expressed as a fraction of the source compartment size.
            - If True, the fractional flow rate for each link will be computed by dividing the net flow by the sum of source compartment sizes
            - If None, it will be set to the same value as par.is_fraction for each parameter requested
        For each parameter label, the flow from all links deriving from the parameter will be summed within requested populations. Thus if there are multiple links
        with the same link_label (e.g. a 'doth' link for death by other causes, for which there may be one for every compartment)
        then these will be aggregated

        OUTPUT
        - datapoints : dictionary `datapoints[pop_label]` has an array of flow rate in requested units
        - units : the units of the returned data
        
        Note that the flow rate is linked to the time step by

        popsize[t+1] = popsize[t]+net_flow[t]

        That is, it is the number of people who will transition between t and t+dt (i.e. after the current time)
        """

        if pop_labels is not None:
            if isinstance(pop_labels, list):
                pop_labels = pop_labels
            else:
                pop_labels = [pop_labels]

        datapoints = defaultdict(float)
        source_size = defaultdict(float)

        for pop in self.model.pops:
            if pop_labels is None or pop.label in pop_labels:
                for link in pop.getLinks(link_tag):
                    if link.vals is None:
                        raise OptimaException('Requested flow rate "%s" was not recorded because only partial results were saved' % (link.label))

                    datapoints[pop.label] += link.vals
                    source_size[pop.label] += (link.source.vals if not link.source.is_junction else sum([x.vals for x in link.source.outlinks]))

        # If as_fraction is None, use the same units as the Parameter. All Parameters should have the same units
        # in all populations so can use whichever one is left after the loop above
        if as_fraction is None and (par.units == 'fraction' or par.units == 'proportion'):
            as_fraction = True

        if as_fraction:
            units = 'proportion'
            for pop_label in datapoints:
                datapoints[pop_label] /=  source_size[pop_label]
        else:
            units = 'people'

        # If we need to convert from dt units to annualized units
        # If as_fraction is true, then the quantity is a proportion and no further time units are required
        if annualize and not as_fraction:
            units += '/year'
            for pop_label in datapoints:
                datapoints[pop_label] = datapoints[pop_label] * (1.0 / self.dt)
        elif not as_fraction:
            units += '/timestep'

        return datapoints, units


    def export(self, filestem=None, sep=',', writetofile=True, use_alltimesteps=True):
        """
        Export method for characteristics results obtained from a simulation that should correspond 
        to times of the observed data (i.e. annually). This method is intended for use with runSim 
        currently and will be extended to include optimization and scenario results.
        """
        import os

        if filestem is None:  # Doesn't include extension, hence filestem
            filestem = self.name
        filestem = os.path.abspath(filestem)
        filename = filestem + '.csv'

        # All dumpable items
        keys = [c.label for c in self.model.pops[0].characs] + [p.label for p in self.model.pops[0].pars]

        if use_alltimesteps:
            output = sep.join(['Indicator', 'Population'] + ['%g' % t for t in self.t_step]) # Create header and years
            npts = len(self.t_step)
        else:
            output = sep.join(['Indicator', 'Population'] + ['%g' % t for t in self.t_observed_data]) # Create header and years
            npts = len(self.t_observed_data)

        for key in keys:
            output += '\n' # Add a line break between different indicators
            for pop in self.pop_labels:
                output += '\n'
                try:
                    var = self.model.getPop(pop).getVariable(key)[0]
                except OptimaException as e:
                    if str(e).startswith("Object not found"):
                        continue
                    else:
                        raise

                if use_alltimesteps:
                    data = var.vals
                else:
                    data = var.vals[self.indices_observed_data]

                output += key + sep + pop + sep

                if data is None: # If full_output=False then some fields will be skipped
                    output += 'Not stored - set full_output to True to retain'
                else:
                    for t in range(npts):
                        output += ('%g' + sep) % data[t]

        if writetofile:
            with open(filename, 'w') as f: f.write(output)
            logger.info('Results exported to "%s"' % (filename))
            return None
        else:
            return output


