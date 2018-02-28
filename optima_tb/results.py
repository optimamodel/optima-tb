import logging
logger = logging.getLogger(__name__)

from optima_tb.utils import OptimaException, odict, defaultrepr, objrepr
import optima_tb.settings as project_settings

import numpy as np
from math import ceil, floor
from uuid import uuid4 as uuid
from copy import deepcopy as dcp

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
            outputs            simulated characteristic data points (format: ?)
            m_pops             totals per pop per compartment
            
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
            self.name = 'results:' + parset.name
        else:
            self.name = name
        self.parset_name = parset.name
        self.parset_id = parset.uid

        self.dt = settings.tvec_dt
        self.t_step = model.sim_settings['tvec']
        self.indices_observed_data = np.where(self.t_step % 1.0 == 0)
        self.t_observed_data = self.t_step[self.indices_observed_data]

        self.outputs = model.calculateOutputs(settings=settings)

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
        self.m_pops = model.pops

        self.pop_label_index = {}
        for i, pop in enumerate(self.m_pops):
            self.pop_label_index[pop.label] = i

        self.sim_settings = model.sim_settings
        self.pop_labels = self.outputs[0].keys()
        self.comp_labels = settings.node_names
        self.comp_specs = settings.node_specs
        self.comp_label_names = self.__generateLabelNames(self.comp_specs.keys(), self.comp_labels)
        self.char_labels = self.outputs.keys() # definitely need a better way of determining these
        self.link_labels = [link.label for link in model.pops[0].links]
        self.link_label_ids = {label:[i for i, link in enumerate(model.pops[0].links) if link.label == label] for label in self.link_labels} # For each link label, the index of the link list that had it

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
        tvals = np.array(self.sim_settings['tvec'])
        if year_end is None: year_end = year_init + dt
        idx = (tvals >= year_init - dt / 2) * (tvals <= year_end - dt / 2)

        # Initialise output as appropriate array of zeros.
        output = np.zeros(len(tvals[idx]))

        # Find values in results and add them to the output array per relevant population group.
        # TODO: Semantics need to be cleaned during design review phase.
        #       Link tags are actually stored in link_labels, while link labels could be in char_labels if a transition is marked as a result output.
        if label in self.link_labels:

            values, _, _, units = self.getFlow(link_label=label, pop_labels=pop_labels)  # Does not return link values directly but calculates flows instead.
            values = values[label]

            for pop in values.keys():
                popvalues = values[pop]
                output += popvalues[idx]

        elif label in self.char_labels:

            values, _, _, units = self.getCharacteristicDatapoints(char_label=label, pop_label=pop_labels, use_observed_times=False)
            values = values[label]

            for pop in values.keys():
                popvalues = values[pop]
                output += popvalues[idx]

        elif label in self.comp_label_names.keys():

            values, _, _, units = self.getCompartmentSizes(comp_label=label, pop_labels=pop_labels, use_observed_times=False)
            for pop in values.keys():
                popvalues = values[pop]
                output += values[pop][label].popsize[idx]

        else:
            logger.warn('Unable to find values for label="%s", with no corresponding characteristic, transition or compartment found.' % label)

        # Do a simple integration process if specified by user.
        if integrated:
            output = output.sum() * dt

        return output, tvals[idx], units


    def getValueAt(self, label, year_init, year_end=None, pop_labels=None):
        """
        Returns scalar at or over a specific period for the value of a label corresponding to either a:
            * characteristic
            * flow
            * compartment
        If over a time period, the values are summed to form the cumulative. The scalar of this value is 
        then returned.
        
        Params:
            label        string value of label
            year_init    time period at which to return the value, or if year_end is specified, the beginning of the time period
            year_end     (optional) 
            pop_labels   (optional) list of populations to evaluate over. Default: None, which evaluates over all
        
        Return:
            value        a scalar
        """
        logging.warn("Method deprecated. Please use result.getValuesAt() instead")
        dt = self.dt

        if pop_labels is None:
            pop_labels = self.pop_labels

        values = []

        tvals = np.array(self.sim_settings['tvec'])
        if year_end is not None: # including the dt/2 window to guard against float comparison
            idx = (tvals >= year_init - dt / 2) * (tvals <= year_end + dt / 2)
        else:
            idx = (tvals >= year_init - dt / 2) * (tvals <= year_init + dt / 2)

        val = 0


        if label in self.char_labels:
#             print "is Char"

            values, _, _ = self.getCharacteristicDatapoints(char_label=label, pop_label=pop_labels, use_observed_times=False)
            values = values[label]

            for pop in values.keys():
                popvalues = values[pop]
                val += popvalues[idx].sum()


        elif label in self.link_labels:
#             print "is Link"

            values, _, _ = self.getFlow(link_label=label, pop_labels=pop_labels)
            values = values[label]

            for pop in values.keys():
                popvalues = values[pop]
                val += popvalues[idx].sum()

            if year_end is not None:
                # getFlow returns an annualised rate, which is correct if we're looking at a  value
                # for a point in time. If we're looking at the value over a range, then we should multiply by
                # dt in order to give the correct value over that period.
                val *= dt


        elif label in self.comp_label_names.keys():
#             print "is Comp"

            values, _, _ = self.getCompartmentSizes(comp_label=label, pop_labels=pop_labels, use_observed_times=False)
            for pop in values.keys():
                popvalues = values[pop]
                pp = popvalues[label]
                val += values[pop][label].popsize[idx].sum()

        else:
            logger.warn("Unable to find value for label='%s': no corresponding characteristic, link, or compartment found." % label)

        return val






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

                comp = self.m_pops[p_index].getComp(c_label)
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
                pops = pop_label
            else:
                pops = [pop_label]
        else:
            pops = self.pop_labels

        if char_label is not None:
            if isinstance(char_label, list):
                chars = char_label
            else:
                chars = [char_label]
        else:
            chars = self.char_labels


        # this currently uses odict for the output ...
        datapoints = odict()
        for cj in chars:
            datapoints[cj] = odict()
            for pi in pops:
                if use_observed_times:
                    datapoints[cj][pi] = self.outputs[cj][pi][self.indices_observed_data]
                else:
                    datapoints[cj][pi] = self.outputs[cj][pi]

        units = 'people'

        return datapoints, chars, pops, units


    def getFlow(self, link_label, pop_labels=None,target_flow=False,annualize=True):
        """
        Return the flow at each time point in the simulation

        INPUTS
        - link_label : A string or list of strings of link labels to get flows for. If None, use all links
        - pop_label : A list or list of strings of population labels. If None, use all populations
        - target_flow : By default, the actual flow rates accounting for compartment sizes during integration will be used. If target_flow=True, then the target flow rate will be returned
        - annualize : Boolean which specifies if the number of moved people should be annualized or not. If True, an annual average is computed; if False, the number of people per time step is computed
        
        For each link label, the flow will be summed within requested populations. Thus if there are multiple links
        with the same link_label (e.g. a 'doth' link for death by other causes, for which there may be one for every compartment)
        then these will be aggregated

        OUTPUT
        - datapoints : dictionary `datapoints[link_label][pop_label]` has an array of number of people moved at each point in time
        - link_label : links that were used
        - pops : Populations that were used
        
        Note that the flow rate is linked to the time step by

        popsize[t+1] = popsize[t]+net_flow[t]

        That is, it is the number of people who will transition between t and t+dt (i.e. after the current time)
        """

        if pop_labels is not None:
            if isinstance(pop_labels, list):
                pops = pop_labels
            else:
                pops = [pop_labels]
        else:
            pops = self.pop_labels

        if link_label is not None:
            if isinstance(link_label, list):
                pass
            else:
                link_label = [link_label]
        else:
            link_label = self.link_label

        datapoints = odict()

        for link_lab in link_label:

            datapoints[link_lab] = odict()
            for pi in pops:

                datapoints[link_lab] [pi] = odict()

                p_index = self.pop_label_index[pi]

                flows = []
                for link_index in self.link_label_ids[link_lab]:
                    link = self.m_pops[p_index].links[link_index]
                    if target_flow:
                        flows.append(link.target_flow)
                    else:
                        flows.append(link.flow)

                datapoints[link_lab][pi] = sum(flows)

                if annualize:
                    datapoints[link_lab][pi] /= self.dt

        if annualize:
            units = 'people/year'
        else:
            units = 'people/timestep'

        return datapoints, link_label, pops, units


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


        keys = self.char_labels


        if use_alltimesteps:
            output = sep.join(['Indicator', 'Population'] + ['%g' % t for t in self.t_step]) # Create header and years
            npts = len(self.t_step)
        else:
            output = sep.join(['Indicator', 'Population'] + ['%g' % t for t in self.t_observed_data]) # Create header and years
            npts = len(self.t_observed_data)
        for key in keys:
            output += '\n' # Add a line break between different indicators
            popkeys = self.pop_labels
            for pk, popkey in enumerate(popkeys):
                output += '\n'
                if use_alltimesteps:
                    data = self.outputs[key][popkey]
                else:
                    data = self.outputs[key][popkey][self.indices_observed_data]

                output += key + sep + popkey + sep
                for t in range(npts):
                    output += ('%g' + sep) % data[t]

        if writetofile:
            with open(filename, 'w') as f: f.write(output)
            logger.info('Results exported to "%s"' % (filename))
            return None
        else:
            return output


