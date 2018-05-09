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
        self.link_labels = model.pops[0].link_ids.keys()
        self.link_label_ids = model.pops[0].link_ids

        self.budgets = {} # placeholders
        self.coverages = {}
        if progset is not None and budget_options is not None:
            budget_options = dcp(budget_options)
#            print [p.label for p in progset.progs]
            # we have results for progsets and budget_options
            if budget_options.has_key('alloc_is_coverage') and budget_options['alloc_is_coverage']: # TODO update post-Belarus
                self.coverages = budget_options['init_alloc']
                self.budgets = progset.getBudgets(self.coverages, model)
            else:
                self.budgets = budget_options['init_alloc']
                self.coverages = progset.getCoverages(self.budgets, model)

        self.dw = dcp(parset.dw)
        self.years_lost = dcp(parset.years_lost)

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

    def getValuesAt(self, label, year_init, year_end=None, pop_labels=None, settings=None, integrated=False):
        """
        Derives transition flow rates, characteristic or compartment values for results, according to label provided.
       
        The range of values correspond to timepoints from year_init (inclusive) to year_end (exclusive).
        In the absence of a year_end, the value corresponding to only the year_init timepoint is returned.
        Values are summed across all population groups unless a subset is specified with pop_labels.
        Values can also be optionally integrated to return a scalar.

        settings is only required to determine DALYs and is ignored if label != 'dalys'

        Outputs values as array or integrated scalar, as well as the the corresponding timepoints.
        """

        dt = self.dt

        if pop_labels is None:
            pop_labels = self.pop_labels
        if not isinstance(pop_labels, list):
            pop_labels = [pop_labels]

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

            values = self.getFlow(link_label=label, pop_labels=pop_labels, annualize=False)  # Does not return link values directly but calculates flows instead.
            values = values[label]

            for pop in values.keys():
                popvalues = values[pop]
                output += popvalues[idx]

        elif label in self.char_labels:

            values, _, _ = self.getCharacteristicDatapoints(char_label=label, pop_label=pop_labels, use_observed_times=False)
            values = values[label]

            for pop in values.keys():
                popvalues = values[pop]
                output += popvalues[idx]

        elif label in self.comp_label_names.keys():

            values, _, _ = self.getCompartmentSizes(comp_label=label, pop_labels=pop_labels, use_observed_times=False)
            for pop in values.keys():
                popvalues = values[pop]
                output += values[pop][label].popsize[idx]

        elif label == 'daly':
            if settings is not None:
                output = self.getDALY(settings, year_init, year_end, pop_labels, integrated)
            else:
                raise OptimaException('Cannot calculate DALYs without a Settings object. ' +
                                      'Please provide Project.settings to ResultSet.getValuesAt()')
        else:
            logger.warn('Unable to find values for label="%s", with no corresponding characteristic, '.format(label) +
                        'transition or compartment found.')

        # Do a simple integration process if specified by user.
        if integrated:
            # enforced that return-type is always an iterable
            output = np.array([output.sum() * dt])

        return output, tvals[idx]


    def getValueAt(self, label, year_init, year_end=None, pop_labels=None):
        return self.getValuesAt(label, year_init, year_end, pop_labels, True)[0]


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

        return datapoints, pops, comps


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
        return datapoints, chars, pops


    def getDALY(self, settings, year_start, year_end=None, pop_labels=None, integrated=False):
        """
        Determine disability-adjusted life years (DALY) over a given interval. It is computed by summing up YLL and YLD.

        :param settings: Settings which contains compartment information
        :param year_start: float which specifies the start of the interval in which DALY is determined
        :param year_end: None or float. If float, it specifies the end of the interval in which YLL is determined; if
            it is None, the end of simulation is chosen as the end of interval
        :param pop_labels: None or list of str. If list of str, it contains the population labels for which the YLL is
            determined; if None, the YLL of all populations is determined
        :param integrated: if True, all values are summed up, otherwise values per time step are retained
        :return: an odict of {pop_label: DALY}
        """
        yll = self.getYLL(settings, year_start, year_end, pop_labels, integrated)
        yld = self.getYLD(settings, year_start, year_end, pop_labels, integrated)

        daly = odict()
        for pop in yll:
            daly[pop] = yll[pop] + yld[pop]

        return daly


    def getFlow(self, link_label=None, pop_labels=None, annualize=True):
        """
        Returns the number of people (as absolute number) moved from one compartment to another during the simulation.

        :param link_label: None or a list of str which represent from which links/transitions the number of moved people
            are extracted. In case of None, all are considered.
        :param pop_labels: None or a list of str which specify which populations to consider. In case of None, all are
            considered.
        :param annualize: Boolean which specifies if the number of moved people should be annualized or not. If True,
            an annual average is computed; if False, the number of people per time step is computed
        :return: dict of {link_label: {pop_label: value}}, where value is either a float (annualize==True)
            or a np.array(annualize==False)
        """
        if pop_labels is not None:
            if isinstance(pop_labels, list):
                pops = pop_labels
            else:
                pops = [pop_labels]
        else:
            pops = self.pop_labels

        if link_label is not None:
            if not isinstance(link_label, list):
                link_label = [link_label]
        else:
            link_label = self.link_labels

        datapoints = odict()

        for link_lab in link_label:
            datapoints[link_lab] = odict()

            for pi in pops:
                datapoints[link_lab][pi] = odict()
                p_index = self.pop_label_index[pi]
                num_flow_list = np.zeros(len(self.t_step))

                for link_index in self.link_label_ids[link_lab]:
                    link = self.m_pops[p_index].links[link_index]
                    comp_source = self.m_pops[p_index].comps[link.index_from[1]]
                    # scaled annual flow rate
                    num_flow = link.vals * link.scale_factor

                    was_proportion = False
                    if link.val_format == 'proportion':
                        denom_val = sum(self.m_pops[lid_tuple[0]].links[lid_tuple[-1]].vals
                                        for lid_tuple in comp_source.outlink_ids)
                        num_flow /= denom_val
                        was_proportion = True
                    if link.val_format == 'fraction' or was_proportion is True:
                        if was_proportion is True:
                            num_flow *= comp_source.popsize_old
                        else:
                            # remove TOLERANCE excess from flow values
                            num_flow[np.logical_and(num_flow > 1., num_flow < (1. + project_settings.TOLERANCE))] = 1.
                            num_flow[np.logical_and(num_flow < 0., num_flow > (0. - project_settings.TOLERANCE))] = 0.

                            # convert annual flow to flow per time step
                            num_flow_list = (1. - (1. - num_flow) ** self.dt) * comp_source.popsize

                if annualize:
                    # All timestep-based effective fractional rates must be annualised.
                    datapoints[link_lab][pi] = np.array([np.sum(num_flow_list / self.dt)])
                else:
                    datapoints[link_lab][pi] = num_flow_list

        return datapoints

    def getYLD(self, settings, year_start, year_end=None, pop_labels=None, integrated=False):
        """
        Determine years lost due to disability (YLD) over a given interval.

        :param settings: Settings which contains compartment information
        :param year_start: float which specifies the start of the interval in which YLL is determined
        :param year_end: None or float. If float, it specifies the end of the interval in which YLD is determined; if
            it is None, the end of simulation is chosen as the end of interval
        :param pop_labels: None or list of str. If list of str, it contains the population labels for which the YLL is
            determined; if None, the YLL of all populations is determined
        :param integrated: if True, all values are summed up, otherwise values per time step are retained
        :return: an odict of {pop_label: YLD}
        """
        if len(self.years_lost) == 0:
            raise OptimaException('ERROR: Cannot compute YLD because no ' +
                                  'disability weight is defined in the databook')

        # fix pop_labels: either use all pops available or only the ones passed to the function
        if pop_labels is None:
            pop_labels = self.pop_labels
        elif pop_labels is not None and isinstance(pop_labels, basestring):
            pop_labels = [pop_labels]

        # fix end point of evaluation
        if year_end is None:
            year_end = self.t_step[-1]

        infected_comps = filter(lambda x: 'infected' in settings.node_specs[x] and settings.node_specs[x]['infected'],
                                settings.node_specs)

        yld = odict()
        # for each population, compute how many people are infected in the specified interval
        for pop in pop_labels:
            yld[pop] = np.zeros(int(year_end - year_start))
            for comp in infected_comps:
                infected, _ = self.getValuesAt(comp, year_start, year_end, pop, None, integrated)
                yld[pop] += infected
            yld[pop] *= self.dw[pop]

        return yld

    def getYLL(self, settings, year_start, year_end=None, pop_labels=None, integrated=False):
        """
        Determine years of life lost due to premature mortality (YLL) over a given interval.

        :param settings: Settings which contains compartment information
        :param year_start: float which specifies the start of the interval in which YLL is determined
        :param year_end: None or float. If float, it specifies the end of the interval in which YLL is determined; if
            it is None, the end of simulation is chosen as the end of interval
        :param pop_labels: None or list of str. If list of str, it contains the population labels for which the YLL is
            determined; if None, the YLL of all populations is determined
        :param integrated: if True, all values are summed up, otherwise values per time step are retained
        :return: an odict of {pop_label: YLL}
        """
        if len(self.years_lost) == 0:
            raise OptimaException('ERROR: Cannot compute YLL because no ' +
                                  'average life expectancy is defined in the databook')

        # fix pop_labels: either use all pops available or only the ones passed to the function
        if pop_labels is None:
            pop_labels = self.pop_labels
        elif pop_labels is not None and isinstance(pop_labels, basestring):
            pop_labels = [pop_labels]

        # fix end point of evaluation
        if year_end is None:
            year_end = self.t_step[-1]

        # obtain relevant death compartments
        # introduced a new tag 'death by disease' which indicates yll-relevant deaths
        death_comps = filter(
            lambda x: 'dbd' in settings.node_specs[x] and settings.node_specs[x]['dbd'], settings.node_specs)
        # obtain transitions leading to the relevant death compartments
        death_trans = filter(lambda x: any([item[1] in death_comps for item in settings.links[x]]), settings.links)

        yll = odict()
        # for each population, compute how many people have died in the specified interval
        for pop in pop_labels:
            yll[pop] = np.zeros(int(year_end - year_start))
            for trans in death_trans:
                deaths, _ = self.getValuesAt(trans, year_start, year_end, pop, None, integrated)
                yll[pop] += deaths
            yll[pop] *= self.years_lost[pop]

        return yll


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


