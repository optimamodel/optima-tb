import logging
logger = logging.getLogger(__name__)

from optima_tb.utils import OptimaException, odict, defaultrepr, objrepr

import numpy as np
from math import ceil, floor
from uuid import uuid4 as uuid

#%% Resultset class that contains one set of results
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
    
    def __init__(self, model, parset, settings, name=None):
        
        self.uuid = uuid()
        if name is None:
            self.name = 'results:' + parset.name
        else:
            self.name = name
        self.parset_name = parset.name
        self.parset_id  = parset.uid
        
        self.dt = settings.tvec_dt
        self.t_step = model.sim_settings['tvec']
        self.indices_observed_data = np.where(self.t_step%1.0==0)
        self.t_observed_data = self.t_step[self.indices_observed_data]
        
        self.outputs = model.calculateOutputs(settings = settings)
        
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
        self.comp_label_names = self.__generateLabelNames(self.comp_specs.keys(),self.comp_labels)
        self.char_labels = self.outputs.keys() # definitely need a better way of determining these
        # /work-in-progress
    
    def __repr__(self):
        ''' Print out useful information when called'''
        output = '====================================================================================\n'
        output += objrepr(self)
        return output
        
        
    def extractSimulationData(self):
        """
        
        Currently, only extract the simulated data for characteristics, as we should only use
        this in fitting for calibration.
        """
        pass
        
    def __generateLabelNames(self,labels,names):
        """
        
        Note: this could potentially go in utils or another location, as it's not really 
        specific to results.py. 
        """
        label_names = odict()
        for (i,label) in enumerate(labels):
            label_names[label] = names[i]
        return label_names
        
        
        
    def getCompartmentSizes(self, pop_labels = None, comp_label = None, use_observed_times=False):
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
            
            for cj in range(len(self.comp_specs.keys())):#comps:
                # TODO: this is an inefficient implementation and needs to be improved; possibly with a rethink of the data structures
                if not self.m_pops[p_index].comps[cj].label in comps:
                    continue
        
                comp_labels.append(self.m_pops[p_index].comps[cj].label)
        
                if use_observed_times:
                    datapoints[pi][self.m_pops[p_index].comps[cj].label] = self.m_pops[p_index].comps[cj][self.indices_observed_data]
                else:
                    datapoints[pi][self.m_pops[p_index].comps[cj].label] = self.m_pops[p_index].comps[cj]
                
        
        return datapoints, pops, comp_labels
    
    
    def getCharacteristicDatapoints(self,pop_label = None, char_label = None, use_observed_times=False):
        """
        Returns the data points for the simulation, for each characteristics, 
        that should correspond to times of the observed data. This method is intended for
        use with calibration. 
        
        [Intended]
        If pop_id or char_id are specified, then the method returns the simulated data points for
        just those populations or compartments.
        Otherwise, the simulated data points for all popuations and compartments are returned.
        
        @param pop_id: population id, either as single id or as list
        @param char_id: characteristic id, either as single label or as list of labels
        
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
      
      
    def getFlow(self, comp, inflows = True, outflows = True, invert = False, link_legend = None, exclude_transfers = False):
        """
        TODO finish method (moved across from plotting.py)
        
        """
        all_labels = []
        all_rates = []
        all_tvecs = []
        
        for in_out in xrange(2): # in, then out
            if (in_out == 0 and plot_inflows) or (in_out == 1 and plot_outflows):
                comp_link_ids = [comp.inlink_ids, comp.outlink_ids][in_out]
                label_tag = ['In: ','Out: '][in_out]
                for link_tuple in comp_link_ids:
                    link = self.m_pops[link_tuple[0]].links[link_tuple[1]]
                    #print link.label, link_labels, include_link_not_exclude
                    if link_labels is None or (invert and link.label in link_labels) or (not invert and link.label not in link_labels):
                        try: 
                            legend_label = label_tag + settings.linkpar_specs[link.label]['name']
                        except: 
                            if exclude_transfers: continue
                            else: legend_label = label_tag + link.label
                        if link.label in link_legend:
                            legend_label = link_legend[link.label]    # Overwrite legend for a label.
                        num_flow = dcp(link.vals)
                        if in_out == 0:
                            comp_source = self.m_pops[link.index_from[0]].comps[link.index_from[1]]
                        else:
                            comp_source = comp
                        was_proportion = False
                        if link.val_format == 'proportion':
                            denom_val = sum(self.m_pops[lid_tuple[0]].links[lid_tuple[-1]].vals for lid_tuple in comp_source.outlink_ids)
                            num_flow /= denom_val
                            was_proportion = True
                        if link.val_format == 'fraction' or was_proportion is True:
                            if was_proportion is True:
                                num_flow *= comp_source.popsize_old
                            else:
                                num_flow = 1 - (1 - num_flow) ** self.dt     # Fractions must be converted to effective timestep rates.
                                num_flow *= comp_source.popsize
                            num_flow /= self.dt      # All timestep-based effective fractional rates must be annualised.
                            
                        all_labels.append(legend_label)
                        all_rates.append(num_flow)
                        all_tvecs.append(tvec)
        
        return all_rates, all_labels
        
                
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
            output = sep.join(['Indicator','Population'] + ['%g'%t for t in self.t_step]) # Create header and years 
            npts = len(self.t_step)
        else:
            output = sep.join(['Indicator','Population'] + ['%g'%t for t in self.t_observed_data]) # Create header and years
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
                
                output += key+sep+popkey+sep
                for t in range(npts):
                    output += ('%g'+sep) % data[t]
            
        if writetofile: 
            with open(filename, 'w') as f: f.write(output)
            logger.info('Results exported to "%s"' % (filename))
            return None
        else:
            return output
        
        
