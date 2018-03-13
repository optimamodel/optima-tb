# %% Imports
import logging
from matplotlib.pyplot import plot

logger = logging.getLogger(__name__)

from optima_tb.utils import odict, OptimaException
from optima_tb.results import ResultSet
import numpy as np

import pylab as pl
from copy import deepcopy as dcp
from copy import copy as ndcp
import matplotlib.cm as cmx
import matplotlib.colors as colors
from random import shuffle
import numbers
import textwrap
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib

def plotSeries(proj,results,outputs,pops=None,axis='outputs',output_aggregation='sum',pop_aggregation='sum',plot_type='line',use_full_labels=True):
    # This function plots a time series for a model output quantities
    #
    # INPUTS
    # - proj - project object
    # - results - A ResultSet or a dict of ResultSets with key corresponding
    #   to name
    # - outputs - The name of an output compartment, characteristic, or
    #   parameter, or list of names. Inside a list, a dict can be given to
    #   specify an aggregation e.g. outputs=['sus',{'total':['sus','vac']}]
    #   where the key is the new name
    # - pops - The name of an output population, or list of names. Like
    #   outputs, can specify a dict with a list of pops to aggregate over them
    # - axis - Display one of results, outputs, or pops as different coloured
    #   lines on the plot. A new figure will be generated for all other
    #   combinations of the remaining quantities
    # - output_aggregation - If an output aggregation is requested, combine
    #   the outputs listed using one of
    #       - 'sum' - just add values together
    #       - 'average' - unweighted average of quantities
    #       - 'weighted' - weighted average where the weight is the
    #         compartment size, characteristic value, or link source
    #         compartment size (summed over duplicate links). 'weighted'
    #         method cannot be used with non-transition parameters and a
    #         KeyError will result in that case
    # - pop_aggregation - Same as output_aggregation, except that 'weighted'
    #   uses population sizes. Note that output aggregation is performed
    #   before population aggregation. This also means that population
    #   aggregation can be used to combine already aggregated outputs (e.g.
    #   can first sum 'sus'+'vac' within populations, and then take weighted
    #   average across populations)
    # - plot_type - 'line' or 'stacked'
    # - use_full_labels - Attempt to convert population and result labels into
    #   full names. If no matches are found, the provided name will be used
    #   directly


    if isinstance(results,ResultSet):
        results = {results.name:results}

    if pops is None:
        pops = results[results.keys()[0]].pop_labels

    if isinstance(outputs,str):
        outputs = [outputs]

    assert isinstance(pops, list), 'Populations need to be specified as a list'
    assert isinstance(outputs,list), 'Outputs need to be specified as a list or a string'

    def extract_labels(l):
        # Flatten the input arrays to extract all requested pops and outputs
        # e.g. ['vac',{'a':['vac','sus']}] -> ['vac','vac','sus'] -> set(['vac','sus'])
        # Also returns the labels for legends etc. i.e. ['vac','a']
        out = []
        for x in l:
            if isinstance(x,dict):
                k = x.keys()
                assert len(k) == 1, 'Aggregation dict can only have one key'
                out += x[k[0]]
            else:
                out.append(x)
        return set(out)

    # First, get all of the pops and outputs requested by flattening the lists
    pops_required = extract_labels(pops)
    outputs_required = extract_labels(outputs)

    final_outputs = dict() 
    tvecs = dict()

    for result_label,result in results.items(): # For each result

        final_outputs[result_label] = defaultdict(dict) # This final_outputs[result_label][agg_pop_label][agg_output_label] - i.e. the final data for plotting
        tvecs[result_label] = result.model.sim_settings['tvec']
        dt = result.model.sim_settings['tvec_dt']

        aggregated_outputs = defaultdict(dict) # Dict with aggregated_outputs[pop_label][aggregated_output_label]
        compsize = dict()
        popsize = dict()

        # Assemble the final output dicts for each population
        for pop_label in pops_required:
            pop = result.model.getPop(pop_label)
            popsize[pop_label] = pop.popsize()
            data_dict = dict()  # Temporary storage for raw outputs

            # First pass, extract the original output quantities
            for output_label in outputs_required:
                if output_label in pop.comp_ids:
                    data_dict[output_label] = pop.getComp(output_label).vals
                    compsize[output_label] = data_dict[output_label]
                elif output_label in pop.charac_ids:
                    data_dict[output_label] = pop.getCharac(output_label).vals
                    compsize[output_label] = data_dict[output_label]
                elif output_label in pop.par_ids:
                    par = pop.getPar(output_label)
                    if par.links: # If this is a transition parameter, use getFlow to get the flow rate
                        data_dict[output_label] = np.zeros(tvecs[result_label].shape)
                        compsize[output_label] = np.zeros(tvecs[result_label].shape)
                        for link in par.links:
                            data_dict[output_label] += link.vals/dt
                            compsize[output_label] += (link.source.vals if not link.source.is_junction else link.source.vals_old)

            # Second pass, aggregate them according to any aggregations present
            for output in outputs: # For each final output
                if isinstance(output,dict):
                    output_name = output.keys()[0]
                    labels = output[output_name]
                    if output_aggregation == 'sum': 
                        aggregated_outputs[pop_label][output_name] = sum(data_dict[x] for x in labels) # Add together all the outputs
                    elif output_aggregation == 'average': 
                        aggregated_outputs[pop_label][output_name] = sum(data_dict[x] for x in labels) # Add together all the outputs
                        aggregated_outputs[pop_label][output_name] /= len(labels)
                    elif output_aggregation == 'weighted':
                        aggregated_outputs[pop_label][output_name] = sum(data_dict[x]*compsize[x] for x in labels) # Add together all the outputs
                        aggregated_outputs[pop_label][output_name] /= sum([compsize[x] for x in labels])
                    else:
                        raise OptimaException('Unknown output aggregation type')
                else:
                    aggregated_outputs[pop_label][output] = data_dict[output]

        # Now aggregate over populations
        # If we have requested a reduction over populations, this is done for every output present
        for pop in pops: # This is looping over the population entries
            for output_name in aggregated_outputs[aggregated_outputs.keys()[0]].keys():
                if isinstance(pop,dict):
                    pop_name = pop.keys()[0]
                    pop_labels = pop[pop_name]
                    if pop_aggregation == 'sum': 
                        final_outputs[result_label][pop_name][output_name] = sum(aggregated_outputs[x][output_name] for x in pop_labels) # Add together all the outputs
                    elif pop_aggregation == 'average': 
                        final_outputs[result_label][pop_name][output_name] = sum(aggregated_outputs[x][output_name] for x in pop_labels) # Add together all the outputs
                        final_outputs[result_label][pop_name][output_name] /= len(labels)
                    elif pop_aggregation == 'weighted':
                        final_outputs[result_label][pop_name][output_name] = sum(aggregated_outputs[x][output_name]*popsize[x] for x in pop_labels) # Add together all the outputs
                        final_outputs[result_label][pop_name][output_name] /= sum([popsize[x] for x in pop_labels])
                    else:
                        raise OptimaException('Unknown output aggregation type')
                else:
                    final_outputs[result_label][pop][output_name] = aggregated_outputs[pop][output_name]


    # Now, we have completed all aggregations. The final_outputs dict is now ready to be plotted across the required axis
    final_result_labels = results.keys()
    final_pop_labels = final_outputs[final_result_labels[0]].keys()
    final_output_labels = final_outputs[final_result_labels[0]][final_pop_labels[0]].keys()

    figs = []

    if use_full_labels:
        name = getFullName
    else:
        name = lambda x,y: x
    
    if axis == 'results':
        for pop in final_pop_labels:
            for output in final_output_labels:
                figs.append(plt.figure())
                plt.xlabel('Year')
                plt.autoscale(enable=True, axis='x', tight=True)
                plt.ylabel(name(output,proj))
                plt.title('%s-%s' % (pop,output))
                if plot_type == 'stacked':
                    y = [final_outputs[result][pop][output] for result in final_result_labels]
                    labels = [name(result,proj) for result in final_result_labels]
                    plt.stackplot(tvecs[result],y,labels=labels)
                    plt.legend(labelspacing=-matplotlib.rcParams["legend.labelspacing"]-2)
                else:
                    for result in final_result_labels:
                        plt.plot(tvecs[result],final_outputs[result][pop][output],label=name(result,proj))
                    plt.legend()

    elif axis == 'pops':
        for result in final_result_labels:
            for output in final_output_labels:
                figs.append(plt.figure())
                plt.xlabel('Year')
                plt.autoscale(enable=True, axis='x', tight=True)
                plt.ylabel(name(output,proj))
                plt.title('%s-%s' % (result,output))
                if plot_type == 'stacked':
                    y = [final_outputs[result][pop][output] for pop in final_pop_labels]
                    labels = [name(pop,proj) for pop in final_pop_labels]
                    plt.stackplot(tvecs[result],y,labels=labels)
                    plt.legend(labelspacing=-matplotlib.rcParams["legend.labelspacing"]-2)
                else:
                    for pop in final_pop_labels:
                        plt.plot(tvecs[result],final_outputs[result][pop][output],label=name(pop,proj))
                    plt.legend()

    elif axis == 'outputs':
        for result in final_result_labels:
            for pop in final_pop_labels:
                figs.append(plt.figure())
                plt.xlabel('Year')
                plt.autoscale(enable=True, axis='x', tight=True)
                plt.ylabel('Mixed')
                plt.title('%s-%s' % (result,pop))
                if plot_type == 'stacked':
                    y = [final_outputs[result][pop][output] for output in final_output_labels]
                    labels = [name(output,proj) for output in final_output_labels]
                    plt.stackplot(tvecs[result],y,labels=labels)
                    plt.legend(labelspacing=-matplotlib.rcParams["legend.labelspacing"]-2)
                else:
                    for output in final_output_labels:
                        plt.plot(tvecs[result],final_outputs[result][pop][output],label=name(output,proj))
                    plt.legend()

    return figs


def getFullName(output_id, proj):
    """
    For a given output_id, returns the user-friendly version of the name. 
    """
    if output_id in proj.settings.charac_specs: # characteristic
        output_id = proj.settings.charac_specs[output_id]['name']
    elif output_id in proj.settings.linkpar_specs: # parameter
        output_id = proj.settings.linkpar_specs[output_id]['name']
    elif output_id in proj.settings.node_specs: # compartment
        output_id = proj.settings.node_specs[output_id]['name']
    elif output_id in proj.data['pops']['label_names'].keys(): # population label
        output_id = proj.data['pops']['label_names'][output_id]
    return output_id





