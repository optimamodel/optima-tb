# %% Imports
import logging
from matplotlib.pyplot import plot

logger = logging.getLogger(__name__)

from optima_tb.utils import odict, OptimaException
from optima_tb.results import ResultSet
import numpy as np
from optima_tb.project import Project
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
import textwrap
from optima_tb.plotting import gridColorMap
from matplotlib.ticker import FuncFormatter

def compute_aggregations(results,outputs,pops,output_aggregation,pop_aggregation):
    # Return a nested dict
    # final_outputs[result][pop][output] = vals
    # where the keys are the aggregated labels (if specified)
    # Aggregation computation is performed here
    # Input must be a list either containing a list of raw output labels or a dict with a single
    # key where the value is a list of raw output labels e.g.
    # outputs = ['vac',{'total':['vac','sus']}]

    assert isinstance(pops, list), 'Populations need to be specified as a list'
    assert isinstance(outputs,list), 'Outputs need to be specified as a list or a string'

    assert output_aggregation in ['sum','average','weighted']
    assert pop_aggregation in ['sum','average','weighted']

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
    final_tvecs = dict()

    for result_label,result in results.items(): # For each result

        final_outputs[result_label] = defaultdict(dict) # This final_outputs[result_label][agg_pop_label][agg_output_label] - i.e. the final data for plotting
        final_tvecs[result_label] = result.model.sim_settings['tvec']
        dt = result.model.sim_settings['tvec_dt']

        aggregated_outputs = defaultdict(dict) # Dict with aggregated_outputs[pop_label][aggregated_output_label]
        output_units = dict()
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
                    output_units[output_label] = pop.getComp(output_label).units
                elif output_label in pop.charac_ids:
                    data_dict[output_label] = pop.getCharac(output_label).vals
                    compsize[output_label] = data_dict[output_label]
                    output_units[output_label] = pop.getCharac(output_label).units
                elif output_label in pop.par_ids:
                    par = pop.getPar(output_label)
                    if par.links: # If this is a transition parameter, use getFlow to get the flow rate
                        data_dict[output_label] = np.zeros(final_tvecs[result_label].shape)
                        compsize[output_label] = np.zeros(final_tvecs[result_label].shape)
                        for link in par.links:
                            data_dict[output_label] += link.vals/dt
                            compsize[output_label] += (link.source.vals if not link.source.is_junction else link.source.vals_old)
                        output_units[output_label] = link.units
                    else:
                        data_dict[output_label] = pop.getPar(output_label).vals
                        output_units[output_label] = pop.getPar(output_label).units
                else:
                    raise OptimaException('Output "%s" not found in pop "%s"' % (output_label,pop_label))

            # Second pass, aggregate them according to any aggregations present
            for output in outputs: # For each final output
                if isinstance(output,dict):
                    output_name = output.keys()[0]
                    labels = output[output_name]
                    if len(set([output_units[x] for x in labels])) > 1:
                        logger.warn('Warning - aggregation for output "%s" is mixing units, this is almost certainly not desired' % (output_name))
                    if output_aggregation == 'sum': 
                        aggregated_outputs[pop_label][output_name] = sum(data_dict[x] for x in labels) # Add together all the outputs
                    elif output_aggregation == 'average': 
                        aggregated_outputs[pop_label][output_name] = sum(data_dict[x] for x in labels) # Add together all the outputs
                        aggregated_outputs[pop_label][output_name] /= len(labels)
                    elif output_aggregation == 'weighted':
                        aggregated_outputs[pop_label][output_name] = sum(data_dict[x]*compsize[x] for x in labels) # Add together all the outputs
                        aggregated_outputs[pop_label][output_name] /= sum([compsize[x] for x in labels])
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
                    final_outputs[result_label][pop][output_name] = aggregated_outputs[pop][output_name]

    # Now, we have completed all aggregations. The final_outputs dict is now ready to be plotted across the required axis
    final_labels = dict()
    final_labels['results'] = results.keys()
    final_labels['pops'] = [x.keys()[0] if isinstance(x,dict) else x for x in pops]
    final_labels['outputs'] = [x.keys()[0] if isinstance(x,dict) else x for x in outputs]

    return final_outputs, final_tvecs, final_labels

def plotSeries(proj,results,outputs=None,pops=None,axis='outputs',output_aggregation='sum',pop_aggregation='sum',plot_type='line',use_full_labels=True,separate_legend=False,plot_observed_data=True,colors=None):
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
    # - colors should be the same length as whichever quantity is being plotted
    assert isinstance(proj,Project)
    if isinstance(results,ResultSet):
        results = {results.name:results}
    if pops is None:
        pops = [pop.label for pop in results[results.keys()[0]].model.pops]
    elif pops == 'All':
        pops = [{'All':[pop.label for pop in results[results.keys()[0]].model.pops]}]
    if outputs is None:
        outputs = [comp.label for comp in results[results.keys()[0]].model.pops[0].comps if not (comp.tag_birth or comp.tag_dead or comp.is_junction)]
    if isinstance(pops,str):
        pops = [pops]
    if isinstance(outputs,str):
        outputs = [outputs]

    assert axis in ['outputs','results','pops']

    if axis == 'outputs':
        if colors is not None:
            assert len(colors) == len(outputs), 'Number of colors does not match number of outputs'
        else:
            colors = gridColorMap(len(outputs))
    elif axis == 'results':
        if colors is not None:
            assert len(colors) == len(results), 'Number of colors does not match number of results'
        else:
            colors = gridColorMap(len(results))
    elif axis == 'pops':
        if colors is not None:
            assert len(colors) == len(pops), 'Number of colors does not match number of pops'
        else:
            colors = gridColorMap(len(pops))

    final_outputs,final_tvecs,final_labels = compute_aggregations(results,outputs,pops,output_aggregation,pop_aggregation)
    figs = []

    if use_full_labels:
        name = getFullName
    else:
        name = lambda x,y: x
    
    if axis == 'results':
        for pop in final_labels['pops']:
            for output in final_labels['outputs']:
                figs.append(plt.figure())
                plt.ylabel(name(output,proj))
                plt.title('%s' % (pop))
                if plot_type in ['stacked','proportion']:
                    y = np.stack([final_outputs[result][pop][output] for result in final_labels['results']])
                    y = y/np.sum(y,axis=0) if plot_type == 'proportion' else y
                    labels = [name(result,proj) for result in final_labels['results']]
                    plt.stackplot(final_tvecs[result],y,labels=labels,colors=colors)
                else:
                    for result,color in zip(final_labels['results'],colors):
                        plt.plot(final_tvecs[result],final_outputs[result][pop][output],color=color,label=name(result,proj))
                        if plot_observed_data:
                            plot_data(proj, pop, output,name,color)
                apply_formatting(plot_type)
                render_legend(plot_type,separate_legend)

    elif axis == 'pops':
        for result in final_labels['results']:
            for output in final_labels['outputs']:
                figs.append(plt.figure())
                plt.ylabel(name(output,proj))
                plt.title('%s' % (result))
                if plot_type in ['stacked','proportion']:
                    y = np.stack([final_outputs[result][pop][output] for pop in final_labels['pops']])
                    y = y/np.sum(y,axis=0) if plot_type == 'proportion' else y
                    labels = [name(pop,proj) for pop in final_labels['pops']]
                    plt.stackplot(final_tvecs[result],y,labels=labels,colors=colors)
                else:
                    for pop,color in zip(final_labels['pops'],colors):
                        plt.plot(final_tvecs[result],final_outputs[result][pop][output],color=color,label=name(pop,proj))
                        if plot_observed_data:
                            plot_data(proj, pop, output,name,color)
                apply_formatting(plot_type)
                render_legend(plot_type, separate_legend)

    elif axis == 'outputs':
        for result in final_labels['results']:
            for pop in final_labels['pops']:
                figs.append(plt.figure())
                # plt.ylabel('Mixed')
                plt.title('%s-%s' % (result,name(pop,proj)))
                if plot_type in ['stacked','proportion']:
                    y = np.stack([final_outputs[result][pop][output] for output in final_labels['outputs']])
                    y = y/np.sum(y,axis=0) if plot_type == 'proportion' else y
                    labels = [name(output,proj) for output in final_labels['outputs']]
                    plt.stackplot(final_tvecs[result],y,labels=labels,colors=colors)
                else:
                    for output,color in zip(final_labels['outputs'],colors):
                        plt.plot(final_tvecs[result],final_outputs[result][pop][output],color=color,label=name(output,proj))
                        if plot_observed_data:
                            plot_data(proj, pop, output,name,color)
                apply_formatting(plot_type)
                render_legend(plot_type,separate_legend)
    return figs

def plot_data(proj,pop,output,name,color):
    data = proj.data

    if output in data['characs'].keys():
        d = data['characs'][output]
    elif output in data['linkpars'].keys():
        d = data['linkpars'][output]
    else:
        return

    if pop in d:
        y = d[pop]['y']
        t = d[pop]['t']
    else:
        return

    plt.scatter(t,y,marker='o', s=40, linewidths=3, facecolors='none',color=color)#label='Data %s %s' % (name(pop,proj),name(output,proj)))

def KMSuffixFormatter(x, pos):
    'The two args are the value and tick position'
    if x >= 1e6:
        return '%1.1fM' % (x * 1e-6)
    # elif x >= 1e3:
    #     return '%1.1fK' % (x * 1e-3)
    else:
        return '%g' % x

def apply_formatting(plot_type):
    ax = plt.gca()
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.xlabel('Year')
    ax.set_ylim(ymin=0)
    if plot_type == 'proportion':
        ax.set_ylim(ymax=1)
        ax.set_ylabel('Proportion ' + ax.get_ylabel())
    else:
        ax.set_ylim(ymax=ax.get_ylim()[1] * 1.05)
    ax.yaxis.set_major_formatter(FuncFormatter(KMSuffixFormatter))


def render_legend(plot_type,separate_legend):
    ax = plt.gca()
    handles, labels_leg = ax.get_legend_handles_labels()

    if separate_legend:
        fig = plt.figure()
        legendsettings = {'loc': 'center', 'bbox_to_anchor': None,'frameon':False}
        from copy import copy
        handles = [copy(x) for x in handles]
        for h in handles:
            h.figure = None
        if plot_type == 'stacked':
            fig.legend(handles=handles[::-1], labels=labels_leg[::-1], **legendsettings)
        else:
            fig.legend(handles=handles, labels=labels_leg, **legendsettings)
    else:
        legendsettings = {'loc': 'center left', 'bbox_to_anchor': (1.05, 0.5), 'ncol': 1}
        labels_leg = [textwrap.fill(label, 16) for label in labels_leg]
        if plot_type == 'stacked':
            ax.legend(handles=handles[::-1], labels=labels_leg[::-1], **legendsettings)
        else:
            ax.legend(handles=handles, labels=labels_leg,**legendsettings)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])


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
