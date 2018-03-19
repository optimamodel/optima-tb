# %% Imports
import logging
from matplotlib.pyplot import plot

logger = logging.getLogger(__name__)

from optima_tb.utils import odict, OptimaException, nestedLoop
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
import textwrap
from optima_tb.plotting import gridColorMap
from matplotlib.ticker import FuncFormatter
import os
from itertools import product
from matplotlib.patches import Rectangle, Patch
from matplotlib.collections import PatchCollection

def save_figs(figs,path = '.',prefix = '',fnames=None):
    #
    if not isinstance(figs,list):
        figs = [figs]
    if fnames is not None:
        if not isinstance(fnames,list):
            fnames = [fnames]
        assert len(fnames) == len(figs), 'Number of figures must match number of specified filenames'

    for i,fig in enumerate(figs):
        if fnames is not None:
            fname = prefix+fnames[i] + '.png'
        else:
            if fig.get_label() == '':
                continue
            fname = prefix+fig.get_label() + '.png'
        fig.savefig(os.path.join(path,fname))
        logger.info('Saved figure "%s"' % fname)


class PlotData(object):
    # This is what gets passed into a plotting function, which displays a View of the data
    # Conceptually, we are applying visuals to the data.
    # But we are performing an extraction step rather than doing it directly because things like
    # labels, colours, groupings etc. only apply to plots, not to results, and there could be several
    # different views of the same data.

    def __init__(self,results,outputs=None,pops=None,output_aggregation='sum',pop_aggregation='sum',project=None):
        # TODO - Add temporal aggregation
        # Construct a PlotData instance from a Results object by selecting data and optionally
        # specifying desired aggregations
        #
        # final_outputs[result][pop][output] = vals
        # where the keys are the aggregated labels (if specified)
        # Aggregation computation is performed here
        # Input must be a list either containing a list of raw output labels or a dict with a single
        # key where the value is a list of raw output labels e.g.
        # outputs = ['vac',{'total':['vac','sus']}]

        #
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

        # Validate inputs
        if isinstance(results,odict):
            results = [result for _,result in results.items()]
        elif isinstance(results,ResultSet):
            results = [results]

        if pops is None:
            pops = [pop.label for pop in results[0].model.pops]
        elif pops == 'All':
            pops = [{'All':[pop.label for pop in results[0].model.pops]}]
        elif isinstance(pops,str):
            pops = [pops]

        if outputs is None:
            outputs = [comp.label for comp in results[0].model.pops[0].comps if not (comp.tag_birth or comp.tag_dead or comp.is_junction)]
        elif isinstance(outputs,str):
            outputs = [outputs]

        assert isinstance(results,list), 'Results should be specified as a Result, list, or odict'
        assert isinstance(pops, list), 'Populations need to be specified as a string or list'
        assert isinstance(outputs,list), 'Outputs need to be specified as a string or list'

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

        self.series = []
        tvecs = dict()

        for result in results: # For each result

            result_label = result.name
            tvecs[result_label] = result.model.sim_settings['tvec']
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
                            data_dict[output_label] = np.zeros(tvecs[result_label].shape)
                            compsize[output_label] = np.zeros(tvecs[result_label].shape)
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
                            vals = sum(aggregated_outputs[x][output_name] for x in pop_labels) # Add together all the outputs
                        elif pop_aggregation == 'average': 
                            vals = sum(aggregated_outputs[x][output_name] for x in pop_labels) # Add together all the outputs
                            vals /= len(labels)
                        elif pop_aggregation == 'weighted':
                            vals = sum(aggregated_outputs[x][output_name]*popsize[x] for x in pop_labels) # Add together all the outputs
                            vals /= sum([popsize[x] for x in pop_labels])
                        self.series.append(Plottable(tvecs[result_label],vals,result_label,pop_name,output_name))
                    else:
                        vals = aggregated_outputs[pop][output_name]
                        self.series.append(Plottable(tvecs[result_label],vals,result_label,pop,output_name))

        self.results = [x.name for x in results] # NB. These are lists that thus specify the order in which plotting takes place
        self.pops = [x.keys()[0] if isinstance(x,dict) else x for x in pops]
        self.outputs = [x.keys()[0] if isinstance(x,dict) else x for x in outputs]

        if project is not None:
            for s in self.series:
                s.result = getFullName(s.result,project)
                s.pop = getFullName(s.pop,project)
                s.output = getFullName(s.output,project)
            self.results = [getFullName(r,project) for r in self.results]
            self.pops = [getFullName(r,project) for r in self.pops]
            self.outputs = [getFullName(r,project) for r in self.outputs]

    def __repr__(self):
        s = 'PlotData\n'
        s += 'Results: %s\n' % (self.results)
        s += 'Pops: %s\n' % (self.pops)
        s += 'Outputs: %s\n' % (self.outputs)
        return s

    def tvals(self):
        assert len(set([len(x.tvec) for x in self.series])) == 1, 'All series must have the same number of time points' # All series must have the same number of timepoints
        tvec = self.series[0].tvec
        for i in xrange(1,len(self.series)):
            assert(all(self.series[i].tvec == tvec)), 'All series must have the same time points'
        return tvec

    def time_aggregate(self,t_bins = None, method='sum'):
        # t_bins is a vector of bin edges
        # For time points where the time is >= the lower bin value and < upper bin value
        # values will be aggregated according to the method
        # If the method is 'sum' then the output time points will be the upper bin edges
        # If the method is 'average' then the output time points will be the bin centres
        assert method in ['sum','average']

        if t_bins is None:
            t_out = np.zeros((1,))
            lower = [-np.inf]
            upper = [np.inf]
            assert method=='sum', 'Cannot average data over all time'
        else:
            lower = t_bins[0:-1]
            upper = t_bins[1:]
            if method == 'sum':
                t_out = upper
            elif method == 'average':
                t_out = (lower+upper)/2.0

        for s in self.series:
            tvec = []
            vals = []
            for i,low,high,t in zip(range(len(lower)),lower,upper,t_out):
                flt = (s.tvec >= low) & (s.tvec < high)
                tvec.append(t)
                if method == 'sum':
                    vals.append(np.sum(s.vals[flt]))
                elif method == 'average':
                    vals.append(np.average(s.vals[flt]))
            s.tvec = np.array(tvec)
            s.vals = np.array(vals)


    def __getitem__(self,key):
        # key is a tuple of (result,pop,output)
        # retrive a single Series e.g. plotdata['default','0-4','sus']
        for s in self.series:
            if s.result == key[0] and s.pop == key[1] and s.output == key[2]:
                return s
        raise OptimaException('Series %s-%s-%s not found' % (key[0],key[1],key[2]))

    def set_colors(self,colors,results=[],pops=[],outputs=[],overwrite=False):
        # Set the colours for plottables matching the filters applied by results, pops, or outputs
        # Colours will be applied across all other dimensions
        # Set the colours in the Plottables manually for further granularity
        # If overwrite=False, then colours will only be set if the existing color is None
        # That way, None colours can be replaced by default colours while still preserving any
        # user-set colors

        assert len(colors) in [len(results),len(pops),len(outputs)]
            
        for i,color in enumerate(colors):
            
            if results:
                series = [x for x in self.series if x.result == results[i]] 
            elif pops:
                series = [x for x in self.series if x.pop == pops[i]] 
            else:
                series = [x for x in self.series if x.output == outputs[i]] 
            
            for s in series:
                s.color = color if (s.color is None or overwrite==True) else s.color
        
class Plottable(object):
    def __init__(self,tvec,vals,result='default',pop='default',output='default',color=None,label=None):
        self.tvec = np.copy(tvec)
        self.vals = np.copy(vals)
        self.result = result
        self.pop = pop
        self.output = output
        self.color = color # Automatically decide color at runtime (e.g. in plotSeries this could vary depending on whether axis='outputs' or 'pops etc)
        self.label = label # Automatically decide label at runtime (e.g. in plotSeries this could be the output name or the pop name depending on the axis)

def plotBars(plotdata,stack_pops=None,stack_outputs=None,adjacent='times',separate_legend=False):
    # We have a collection of bars - one for each Result, Pop, Output, and Timepoint.
    # Any aggregations have already been done. But _groupings_ have not. Let's say that we can group
    # pops and outputs but we never want to stack results. At least for now. 
    # But, the quantities could still vary over time. So we will have
    # - A set of arrays where the number of arrays is the number of things being stacked and
    #   the number of values is the number of bars - could be time, or could be results?
    # - As many sets as there are ungrouped bars

    plotdata = dcp(plotdata)

    # Note - all of the tvecs must be the same
    tvals = plotdata.tvals() # We have to iterate over these, with offsets, if there is more than one

    # If split_time = True, then different timepoints will be on different figures
    # If split_results = True, then different results will be on different figures

    # First - we should come up with the time/result grouped bars
    # The bar specification is - for each synchronized element of stack_pops and stack_outputs, it will
    # group those quantities into a bar. There could be multiple bars. 
    # If stack pops is not specified, it will default to all populations
    #
    # The rule is - if you want quantities to appear as a single color, they should be 
    # aggregated in plotdata. If you want quantities to appear as separate colors, 
    # they should be stacked in plotBars

    def get_unique(x):
        o = set()
        for y in x:
            if isinstance(y,list):
                o.update(y)
            else:
                o.add(y)
        return o

    if stack_pops is None:
        stack_pops = plotdata.pops
    else:
        stack_pops += list(set(plotdata.pops)-get_unique(stack_pops))

    if stack_outputs is None:
        stack_outputs = plotdata.outputs
    else:
        stack_outputs += list(set(plotdata.outputs)-get_unique(stack_outputs))

    bar_pops = []
    bar_outputs = []

    # TODO - The order of this loop will affect the order the bars are plotted, could
    # have a switch to manage this order
    for pop in stack_pops:
        for output in stack_outputs:
            bar_pops.append(pop if isinstance(pop,list) else [pop])
            bar_outputs.append(output if isinstance(output,list) else [output])


    # Set colours on a per-output basis by default
    plotdata.set_colors(gridColorMap(len(plotdata.outputs)),outputs=plotdata.outputs)

    # now bar_pops and bar_outputs are lists of the same length, with length
    # equal to the total number of bars for a single timepoint-result unit
    # So a set of these is rendered for every result-output combination
    # Then, we need to calculate the offsets if we are putting different
    # results and timepoints on different plots or not
    # Or the alternative is to say that we want separate plots...?
    # i.e. force the user to specify individual ones? No, that would be annoying?
    width = 1.0
    gaps = (0.1,0.5,2.0) # Spacing within blocks, across blocks, and across major blocks

    block_width = len(bar_pops)*(width+gaps[0])

    if adjacent == 'results':
        result_offset = block_width+gaps[1]
        tval_offset = len(plotdata.results)*(block_width+gaps[1])+gaps[2]
        iterator = nestedLoop([range(len(plotdata.results)),range(len(tvals))],[0,1])
    else:
        result_offset = len(tvals)*(block_width+gaps[1])+gaps[2]
        tval_offset = block_width+gaps[1]
        iterator = nestedLoop([range(len(plotdata.results)),range(len(tvals))],[1,0])

    figs = []
    fig,ax = plt.subplots()
    figs.append(fig)

    rectangles = defaultdict(lambda: defaultdict(list)) # Accumulate the list of rectangles for each colour (pop-output combination)
    # NOTE
    # pops, output - colour separates them. To merge colours, aggregate the data first
    # results, time - spacing separates them. Can choose to group by one or the other

    # Now, there are three levels of ticks
    # There is the within-block level, the inner group, and the outer group
    # 
    xticks = []

    for r_idx,t_idx in iterator:
        base_offset = r_idx*result_offset + t_idx*tval_offset
        block_offset = 0.0
        

        for bar_pop,bar_output in zip(bar_pops,bar_outputs): # For each bar within the bar collection that we are going to be plotting
            # pop is something like ['0-4','5-14'] or ['0-4']
            # output is something like ['sus','vac'] or ['0-4'] depending on the stack
            y0 = 0
            for pop in bar_pop:
                for output in bar_output:
                    y = plotdata[plotdata.results[r_idx],pop,output].vals[t_idx]
                    rectangles[pop][output].append(Rectangle((base_offset+block_offset,y0), width, y))
                    y0 += y
            block_offset += width + gaps[0]

    legend_patches = []
    for pop,x in rectangles.items():
        for output,y in x.items():
            if len(get_unique(stack_pops)) == 1:
                label = output
            elif len(get_unique(stack_outputs)) == 1:
                label = pop
            else:
                label = '%s-%s' % (pop,output)
            pc = PatchCollection(y, facecolor=plotdata[plotdata.results[0],pop,output].color)
            ax.add_collection(pc)
            legend_patches.append(Patch(facecolor=plotdata[plotdata.results[0],pop,output].color,label=label))


    if separate_legend:
        figs.append(render_separate_legend(ax),handles=legend_patches)
    else:
        render_legend(ax,handles=legend_patches)

    ax.autoscale()
    ax.set_ylim(ymin=0)
    return fig


def plotSeries(plotdata,plot_type='line',axis='outputs',separate_legend=False,plot_observed_data=True,project=None):
    # TODO -
    # - Clean up doing aggregation separately
    # - Implement separate figures as everything starting out on the same plot and then
    #   distributing them across figures if requested? But does this still have those problems associated
    #   with colours? Maybe the colours should be specified in a dict the same as the data?
    # This function plots a time series for a model output quantities
    #
    # INPUTS
    # - plotdata - a PlotData instance to plot
    # - plot_type - 'line', 'stacked', or 'proportion' (stacked, normalized to 1)
    # - project - Draw scatter points for data wherever the output label matches
    #   a data label. Only draws data if the plot_type is 'line'
    # - separate_legend - Show the legend in a separate figure,

    assert axis in ['outputs','results','pops']

    figs = []

    plotdata = dcp(plotdata)

    # Update colours with defaults, if they were not set

    if axis == 'results':
        plotdata.set_colors(gridColorMap(len(plotdata.pops)),pops=plotdata.pops)
    elif axis == 'pops':
        plotdata.set_colors(gridColorMap(len(plotdata.pops)),pops=plotdata.pops)
    elif axis == 'outputs':
        plotdata.set_colors(gridColorMap(len(plotdata.outputs)),outputs=plotdata.outputs)

    if axis == 'results':
        for pop in plotdata.pops:
            for output in plotdata.outputs:
                fig,ax = plt.subplots()
                fig.set_label('%s_%s' % (pop,output))
                figs.append(fig)
                ax.set_ylabel(name(output,proj))
                ax.set_title('%s' % (pop))
                if plot_type in ['stacked','proportion']:
                    y = np.stack([plotdata[result,pop,output].vals for result in plotdata.results])
                    y = y/np.sum(y,axis=0) if plot_type == 'proportion' else y
                    ax.stackplot(plotdata[plotdata.results[0],pop,output].tvec,y,labels=plotdata.results,colors=[plotdata[result,pop,output].color for result in plotdata.results])
                else:
                    for result in plotdata.results:
                        ax.plot(plotdata[result,pop,output].tvec,plotdata[result,pop,output].vals,color=plotdata[result,pop,output].color,label=result)
                        if project is not None:
                            render_data(ax,project, pop, output,plotdata[result,pop,output].color)
                apply_series_formatting(ax,plot_type)
                if not separate_legend:
                    render_legend(ax,plot_type)

    elif axis == 'pops':
        for result in plotdata.results:
            for output in plotdata.outputs:
                fig,ax = plt.subplots()
                fig.set_label('%s_%s' % (result,output))
                figs.append(fig)
                ax.set_title('%s' % (result))
                if plot_type in ['stacked','proportion']:
                    y = np.stack([plotdata[result,pop,output].vals for pop in plotdata.pops])
                    y = y/np.sum(y,axis=0) if plot_type == 'proportion' else y
                    ax.stackplot(plotdata[result,plotdata.pops[0],output].tvec,y,labels=plotdata.pops,colors=[plotdata[result,pop,output].color for pop in plotdata.pops])
                else:
                    for pop in plotdata.pops:
                        ax.plot(plotdata[result,pop,output].tvec,plotdata[result,pop,output].vals,color=plotdata[result,pop,output].color,label=pop)
                        if project is not None:
                            render_data(ax,project, pop, output,plotdata[result,pop,output].color)
                apply_series_formatting(ax,plot_type)
                if not separate_legend:
                    render_legend(ax,plot_type,)

    elif axis == 'outputs':
        for result in plotdata.results:
            for pop in plotdata.pops:
                fig,ax = plt.subplots()
                fig.set_label('%s_%s' % (result,pop))
                figs.append(fig)
                # plt.ylabel('Mixed')
                ax.set_title('%s-%s' % (result,pop))
                if plot_type in ['stacked','proportion']:
                    y = np.stack([plotdata[result,pop,output].vals for output in plotdata.outputs])
                    y = y/np.sum(y,axis=0) if plot_type == 'proportion' else y
                    ax.stackplot(plotdata[result,pop,plotdata.outputs[0]].tvec,y,labels=plotdata.outputs,colors=[plotdata[result,pop,output].color for output in plotdata.outputs])
                else:
                    for output in plotdata.outputs:
                        ax.plot(plotdata[result,pop,output].tvec,plotdata[result,pop,output].vals,color=plotdata[result,pop,output].color,label=output)
                        if project is not None:
                            render_data(ax,project, pop, output,plotdata[result,pop,output].color)
                apply_series_formatting(ax,plot_type)
                if not separate_legend:
                    render_legend(ax,plot_type)

    if separate_legend:
        figs.append(render_separate_legend(ax,plot_type))

    return figs

def render_data(ax,proj,pop,output,color):
    # This function renders a scatter plot for a single variable (in a single population)
    # The scatter plot is drawn in the current axis
    # INPUTS
    # proj - Project object
    # pop - name of a population (str)
    # output - name of an output (str)
    # name - The name-formatting function to retrieve full names (currently unused)
    # color - The color of the data points to use
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

    ax.scatter(t,y,marker='o', s=40, linewidths=3, facecolors='none',color=color)#label='Data %s %s' % (name(pop,proj),name(output,proj)))

def KMSuffixFormatter(x, pos):
    # This function specifies the formatting for the Y-axis labels

    'The two args are the value and tick position'
    if x >= 1e6:
        return '%1.1fM' % (x * 1e-6)
    # elif x >= 1e3:
    #     return '%1.1fK' % (x * 1e-3)
    else:
        return '%g' % x

def apply_series_formatting(ax,plot_type):
    # This function applies formatting that is common to all Series plots
    # (irrespective of the 'axis' setting)
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.set_xlabel('Year')
    ax.set_ylim(ymin=0)
    if plot_type == 'proportion':
        ax.set_ylim(ymax=1)
        ax.set_ylabel('Proportion ' + ax.get_ylabel())
    else:
        ax.set_ylim(ymax=ax.get_ylim()[1] * 1.05)
    ax.yaxis.set_major_formatter(FuncFormatter(KMSuffixFormatter))

def render_separate_legend(ax,handles=None,plot_type=None):
    if handles is None:
        handles, labels_leg = ax.get_legend_handles_labels()
    else:
        labels = [h.get_label() for h in handles]

    fig = plt.figure()
    legendsettings = {'loc': 'center', 'bbox_to_anchor': None,'frameon':False}

    # A legend renders the line/patch based on the object handle. However, an object
    # can only appear in one figure. Thus, if the legend is in a different figure, the
    # object cannot be shown in both the original figure and in the legend. Thus we need
    # to copy the handles, and use the copies to render the legend
    from copy import copy
    handles = [copy(x) for x in handles]

    # Stop the figures from being rendered in the original figure, which will allow them to
    # then be rendered in the legend figure
    for h in handles:
        h.figure = None

    if plot_type in ['stacked', 'proportion']:
        fig.legend(handles=handles[::-1], labels=labels[::-1], **legendsettings)
    else:
        fig.legend(handles=handles, labels=labels, **legendsettings)

    return fig

def render_legend(ax,handles=None,plot_type=None):
    # This function renders a legend
    # INPUTS
    # - plot_type - Used to decide whether to reverse the legend order for stackplots
    if handles is None:
        handles, labels_leg = ax.get_legend_handles_labels()
    else:
        labels = [h.get_label() for h in handles]

    legendsettings = {'loc': 'center left', 'bbox_to_anchor': (1.05, 0.5), 'ncol': 1}
    labels = [textwrap.fill(label, 16) for label in labels]

    if plot_type in ['stacked', 'proportion']:
        ax.legend(handles=handles[::-1], labels=labels[::-1], **legendsettings)
    else:
        ax.legend(handles=handles, labels=labels,**legendsettings)
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
