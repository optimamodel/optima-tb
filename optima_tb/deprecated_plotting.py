import logging
from matplotlib.pyplot import plot

logger = logging.getLogger(__name__)

from optima_tb.utils import odict, OptimaException

import numpy as np

import pylab as pl
from copy import deepcopy as dcp
from copy import copy as ndcp
import matplotlib.cm as cmx
import matplotlib.colors as colors
from random import shuffle
import numbers

logging.warn("Deprecated plotting functions are scheduled for removal")

from optima_tb.plotting import * # Import all public functions
from optima_tb.plotting import _getColormapRange,_turnOffBorder,_calcRelativeDatapoint,_extractDatapoint,_convertPercentage,_plotYIntercept,_plotTrends,_plotBars # Import current set of private functions

def plotPopulation(results, data, pop_labels=None, title='', colormappings=None,
                   plot_observed_data=True, plot_observed_label="alive", save_fig=False, fig_name=None,
                   use_full_labels=True, comp_labels=None, plotdict=None, cat_labels=None):
    """

    Plot all compartments for all populations

    Params:
        results         Results object
        data            project data
        pop_labels      list of populations to be plotted.
                        Default: None, all populations will be plotted
        comp_labels     list of compartments to be plotted. Note that only plottable compartments (as defined in isPlottableComp())
                        are plotted.
                        Default: None, all compartments will be plotted
        title           String, for title of plots
        colormappings    odict with comp.labels as keys and corresponding rgba colors as values
        cat_labels       m
        plot_observed_data    boolean flag: whether observed data should be plotted
        plot_observed_label
        save_fig        boolean flag, whether to save figure as "Full_compartment_<pop_label>"
        use_full_labels     boolean flag: whether to use full compartment label (i.e. 'Susceptible') or short label version (i.e. 'sus')
                            The latter is easier for debugging purposes
        plotdict


    Example usage:

        # Plot all compartments for 15-64, using default color mapping, with no observed data

        pop_labels = ['15-64']
        plotPopulation(results, data, pop_labels, plot_observed_data=False)



        # Plot only active compartments for all populations, using a specific colormap and plotting
        # observed data (within the data object, with the label 'ac_inf')

        cat_list = odict()
        cat_list['#005B9A'] = ['sus']
        cat_list['#0191C8'] = ['vac']
        cat_list['Purples'] = ['lteu', 'ltlu','lted','ltet','ltld','ltlt']
        cat_list['Oranges'] = ['spdu', 'spdd', 'spdt', 'spmu', 'spmd', 'spmt', 'spxu', 'spxd', 'spxt']
        cat_list['Reds']    = ['sndu', 'sndd', 'sndt', 'snmu', 'snmd', 'snmt', 'snxu', 'snxd', 'snxt']
        cat_list['Greens']  = ['acr','ltr']
        labels = ['Susceptible','Vaccinated','Latent TB','Active TB (S+)','Active TB (S-)','Recovered']
        plot_comp_labels = ['spdu', 'sndu', 'spdd', 'sndd', 'spdt', 'sndt',
                        'spmu', 'snmu', 'spmd', 'snmd', 'spmt', 'snmt',
                        'spxu', 'spxd', 'spxt']

        plotPopulation(results=results,
                   data=proj.data,
                   pop_labels=None,
                   title="Active TB",
                   colormappings=cat_list,
                   cat_labels=labels,
                   plot_observed_label="ac_inf",
                   comp_labels=plot_comp_labels)

    """
    
    # setup data structures
    if pop_labels is None:
        pop_labels = results.pop_labels
    tvec = results.sim_settings['tvec']
    year_inc = 5.  # TODO: move this to setting
    if plotdict is not None and 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]

    yr_range = np.arange(start_year, end_year + 0.1, year_inc, dtype=int)
    mpops = results.model.pops
    sim_settings = results.sim_settings
    save_figname = None
    dataobs = None  # placeholder for observed data


    if use_full_labels:
        ncol = 1
    else:
        ncol = 2

    # setup: determine compartment indices to be plotted
    if comp_labels is None:
        comp_labels = sorted(mpops[0].comp_ids, key=mpops[0].comp_ids.get)

    # select only compartments that are plottable
    comp_labels = [comp_label for comp_label in comp_labels if isPlottableComp(comp_label, sim_settings, results.comp_specs)]


    # setup: determine colors to be used
    colors = []
    if colormappings is not None:
        colors_dict, cat_colors = getCategoryColors(colormappings, plotdict['colormapping_order'])
        # reorder so that colors are same as expected for plotting the population
        for (j, comp_label) in enumerate(comp_labels):
            if isPlottableComp(comp_label, sim_settings, results.comp_specs):
                colors.append(colors_dict[comp_label])

    # setup: plotting dict structures
    if plotdict is None:
        plotdict = {}


    # preselect for populations

    y_values, pop_labels, labels, dataobs = extractCompartment(results, data,
                                                     pop_labels=pop_labels,
                                                     comp_labels=comp_labels,
                                                     plot_observed_data=plot_observed_data,
                                                     plot_observed_label=plot_observed_label,
                                                     use_full_labels=True)

    if dataobs is not None:
        that, yhat = dataobs

    # iterate for each key population group
    for (i, poplabel) in enumerate(pop_labels):

        pl_title = title + ' Population: %s' % (poplabel)
        if save_fig:
            save_figname = fig_name + "_compartments_%s" % poplabel

        pdict = {  'ymin': 0,  # required
                  'xlabel': 'Year',
                  'year_inc' :  5.,
                  'ylabel': 'People',
                  'mec' : 'k',
                  'title' : pl_title,
                  'x_ticks' : (yr_range, yr_range),
                  'colors': colors,
                  'save_figname' : save_figname,
                  'save_fig': save_fig,
                  }

        if dataobs is not None:
            pdict['datapoints'] = (that[i], yhat[i])

        pdict.update(plotdict)

        legendsettings = {'loc':'center left',
                           'bbox_to_anchor':(1.05, 0.5),
                           'ncol':ncol}

        _plotStackedCompartments(tvec, y_values[i][:], labels,
#                                  legendsettings=legendsettings,
                                 catlabels=cat_labels, catcolors=colors,
                                 # save_fig=save_fig, save_figname=save_figname,
                                 **pdict)


    if pdict.has_key('legend_off') and pdict['legend_off']:
        # Do this separately to main iteration so that previous figure are not corrupted
        # Note that colorlist may be different to colors, as it represents
        # classes of compartments
        separateLegend(labels=cat_labels, colors=cat_colors, fig_name=fig_name, reverse_order=True, **legendsettings)





def _calculateDatapoints(data, data_labels, pop):
    

    yvals = None
    for (i, label) in enumerate(data_labels):
        ys = data['characs'][label][pop]['y']
        ts = data['characs'][label][pop]['t']

        if i == 0:
            yvals = ys
        else:
            yvals += ys
    return yvals, ts


def extractCompartment(results, data, pop_labels=None, comp_labels=None,
                       plot_observed_data=True, plot_observed_label="alive", use_full_labels=False):
    """
    Wrapper method to extract compartments for a subset of populations and compartments:

    Params:
        results            results
        data                project data
        pop_labels          list of populations to extrac. Default: None, which returns all populations
        comp_labels         list of compartments to extract. Default: None, which returns all compartments
        plot_observed_data    flag indicating whether to extract observed datapoints
        plot_observed_label   label of location of observed datapoints
        use_full_labels
    """
    

    datapoints, pop_labels, comp_labels = results.getCompartmentSizes(pop_labels=pop_labels, comp_label=comp_labels, use_observed_times=False)
    yhat, that = [], []

    if use_full_labels:
        labels = [comp.label for comp in results.model.pops[0].comps]
    else:
        labels = comp_labels

    if plot_observed_data:
        for pop in pop_labels:

            if isinstance(plot_observed_label, basestring):
                ys, ts = _calculateDatapoints(data, [plot_observed_label], pop)
            elif isinstance(plot_observed_label, list):
                ys, ts = _calculateDatapoints(data, plot_observed_label, pop)
            else:
                logger.warn("Unknown data characteristic: ")
                logger.warn(plot_observed_label)

            yhat.append(ys)
            that.append(ts)

        dataobs = (that, yhat)
    else:
        dataobs = None

    return datapoints, pop_labels, comp_labels, dataobs



def extractCharacteristic(results, data, charac_specs, charac_labels=None, pop_labels=None, plot_observed_data=True, plot_total=False):
    """
    Wrapper method to extract characteristics for a subset of populations:

    Params:
        results            results
        data                project data
        charac_specs
        pop_labels          list of populations to extrac. Default: None, which returns all populations
        charac_labels         list of characteristics to extract. Default: None, which returns all characteristics
        plot_observed_data    flag indicating whether to extract observed datapoints
        plot_total            flag indicating whether to sum across populations

    """
    

    datapoints = results.getCharacteristicDatapoints(pop_label=pop_labels, char_label=charac_labels, use_observed_times=False)[0]

    unit_tags = []
    dataobs = None
    yhat, that = [], []
    labels = pop_labels

    # some post processing to make our lives easier (and plots prettier):
    # 1) prettify percentages format i.e. from 0.12 --> 12%

    for k, output_id in enumerate(charac_labels):

        if output_id in charac_specs and 'plot_percentage' in charac_specs[output_id].keys():
            y_values = [datapoints[output_id][i] for i, p in enumerate(pop_labels)]
            y_values = np.array(y_values)
            y_values *= 100
            datapoints[output_id] = y_values
            unit_tags.append(' (%)')
        else:
            unit_tags.append('')


        # 2) plot observed data. Note that we don't want to returned observed data for totals
        if plot_observed_data and not plot_total and output_id in data['characs'].keys():

            ys = [data['characs'][output_id][poplabel]['y'] for poplabel in pop_labels]
            ts = [data['characs'][output_id][poplabel]['t'] for poplabel in pop_labels]

            if 'plot_percentage' in charac_specs[output_id].keys():
                ys *= 100

        else:  # For the case when plottable characteristics were not in the databook and thus not converted to data.
            ys = []
            ts = []

        yhat.append(ys)
        that.append(ts)

        # 3) plot as total for a characteristic across populations. Note that we shouldn't plot
        #    totals for percentages, but we won't enforce this (for the moment)
        if plot_total:

            y_values = [datapoints[output_id][i] for i, p in enumerate(pop_labels)]
            y_values = np.array(y_values)
            y_values = [y_values.sum(axis=0)]
            datapoints[output_id] = y_values

            labels = ['Total']

    dataobs = (that, yhat)


    return datapoints, labels, unit_tags, dataobs


def extractFlows(pop_labels, comp_label, results, settings, tvec, link_labels=None, include_link_not_exclude=True, link_legend=None,
                  plot_inflows=True, plot_outflows=True, exclude_transfers=False, sum_total=False, sum_population=False):
    """
    Wrapper method to extract flows for a subset of populations and compartments:

    Params:
        pop_labels         list of populations to extrac.
        comp_label         list of characteristics to extract.
        results            results
        setting            project setting
        tvec
        sum_total            flag indicating whether to sum across the flows
        sum_population       flag indicating whether to sum across the populations

    """
    

    all_rates = []
    all_tvecs = []

    for (j, pid) in enumerate(pop_labels):

        comp = results.model.pops[pid].getComp(comp_label)

        all_labels = []
        pop_rates = []
        pop_tvecs = []
        for in_out in xrange(2):
            if (in_out == 0 and plot_inflows) or (in_out == 1 and plot_outflows):
                comp_link_ids = [comp.inlink_ids, comp.outlink_ids][in_out]
                label_tag = ['In: ', 'Out: '][in_out]
                for link_tuple in comp_link_ids:
                    link = results.model.pops[link_tuple[0]].links[link_tuple[1]]
                    # print link.label, link_labels, include_link_not_exclude
                    if link_labels is None or (include_link_not_exclude and link.label in link_labels) or (not include_link_not_exclude and link.label not in link_labels):
                        try:
                            legend_label = label_tag + settings.linkpar_specs[link.label]['name']
                        except:
                            if exclude_transfers: continue
                            else: legend_label = label_tag + link.label
                        if link.label in link_legend:
                            legend_label = link_legend[link.label]  # Overwrite legend for a label.
                        num_flow = dcp(link.vals)
                        if in_out == 0:
                            comp_source = results.model.pops[link.index_from[0]].comps[link.index_from[1]]
                        else:
                            comp_source = comp
                        was_proportion = False
                        if link.val_format == 'proportion':
                            denom_val = sum(results.model.pops[lid_tuple[0]].links[lid_tuple[-1]].vals for lid_tuple in comp_source.outlink_ids)
                            num_flow /= denom_val
                            was_proportion = True
                        if link.val_format == 'fraction' or was_proportion is True:
                            if was_proportion is True:
                                num_flow *= comp_source.vals_old
                            else:
                                num_flow[num_flow > 1.] = 1.
                                num_flow = 1 - (1 - num_flow) ** results.dt  # Fractions must be converted to effective timestep rates.
                                num_flow *= comp_source.vals
                            num_flow /= results.dt  # All timestep-based effective fractional rates must be annualised.

                        all_labels.append(legend_label)
                        pop_rates.append(num_flow)
                        pop_tvecs.append(tvec)

        if sum_total:
            pop_tvecs = pop_tvecs[:1]
            all_labels = ['Total summed movement']
            pop_rates_tmp = np.array(pop_rates)
            pop_rates = [pop_rates_tmp.sum(axis=0)]

        all_rates.append(pop_rates[0])
        all_tvecs.append(pop_tvecs[0])



    if sum_population:
        all_tvecs = all_tvecs[:1]
        all_labels = ['Total']
        all_rates_tmp = np.array(all_rates)
        all_rates = [all_rates_tmp.sum(axis=0)]

    return all_rates, all_tvecs, all_labels



def plotSingleCompartmentFlow(results, settings, comp_labels=None, comp_titles=None, plot_pops=None, pop_labels=None, pop_titles=None,
              link_labels=None, include_link_not_exclude=True, link_legend=None, sum_total=False,
              plot_inflows=True, plot_outflows=True, exclude_transfers=False, observed_data=None,
              save_fig=False, fig_name=None, colors=None):
    """
    Plot flows rates in and out of a compartment.

    # TODO complete

    """
    

    plotdict = settings.plot_settings
    year_inc = 5.  # TODO remove hardcoded ref
    tvec = results.sim_settings['tvec']
    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]
    yr_range = np.arange(start_year, end_year + 0.1, year_inc, dtype=int)


    if pop_labels is None:
        pop_labels = results.pop_labels


    if link_legend is None: link_legend = dict()

    if plot_pops is not None:
        plot_pids = getPIDs(results, pop_labels)
    else:
        plot_pids = range(len(results.model.pops))
        plot_pops = [pop.label for pop in results.model.pops]

    if comp_labels is None:
        logger.info("No compartments have been selected for flow-plots.")
        comp_labels = []

    if comp_titles is not None and len(comp_titles) != len(comp_labels):
        logger.error("Flow-plot failure due to the number of compartment plot titles not matching the number of compartments to analyse.")
    if pop_titles is not None and len(pop_titles) != len(pop_labels):
        logger.error("Flow-plot failure due to the number of population plot titles not matching the number of populations to analyse.")


    for (i, comp_label) in enumerate(comp_labels):


        for (j, pid) in enumerate(plot_pids):

            plot_label = plot_pops[j]

            comp = results.model.pops[pid].getComp(comp_label)

            all_rates, all_tvecs, all_labels = extractFlows(comp=comp,
                                                            results=results,
                                                            settings=settings,
                                                            tvec=tvec,
                                                            link_labels=link_labels,
                                                            include_link_not_exclude=include_link_not_exclude,
                                                            link_legend=link_legend,
                                                            plot_inflows=plot_inflows,
                                                            plot_outflows=plot_outflows,
                                                            sum_total=sum_total,
                                                            exclude_transfers=exclude_transfers)

            if comp_titles is not None:
                title_comp = comp_titles[i]
            else:
                title_comp = 'Compartment: "%s"' % settings.node_specs[comp_label]['name']
            if pop_titles is not None:
                title_pop = plot_pops[j]
            else:
                title_pop = '\nPopulation: "%s"' % pop_labels[j]


            final_dict = {
              'ylim' : 0,
              'xlabel':'Year',
              'ylabel': 'Number of People',
              'x_ticks' : (yr_range, yr_range),
              'title': title_comp + title_pop,
              'save_figname': '%s_FlowComparision_%s_%s' % (fig_name, comp_label, plot_label)
              }

            if observed_data is not None:
                final_dict['y_hat'] = [observed_data[0]]
                final_dict['t_hat'] = [observed_data[1]]
            final_dict.update(plotdict)


            if len(all_rates) > 0:
                _plotLine(ys=all_rates, ts=all_tvecs, labels=all_labels, colors=colors, save_fig=save_fig, **final_dict)
            else:
                logger.warn("No flows selected for plotting")


def plotPopulationFlows(results, settings, comp_labels=None, comp_titles=None, pop_labels=None,
              link_labels=None, include_link_not_exclude=True, link_legend=None, sum_total=False, sum_population=False,
              plot_inflows=True, plot_outflows=True, exclude_transfers=False, observed_data=None,
              save_fig=False, fig_name=None, colors=None, colormappings=None, linestyles=None, **kwargs):
    """
    Plot flows rates in and out of a compartment, across populations. Intended usage is
    for plots such as total new infections, deaths, etc. where the net flow is required.

    Params:
        results
        settings
        comp_labels        list of compartments
        comp_titles        list of titles for the plots. Should be the same length as comp_labels
        pop_labels         list of population labels
        link_labels        ? list of labels for each flows
        include_link_not_exclude    flag indicating whether to invert choice
        link_legend          list with y-labels
        sum_total            flag indicating whether to sum across all outflows / inflows
        sum_population       flag indicating whether to sum across all populations
        plot_inflows         flag indicating whether to plot inflows
        plot_outflows        flag indicating whether to plot outflows
        exclude_transfers    flag indicating whether to consider transfer from/to other populations
        observed_data        tuple of list of observed values (phew) with ([y_values],[x_values])
        save_fig             flag indicating whether to save
        fig_name             name of filename
        colors               list of colors, either as hex format

    Returns:



    Example:
        pop_labels = ['15-64','65+']
        single_color = ['#562387']
        include_link_not_exclude = True
        plot_inflows = True
        plot_outflows = False
        exclude_transfers = True
        comp_labels = ['spxu']
        comp_titles = ['New XDR infections']
        link_legend = ['Number of new cases']
        link_labels = None
        sum_total = True

        # plot new XDR infections for 15-64, and 65+ populations
        plotPopulationFlows(results=results, settings=proj.settings,
                      comp_labels=[comp_label], comp_titles=[comp_titles[i]+' per population'],
                      pop_labels=pop_labels,
                      link_labels=link_labels, include_link_not_exclude=include_link_not_exclude,
                      link_legend=link_legend, sum_total=sum_total, sum_population = False,
                      plot_inflows=plot_inflows, plot_outflows=plot_outflows,
                      exclude_transfers=exclude_transfers, colors=single_color,
                      save_fig=True, fig_name=fig_name+"_XDRIncidencePop")

       pop_labels = None
       # plot total new XDR infections across populations
        plotPopulationFlows(results=results, settings=proj.settings,
                      comp_labels=[comp_label], comp_titles=[comp_titles[i]+' total'],
                      pop_labels=pop_labels,
                      link_labels=link_labels, include_link_not_exclude=include_link_not_exclude,
                      link_legend=link_legend, sum_total=sum_total, sum_population = True,
                      plot_inflows=plot_inflows, plot_outflows=plot_outflows,
                      exclude_transfers=exclude_transfers, colors=single_color,
                      save_fig=True, fig_name=fig_name+"_XDRIncidencePop")

    """
    

    plotdict = settings.plot_settings
    year_inc = 5.  # TODO remove hardcoded ref
    tvec = results.sim_settings['tvec']
    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]
    yr_range = np.arange(start_year, end_year + 0.1, year_inc, dtype=int)

    if link_legend is None: link_legend = dict()

    if pop_labels is not None:
        plot_pids = getPIDs(results, pop_labels)
    else:
        pop_labels = [pop.label for pop in results.model.pops]
        plot_pids = range(len(results.model.pops))

    if comp_labels is None:
        logger.error("No compartments have been selected for flow-plots.")
        comp_labels = []

    if comp_titles is not None and len(comp_titles) != len(comp_labels):
        logger.error("Flow-plot failure due to the number of compartment plot titles not matching the number of compartments to analyse.")

    if colors is not None and len(colors) >= len(pop_labels):
        pass  # colors as defined in the args should be used as is
    elif colormappings is not None and colors is None:
        colors = []
        colors_dict, cat_colors = getCategoryColors(colormappings, 'sequential')
        # reorder so that colors are same as expected for plotting the population
        for (j, pop_label) in enumerate(pop_labels):
            colors.append(colors_dict[pop_label])
    else:
        colors = gridColorMap(len(pop_labels))
        logger.info("Plotting: setting color scheme to be default colormap, as not all lines had color assigned")

    if linestyles is not None:
        # convert from odict (key: style) to odict (key: population)
        linestyles = getLinemapping(linestyles)
        linestyles = [linestyles[pop] for pop in pop_labels]



    for (i, comp_label) in enumerate(comp_labels):

        rates, tvecs, all_labels = extractFlows(pop_labels=plot_pids,
                                                comp_label=comp_label,
                                                results=results,
                                                settings=settings,
                                                tvec=tvec,
                                                link_labels=link_labels,
                                                include_link_not_exclude=include_link_not_exclude,
                                                link_legend=link_legend,
                                                plot_inflows=plot_inflows,
                                                plot_outflows=plot_outflows,
                                                sum_total=sum_total,
                                                sum_population=sum_population,
                                                exclude_transfers=exclude_transfers)



        if comp_titles is not None:
            title_comp = comp_titles[i]
        else:
            title_comp = 'Compartment: "%s"' % settings.node_specs[comp_label]['name']

        if sum_population:
            labels = all_labels
        else:
            labels = pop_labels

        final_dict = {
          'ylim' : 0,
          'xlabel':'Year',
          'ylabel': 'Number of People',
          'x_ticks' : (yr_range, yr_range),
          'title': title_comp,
          'save_figname': '%s_FlowComparision_%s' % (fig_name, comp_label)
          }

        if observed_data is not None:
            final_dict['y_hat'] = observed_data[0]
            final_dict['t_hat'] = observed_data[1]
        final_dict.update(plotdict)

        if len(rates) > 0:
            _plotLine(ys=rates, ts=tvecs, labels=labels, colors=colors, linestyles=linestyles,
                      save_fig=save_fig, **final_dict)
        else:
            logger.warn("No flows selected for plotting")

    # TODO: plot separate legend


def plotCharacteristic(results, settings, data, title='', outputIDs=None, y_bounds=None,
                       pop_labels=None, plot_total=False,
                       plot_observed_data=True, save_fig=False, fig_name=None, linestyles=None,
                       colormappings=None, colors=None, plotdict=None, legendsettings=None):
    """
    Plot a characteristic across all populations

    Params:
        results
        settings         project settings
        data             observed data
        title            label for plot
        outputIDs        list of characeristics which will be selectively be plotted
                         Default: None, which causes all labels to be plotted
        pop_labels          list of labels for subset of populations to be plotted. If None, all populations are plotted
        plot_total       flag indicating whether to sum and plot the total.
        plot_observed_data plot observed data points on top of simulated data. Useful for calibration
        save_fig         bool to indicate whether to save the figure to file
        fig_name         name to save figures to
        colors           list of colors for populations
        plotdict         project settings plotting dictionary


    Example usage:
        dict = {'y_hat': yhat,
                't_hat': that,
                'unit_tag': '(%)', # default = ''
                'xlabel':'Year',
                'ylabel': charac_specs[output_id]['name'] + '(%)',
                'title': '%s Outputs - %s' % (title, charac_specs[output_id]['name']),
                'marker': 'o',
                'x_ticks' : ([2000,2030],[2000,2030]),
                'save_figname': 'MyPlot'}

        charac_labels = ['ac_inf','dead']
        plotCharacteristic(results=results, settings=settings, outputIDs=charac_labels, pop_labels=pop_labels, data=data, plot_observed_data=plot_observed_data, save_fig=save_fig, fig_name=fig_name, plotdict=plotdict)



    """
    
    charac_specs = settings.charac_specs
    # setup:
    if outputIDs is None:
        outputIDs = results.outputs.keys()
        # now only select characteristics which are plottable
        outputIDs = [output_id for output_id in outputIDs if isPlottableCharac(output_id, charac_specs)]
    # else we already know which characteristics to plot
#     print "-------", outputIDs


    if pop_labels is None:
        pop_labels = [pop.label for pop in results.model.pops]
#     print "-----2--", pop_labels

    if plotdict is None:
        plotdict = settings.plot_settings

    if legendsettings is None:
        legendsettings = {}

    tvec = results.sim_settings['tvec']

    year_inc = 5.  # TODO: move this to setting
    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[-1]
    yr_range = np.arange(start_year, end_year + 0.1, year_inc, dtype=int)

    # placeholder for y_bounds:
    yb = None

    if colors is not None and len(colors) >= len(pop_labels):
        pass  # colors as defined in the args should be used as is
    elif colormappings is not None and colors is None:
        colors = []
        colors_dict, cat_colors = getCategoryColors(colormappings, 'sequential')
        # reorder so that colors are same as expected for plotting the population
        for (j, pop_label) in enumerate(pop_labels):
            colors.append(colors_dict[pop_label])
    else:
        colors = gridColorMap(len(pop_labels))
        logger.info("Plotting: setting color scheme to be default colormap, as not all lines had color assigned")

    if linestyles is not None:
        # convert from odict (key: style) to odict (key: population)
        linestyles = getLinemapping(linestyles)
        linestyles = [linestyles[pop] for pop in pop_labels]


    # extract all characteristics we're interested in, all at once
    y_values, labels, unit_tags, dataobs = extractCharacteristic(results, data, charac_specs, charac_labels=outputIDs, pop_labels=pop_labels, plot_observed_data=plot_observed_data, plot_total=plot_total)

    if dataobs is not None:
        that, yhat = dataobs

    # now plot through each characteristic
    for i, output_id in enumerate(outputIDs):

        if output_id in charac_specs:
            name = charac_specs[output_id]['name']
        elif output_id in settings.linkpar_specs:
            name = settings.linkpar_specs[output_id]['name']
        else:
            raise OptimaException('ERROR: Attempting to plot characteristic "%s" but cannot locate it in either characteristic or parameter specs.' % output_id)

        final_dict = {
                  'unit_tag': unit_tags[i],
                  'xlabel':'Year',
                  'ylabel': name + unit_tags[i],
                  'x_ticks' : (yr_range, yr_range),
                  'title': '%s\n%s' % (title, name),
                  'save_figname': '%s_characteristic_%s' % (fig_name, name)}

        if dataobs is not None:  # this can be improved
            final_dict['y_hat'] = yhat[i]
            final_dict['t_hat'] = that[i]
        print "/////////////"
        print plotdict
        final_dict.update(plotdict)
        print final_dict
        print "\\\\\\\\\\\\\\\\\\\\\|||"

        if y_bounds is not None:
            yb = y_bounds[i]

        _plotLine(y_values[output_id][:], np.tile(tvec, (len(labels), 1)), labels, y_bounds=yb,
#                 legendsettings=None,
                save_fig=save_fig, colors=colors, # ylim=(0, 1900000),
                  linestyles=linestyles, **final_dict)

    if final_dict.has_key('legend_off') and final_dict['legend_off']:
        # Do this separately to main iteration so that previous figure are not corrupted
        # Note that colorlist may be different to colors, as it represents
        # classes of compartments

        separateLegend(labels=labels, colors=colors, fig_name=fig_name + "_LegendCharac", **legendsettings)



def plotScenarioFlows(scen_results, scen_labels, settings, data,
                      percentage_relative_to=None, y_intercept=None, ylabel=None,
                      comp_labels=None, comp_titles=None, pop_labels=None,
                      link_labels=None, include_link_not_exclude=True, link_legend=None,
                      plot_inflows=True, plot_outflows=True, exclude_transfers=False, sum_total=False, sum_population=False,
                      colormappings=None, plot_observed_data=True, save_fig=False, fig_name=None, colors=None, legendsettings=None):
    """
    Line plots of flows across scenarios. Should be used for flows (instantaneous rates) only, such
    as incidence, deaths, notifications, etc.

    Params
        scen_results        list of results
        scen_labels         list of scenario labels, to be displayed
        settings            project settings
        data                project data

        other params as per plotFlows, including comp_labels, comp_titles, pop_labels



    """
    

    # close all remaining windows
    pl.close("all")
    # setup
    charac_specs = settings.charac_specs
    plotdict = settings.plot_settings
    year_inc = 5.  # TODO: move this to setting
    tvec = scen_results[0].sim_settings['tvec']  # TODO

    if legendsettings is None:
        legendsettings = {}

    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]

    yr_range = np.arange(start_year, end_year + 0.1, year_inc, dtype=int)


    if pop_labels is not None:
        plot_pids = getPIDs(scen_results[0], pop_labels)
    else:
        pop_labels = scen_results[0].pop_labels
        plot_pids = range(len(scen_results[0].m_pops))

    if link_legend is None: link_legend = dict()

    if comp_labels is None:
        logger.error("No compartments have been selected for flow-plots.")
        comp_labels = []
        return

    if comp_titles is not None and len(comp_titles) != len(comp_labels):
        logger.error("Flow-plot failure due to the number of compartment plot titles not matching the number of compartments to analyse.")

    """

    for (i,comp_label) in enumerate(comp_labels):

        rates, tvecs, all_labels = extractFlows(pop_labels=plot_pids,
    """

    for (i, comp_label) in enumerate(comp_labels):


        yvals = []
        tvals = []
        labels = []

        for (k, result_name) in enumerate(scen_results.keys()):

            result = scen_results[result_name]
            all_rates, all_tvecs, all_labels = extractFlows(pop_labels=plot_pids,
                                                            comp_label=comp_label,
                                                            results=result,
                                                            settings=settings,
                                                            tvec=tvec,
                                                            link_labels=link_labels,
                                                            include_link_not_exclude=include_link_not_exclude,
                                                            link_legend=link_legend,
                                                            plot_inflows=plot_inflows,
                                                            plot_outflows=plot_outflows,
                                                            sum_total=sum_total,
                                                            sum_population=sum_population,
                                                            exclude_transfers=exclude_transfers)

            # WARNING: Summing across rates so that multiple populations are combined. Not sure if this is fine for all intended cases.
            yv = sum(all_rates)  # all_rates[0]

            if percentage_relative_to is not None:
                yv /= percentage_relative_to
                yv *= 100.

            if len(all_tvecs) == 0 :
                continue

            yvals.append(yv)
            tvals.append(all_tvecs[0])
            labels.append(scen_labels[k])

        if comp_titles is not None:
            title_comp = comp_titles[i]
        else:
            title_comp = 'Compartment: "%s"' % settings.node_specs[comp_label]['name']


        # title_pop = '\nPopulation: "%s"' % pop_label

        final_dict = {'ylim' : 0,
          'xlabel':'Year',
          'ylabel': 'Number of People',
          'y_intercept': y_intercept,
          'x_ticks' : (yr_range, yr_range),
          'title': '',  # %s %s' % (title_comp,title_pop),
          'save_figname': '%s_ScenarioFlowComparision_%s_%s' % (fig_name, comp_label, all_labels[0])
          }
        if percentage_relative_to is not None:
            final_dict['ylabel'] = ylabel
            final_dict['ylim'] = [0, 105.]
        final_dict.update(plotdict)
        print colors
        _plotLine(ys=yvals, ts=tvals, labels=labels, colors=colors, save_fig=save_fig, reverse_order=True, smooth=True, **final_dict)

        if final_dict.has_key('legend_off') and final_dict['legend_off']:
        # Do this separately to main iteration so that previous figure are not corrupted
        # Note that colorlist may be different to colors, as it represents
        # classes of compartments

            separateLegend(labels=labels, colors=colors, fig_name=fig_name + "_LegendScenFlow", **legendsettings)




def isPlottableCharac(output_id, charac_specs):
    """
    Returns False if specified in cascade spreadsheet as 'n' or 'N' for plot characteristic.
    Else returns True
    """
    

    try:  # Better to try and ask for forgiveness than for permission ...
        if charac_specs[output_id]['plot_characteristic'].lower() == 'n':
            return False
    except:
        return True

def getPIDs(results, poplabels):
    """
    Takes in a list of poplabels and returns the corresponding PIDs

    TODO: this can be improved and made more efficient by either a better implementation of this
    look up, OR (better yet) by improving the data structures of mpops.
    """
    

    pids = []
    for poplabel in poplabels:
        for i, pop in enumerate(results.model.pops):
            if pop.label == poplabel:
                pids.append(i)
    return pids



def plotProjectResults(results, settings, data, title='',
                       colormappings=None, colorlabels=None, pop_colormappings=None,
                       pop_labels=None, plot_comp_labels=None, debug=False, plot_observed_data=True, save_fig=False, fig_name=None):
    """
    Plot all results associated with a project. By default, this is a disease cascade for
    each population, as well as characteristics of interest.
    If debug, then plot auxilliary plots
    associated with internal values of interest.

    Params:
        results         Results object
        charac_specs    ...
        title           project title, displayed in title
        colormapping
        pop_labels      list of population labels. Default: None, which selects all populations
        debug           boolean flag, indicating whether to plot internal variables if True
    """
    

    # close all remaining windows
    pl.close("all")
    # setup
    if pop_labels is None:
        pop_labels = results.pop_labels
    charac_specs = settings.charac_specs
    plotdict = settings.plot_settings

    # plot each disease cascade for every population
    plotPopulation(results=results, data=data, title=title, colormappings=colormappings, cat_labels=colorlabels, pop_labels=pop_labels, plot_observed_data=plot_observed_data, \
                    save_fig=save_fig, fig_name=fig_name, plotdict=plotdict, comp_labels=plot_comp_labels)

    # plot characteristics
    plotCharacteristic(results=results, settings=settings, pop_labels=pop_labels, data=data, colormappings=pop_colormappings,
                       plot_observed_data=plot_observed_data, save_fig=save_fig, fig_name=fig_name, plotdict=plotdict)

    # internal plotting
    if debug:
        plotAllOutflows(results)

def plotScenarios(scen_results, scen_labels, settings, data, plot_charac=None, pop_labels=None,
                  percentage_relative_to=None, y_intercept=None, ylabel=None,
                  colormappings=None, colors=None, plot_observed_data=False, save_fig=False, fig_name=None):
    """
    Line plots of characteristics across scenarios, including alive, total infected population, etc.

    Params
        scen_results        list of results
        scen_labels         list of scenario labels, to be displayed
        settings            project settings
        data                project data
        plot_charac         list of characteristics to be plotted. If None, then the default list from databook is used
        pop_labels          list of populations to be plotted. If None, then the default list from databook
        colormappings
        colors
        plot_observed_data
        save_fig
        fig_name
        percentage_relative_to    scalar (float) that results should be rescaled to, and shown as 100% of.
        y_bar                scalar that plots a dashed grey line at value at.
    """
    

    print "figname", fig_name

    # close all remaining windows
    pl.close("all")
    # setup
    charac_specs = settings.charac_specs
    plotdict = settings.plot_settings
    year_inc = 5.  # TODO: move this to setting
    tvec = scen_results[0].sim_settings['tvec']  # TODO
    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]
    yr_range = np.arange(start_year, end_year + 0.1, year_inc, dtype=int)

    if colors is None or len(colors) < len(scen_labels):
        colors = gridColorMap(len(scen_labels))
        logger.info("Plotting: setting color scheme to be default colormap, as not all lines had color assigned")


    if plot_charac is None:
        plot_charac = results.outputs.keys()
        plot_charac = [output_id for output_id in plot_charac if isPlottableCharac(output_id, charac_specs)]

    print "plot_charac --> ", plot_charac

    if pop_labels is not None:

        plot_pids = getPIDs(scen_results[0], pop_labels)  #####
    else:
        plot_pids = range(len(scen_results[0].m_pops))
        pop_labels = [pop.label for pop in scen_results[0].m_pops]

    # generate plots
    for (i, pid) in enumerate(plot_pids):
        pop_label = pop_labels[i]

        for charac in plot_charac:

            yvals = []
            labels = []
            observed_data = []
            yhat = []
            that = []

            for (ri, result_name) in enumerate(scen_results.keys()):
                result = scen_results[result_name]
                y_values_cur, _, _, dataobs = extractCharacteristic(results=result, charac_labels=[charac], charac_specs=charac_specs, data=data, pop_labels=[pop_label], plot_observed_data=plot_observed_data)

#                 y_values, labels, unit_tags, dataobs = extractCharacteristic(results, data, charac_specs, charac_labels=outputIDs, pop_labels=pop_labels, plot_observed_data=plot_observed_data, plot_total=plot_total)

                ys = y_values_cur[charac][:][0]

                if percentage_relative_to is not None:
                    # normalize, and also change ylabel
                    ys /= percentage_relative_to

                yvals.append(ys)
                labels.append(scen_labels[ri])

                if dataobs is not None:  # this can be improved
                    that_i, yhat_i = dataobs
                    yhat.append(yhat_i)
                    that.append(that_i)

            if plot_observed_data:
                pass  # TODO: include in future, but will require information as to which is the 'Current conditions' / 'BAU' scenario.

            unit_tag = ''
            if 'plot_percentage' in charac_specs[charac].keys():
                # ## we don't vals *= 100, as this is already done in extractCharacteristic()
                unit_tag = ' (%)'

            final_dict = {  # 'y_hat': yhat,
                  # 't_hat': that,
                  'ylim' : 0,
                  'unit_tag': unit_tag,
                  'xlabel':'Year',
                  'ylabel': charac_specs[charac]['name'] + unit_tag,
                  'y_intercept': y_intercept,
                  'x_ticks' : (yr_range, yr_range),
                  'smooth' : True,  # smooth plots for scenarios
                  'title': 'Scenario comparison:\n%s [%s]' % (charac_specs[charac]['name'], pop_label),
                  'save_figname': '%s_ScenarioComparision_%s_%s' % (fig_name, pop_label, charac_specs[charac]['name'])}

            if dataobs is not None:  # this can be improved

                final_dict['y_hat'] = yhat
                final_dict['t_hat'] = that
                print "Observed data = ", len(that)
                print that

            if percentage_relative_to is not None:
                final_dict['ylabel'] = ylabel

            print "/////////////"
            print plotdict
            final_dict.update(plotdict)
            print final_dict
            print "\\\\\\\\\\\\\\\\\\\\\|||"


            figure = _plotLine(ys=yvals, ts=np.tile(tvec, (len(labels), 1)), labels=labels, reverse_order=True,
                               save_fig=save_fig, fig_name=fig_name, colors=colors, **final_dict)  # , y_hat=[final_dict_cur['y_hat'][pid],final_dict_com['y_hat'][pid]], t_hat=[final_dict_cur['t_hat'][pid],final_dict_com['t_hat'][pid]])

    if final_dict.has_key('legend_off') and final_dict['legend_off']:
        # Do this separately to main iteration so that previous figure are not corrupted
        # Note that colorlist may be different to colors, as it represents
        # classes of compartments
        separateLegend(labels=labels, colors=colors, fig_name=fig_name + "_LegendScenarioComparison")

    return figure



def _plotLine(ys, ts, labels, colors=None, y_hat=[], t_hat=[],
             legendsettings=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, y_ticks=None, x_ticks=None,
             y_intercept=None, reverse_order=False, y_bounds=None, linestyles=None,
             smooth=False, symmetric=False, repeats=5, formatter=None,
             alpha=0.3, marker='o', s=40, facecolors='none', linewidth=3, zorder=10,
             save_fig=False, save_figname=None, legend_off=False,
             box_width=0.9, box_offset=0.0, **kwargs):
    """
    Plots multiple lines, with additional option of overlaying observed datapoints

    Params:
        ys        list of values for ys, with each entry corresponding to a line
        ts        list of values for xs, with each entry corresponding to a line
        labels    list of labels for each line
        colors    list of colors for each line
        y_intercept
        reverse_order
        y_bounds    list of array for each ys entry, with format of (tbound, ybound_min, ybound_ymax), thus can be specified independently of ts
        **kwargs    further keyword arguments, such as ylims, legend_off, edgecolors, etc.
    """
    


    if legendsettings is None: legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':1}

    if len(ys) == 0:
        # TODO: error to be fixed
        return

    if xlim is None: xlim = (ts[0][0], ts[0][-1])
    ymin_val = np.min(ys[0])
    indices = (ts[0] >= xlim[0]) * (ts[0] <= xlim[1])
    ymax_val = np.max(ys[0][indices])

    if colors is None or len(colors) < len(ys):
        colors = gridColorMap(len(ys))
        logger.info("Plotting: setting color scheme to be default colormap, as not all lines had color assigned")

    if linestyles is None or len(linestyles) < len(ys):
        try: linestyles = [kwargs['default_linestyle']] * len(ys)
        except: linestyles = ['-'] * len(ys)

        logger.info("Plotting: setting linestyles to be default value, as not all lines had styles assigned")


    fig, ax = pl.subplots()

    if y_intercept is not None:
        # TODO remove hardcoded values
        ax.hlines([y_intercept], np.min(ts[0]), np.max(ts[0]), colors='#AAAAAA', linewidth=0.75 * linewidth, linestyle='--')


    # plot ys, but reversed - and also reverse the labels (useful for scenarios, and optimizations):
    order_ys = range(len(ys))
    if reverse_order:
        logger.info("Reversing order of plot lines")
        order_ys = order_ys[::-1]  # surely there are more elegant ways to do this ...
        labels = labels[::-1]

    for k in order_ys:

        yval = ys[k]

        # if there are confidence bounds, plot using fill_between
        if y_bounds is not None:
            t_bound, y_min_bound, y_max_bound = zip(*y_bounds[k])[0] , zip(*y_bounds[k])[1] , zip(*y_bounds[k])[2]
            ax.fill_between(t_bound, y_min_bound, y_max_bound, facecolor=colors[k], alpha=alpha, linewidth=0.1, edgecolor=colors[k])

        # smooth line
        if smooth:
            yval = smoothfunc(yval, symmetric, repeats)

        # plot line
        ax.plot(ts[k], yval, c=colors[k], ls=linestyles[k])

        # identify ymin and ymax from line
        if np.min(yval) < ymin_val:
            ymin_val = np.min(yval)
        if np.max(yval[indices]) > ymax_val:
            ymax_val = np.max(yval[indices])

        # scatter data points
        if len(y_hat) > 0:
            try:
                if len(y_hat[k]) > 0:  # i.e. we've seen observable data
                    ax.scatter(t_hat[k], y_hat[k], marker=marker, edgecolors=colors[k], facecolors=facecolors, s=s, zorder=zorder, linewidth=linewidth)
                    # update min and max y based on observed datapoints
                    try:
                        if np.min(y_hat[k]) < ymin_val:
                            ymin_val = np.min(y_hat[k])
                    except:
                        pass
                    try:
                        if np.max(y_hat[k]) > ymax_val:
                            ymax_val = np.max(y_hat[k])
                    except:
                        pass
            except:
                logger.debug("No data plottable for index k=%i, data=\n" % k)
                logger.debug(y_hat)

    if formatter is not None:
        ax.yaxis.set_major_formatter(formatter)
    else:
        ax.get_yaxis().get_major_formatter().set_scientific(False)



    # set position
    box = ax.get_position()
    ax.set_position([box.x0 + box.width * box_offset, box.y0, box.width * box_width, box.height])

    # set title
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


    if not legend_off:
        ax.legend(labels, **legendsettings)

#     ymin = ax.get_ylim()[0]
    # Set the ymin to be halfway between the ymin_val and current ymin.
    # This seems to get rid of the worst of bad choices for ylabels[0] = -5 when the real ymin=0
#     tmp_val = (ymin + ymin_val) / 2.
#     ax.set_ylim(ymin=tmp_val, ymax=ymax_val)
    # Temporary choice (?) to enforce ylim = 0
    ax.set_ylim(ymin=0)

    # overwrite with specified choice
    if ylim is not None:
        ax.set_ylim(ylim)
    if xlim is not None:
        ax.set_xlim(xlim)

    ax.set_ylim(ymax=ax.get_ylim()[1] * 1.05)


    if x_ticks is not None:
        ax.set_xticks(x_ticks[0])
        ax.set_xticklabels(x_ticks[1])
    if y_ticks is not None:
        ax.set_yticks(y_ticks[0])
        ax.set_yticklabels(y_ticks[1])

    _turnOffBorder()
    if save_fig:
        fig.savefig('%s' % (save_figname))
        logger.info("Saved figure: '%s'" % save_figname)

    return fig


def _plotStackedCompartments(tvec, comps, labels=None, datapoints=None, title='', ylabel=None, xlabel=None, xlim=None, ymin=None, ylim=None, save_figname=None, legend_off=False,
                            save_fig=False, colors=None, catlabels=None, catcolors=None, year_inc=5, formatter=None,
                            marker='o', edgecolors='k', facecolors='none', s=40, zorder=10, linewidth=3, x_ticks=None, legendsettings={}, **kwargs):
    """
    Plot compartments in time.
    This creates a stacked plot of several compartments, over a given period.
    Observed data points can also be additionally plotted, overlaying the data

    Params:
        tvec        time period
        comps       compartment sizes
        labels      list of labels
        datapoints  observed datapoints, specified as a list of tuples
        **kwargs    further keyword arguments, such as ylims, legend_off, edgecolors, etc.
    """
    

    if colors is None or len(colors) != len(comps):
        if len(colors) != len(comps):
            logger.info("Plotting: setting color scheme to be default colormap, as not all compartments had color assigned")
        colors = gridColorMap(len(comps))

    # setup
    fig, ax = pl.subplots()
    bottom = 0 * tvec
    max_val = 0

    for (k, comp) in enumerate(comps):
        top = bottom + comp.vals
        lw = 0
        if save_fig:
            lw = 0.1

        ax.fill_between(tvec, bottom, top, facecolor=colors[k], alpha=1, lw=lw, edgecolor=colors[k])  # for some reason, lw=0 leads to no plot if we then use fig.savefig()
        reg, = ax.plot((0, 0), (0, 0), color=colors[k], linewidth=10)  # TODO fix this by using xlims and ylims appropriately

        bottom = dcp(top)

    max_val = max(top)

    if datapoints is not None:
        ts, ys = datapoints[0], datapoints[1]
        ax.scatter(ts, ys, marker=marker, edgecolors=edgecolors, facecolors=facecolors, s=s, zorder=zorder, linewidth=linewidth)
        if max(ys) > max_val:
            max_val = max(ys)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if formatter is not None:
        ax.yaxis.set_major_formatter(formatter)
    else:
        ax.get_yaxis().get_major_formatter().set_scientific(False)


    ax.set_ylim(ymax=max_val * 1.05)

    if ylim is not None:
        ax.set_ylim(ylim)
    elif ymin is not None:
        ax.set_ylim(ymin=ymin)

    if x_ticks is not None:
        ax.set_xticks(x_ticks[0])
        ax.set_xticklabels(x_ticks[1])
    else:
        ticks = np.arange(tvec[0], tvec[-1], year_inc)
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks)

    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim(tvec[0], tvec[-1])


    if not legend_off:
        ax.legend(labels, **legendsettings)

    _turnOffBorder()
    pl.suptitle('')

    if save_fig:
        fig.savefig('%s' % (save_figname))
        logger.info("Saved figure: '%s'" % save_figname)


def plotAllOutflows(results, num_subplots=5):
    """
    Visualise outflows for each compartment in each population as fractions of compartment size
    """
    

    mpops = results.model.pops
    sim_settings = results.sim_settings

    legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':2}

    pid = 0
    for pop in mpops:
        num_links = len(pop.links)
        colors = gridColorMap(num_links)
        cid = 0
        for comp in pop.comps:
            plot_id = cid % num_subplots
            patch_labels = []
            if plot_id == 0:
                fig, ax = pl.subplots(nrows=num_subplots, sharex=True)
                pl.suptitle('Population (%i): %s' % (pid, pop.label.title()))
                ax[num_subplots - 1].set_xlabel('Year')
            bottom = 0 * sim_settings['tvec']
            for link_tuple in comp.outlink_ids:
                link_id = link_tuple[-1]
                extra = dcp(pl.nan_to_num(pop.links[link_id].vals))
                top = bottom + extra
                ax[plot_id].fill_between(sim_settings['tvec'], bottom, top, facecolor=colors[link_id], alpha=1, lw=0)
                target_label = pop.comps[pop.links[link_id].index_to[1]].label
                target_pid = pop.links[link_id].index_to[0]
                if target_pid != pid: target_label += '(' + str(target_pid) + ')'
                patch_labels.append(target_label)
                ax[plot_id].plot((0, 0), (0, 0), color=colors[link_id], linewidth=10)
                bottom = dcp(top)

            ax[plot_id].set_ylabel(comp.label)
            ax[plot_id].set_xlim((sim_settings['tvec'][0], sim_settings['tvec'][-1]))
            ax[plot_id].set_ylim((0, max(top)))
            box = ax[plot_id].get_position()
            ax[plot_id].set_position([box.x0, box.y0, box.width * 0.75, box.height])
            ax[plot_id].legend(patch_labels, **legendsettings)

            cid += 1

        while cid % num_subplots != 0:
            plot_id = cid % num_subplots
            box = ax[plot_id].get_position()
            ax[plot_id].set_position([box.x0, box.y0, box.width * 0.75, box.height])
            cid += 1


        pid += 1




def plotScenarioBar (scen_results, scen_labels, settings, data, output_list=None, year=None, pop_labels=None, legendsettings=None,
                  colormappings=None, colors=None, plot_observed_data=False, save_fig=False, fig_name=None) :
    """


    """
    

    xlim = 3
    if len(scen_labels) > 3:
        xlim = len(scen_labels)

    plotdict = settings.plot_settings
    if plotdict is None:
        plotdict = {}

    if legendsettings is None:
        legendsettings = {}

    if year is None:
        pass  # set it to be the last year in the simulation

    # setup: determine colors to be used
    colors = []
    if colormappings is not None:
        colors_dict, cat_colors = getCategoryColors(colormappings, 'sequential')
        # reorder so that colors are same as expected for plotting the population
        for olabel in output_list:
            colors.append(colors_dict[olabel])

    values = [ [scen_results[rname].getValuesAt(output, year)[0] for output in output_list] for rname in scen_results.keys()]
    print values, '\n' * 10
    final_dict = plotdict

    final_dict2 = {'xlim': (0, xlim),
#                       'title':  title,
                  'ylabel': "",
                  'save_figname': fig_name + "_%s" % output}
    final_dict.update(final_dict2)


    _plotBars(values, labels=output, colors=colors, xlabels=scen_results.keys(), legendsettings=legendsettings,
              save_fig=save_fig, **final_dict)


    if final_dict.has_key('legend_off') and final_dict['legend_off']:
        # Do this separately to main iteration so that previous figure are not corrupted
        # Note that colorlist may be different to colors, as it can represent
        # classes of budgets
        separateLegend(labels=output, colors=cat_colors, fig_name=fig_name, reverse_order=True, **legendsettings)


def plotStackedBarOutputs(results, settings, year_list, output_list, output_labels=None,
                          xlabels=None, ylim=None, pop_labels=None,
                          title="", colormappings=None, save_fig=False, fig_name=None, legendsettings=None):
    """


    """
    

    xlim = 3
    if len(xlabels) > 3:
        xlim = len(xlabels)

    plotdict = settings.plot_settings
    if plotdict is None:
        plotdict = {}

    if legendsettings is None:
        legendsettings = {}


    # setup: determine colors to be used
    colors = []
    if colormappings is not None:
        colors_dict, cat_colors = getCategoryColors(colormappings, 'sequential')
        # reorder so that colors are same as expected for plotting the population
        for olabel in output_list:
            colors.append(colors_dict[olabel])

    if output_labels is None:
        output_labels = output_list

    # unfortunately we have to do it this way to ensure that the programs are all extracted in the same order
    values = [[results.getValuesAt(output_label, year, pop_labels=pop_labels)[0]  for output_label in output_list ]  for year in year_list]

#     print plotdict
#     final_dict = dcp(plotdict)

    plotdict['plot_stacked'] = True
    final_dict = {'xlim': (0, xlim),
                   'ylim': ylim,
                  'title':  '',
                  'ylabel': "",
                  'plot_stacked': True,
                  'barwidth': 0.8,
                  'save_figname': fig_name}
#     final_dict.update(final_dict2)

#     final_dict.update(plotdict)
    plotdict.update(final_dict)
#    print plotdict
    _plotBars(values, labels=output_labels, colors=colors, xlabels=xlabels, # legendsettings=legendsettings,
              save_fig=save_fig, **plotdict)


    if plotdict.has_key('legend_off') and plotdict['legend_off']:
        # Do this separately to main iteration so that previous figure are not corrupted
        # Note that colorlist may be different to colors, as it can represent
        # classes of budgets
        separateLegend(labels=output_labels, colors=cat_colors, fig_name=fig_name, reverse_order=True, **legendsettings)


def plotValueBars(results, settings, pop_labels, label, years, title=None, y_intercept=None,
                    additional_labels=None, additional_initial_alpha=0.5, colors=None,
                    colormappings=None, cat_labels=None, use_full_labels=False, full_labels=None,
                    hatchstyles=None,
                    save_fig=False, fig_name=None, legendsettings=None):
    """


    hatch='//'
    """
    

    plotdict = settings.plot_settings
    num_bars = len(pop_labels)
    plotdict['xlim'] = (-0.25, num_bars)
    bar_width = 0.8 / len(years)
    plotdict['barwidth'] = bar_width
    plotdict['bar_offset'] = 0.2
    plotdict['plot_stacked'] = False
    plotdict['ylabel'] = settings.linkpar_specs[label]['name']

    alphas = list(np.linspace(1., 0.6, len(years)))

#     if hatchstyles is not None:
#         # convert from odict (key: style) to odict (key: population)
#         hatchstyles = getLinemapping(hatchstyles)
#         hatchstyles = [hatchstyles[pop] for pop in pop_labels]

    if colors is not None and len(colors) >= len(pop_labels):
        pass  # colors as defined in the args should be used as is
    elif colormappings is not None and colors is None:
        colors = []
        colors_dict, cat_colors = getCategoryColors(colormappings, 'sequential')
        # reorder so that colors are same as expected for plotting the population
        for (j, pop_label) in enumerate(pop_labels):
            colors.append(colors_dict[pop_label])
    else:
        colors = None


    values = []
    inds = []

    for (i, year_init) in enumerate(years):
        pval = []
        for pop in pop_labels:
            pval.append(results.getValuesAt(label, year_init, year_init + 1., pop, integrated=True)[0])
        values.append(pval)

        inds.append(list(np.arange(num_bars) + np.ones(num_bars) * i * bar_width))

    if colors is not None:
        colors = [colors] * len(years)

    plotdict['save_figname'] = fig_name
#     plotdict['']

    _plotBars(values, pop_labels, inds=inds, colors=colors, xlabels=pop_labels, alphas=alphas, save_fig=save_fig, **plotdict)

def plotCareCascade(results, settings, pop_labels, labels, years, title="", normalize=False, y_intercept=None,
                    additional_labels=None, additional_initial_alpha=0.5,
                    colormappings=None, cat_labels=None, use_full_labels=False, full_labels=None,
                    save_fig=False, fig_name=None, legendsettings=None):
    """

    """
    

    xlabels = labels.keys()
    plotdict = dcp(settings.plot_settings)

    yticks = None
    if plotdict is None:
        plotdict = {}

    plotdict['bar_width'] = 0.8 / len(years)
    xlim = (-0.25, len(xlabels))

    # setup: determine colors to be used
    colors = []
    if colormappings is not None:
        colors_dict, cat_colors = getCategoryColors(colormappings, 'sequential')
        # reorder so that colors are same as expected for plotting the population
        for (j, pop_label) in enumerate(pop_labels):
            colors.append(colors_dict[pop_label])


    cat_labels = pop_labels.keys()

    for (i, year_init) in enumerate(years):
        values = []
        inds = []
        for (j, lab_key) in enumerate(labels.keys()):
            lab_set = labels[lab_key]
            pop_vals = []
            for pop_label, pop_set in pop_labels.iteritems():
                count = 0
                for lab in lab_set:
                    count += results.getValuesAt(lab, year_init, year_init + 1., pop_set, integrated=True)[0]
                pop_vals.append(count)
            values.append(pop_vals)
            inds.append(range(len(xlabels)))

        if normalize:
            # determine the max for the first state in the cascade
            init_state = np.sum(values[0])
            y_ticks_pos = np.arange(0.1, 1.01, 0.1) * init_state
            y_ticklabels = ['%g' % tick for tick in np.arange(10, 101, 10)]
            yticks = [y_ticks_pos, y_ticklabels]

            if y_intercept is not None:
                # convert to real numbers
                y_intercept = [yis * init_state for yis in y_intercept]


#        print values
#        print plotdict
#        print colors
#        print "--> convert"
        new_val = np.array(values)
        new_val.transpose()

#         plotdict['bar_offset'] = i * plotdict['bar_width']
        plotdict['xlim'] = xlim
        plotdict['save_figname'] = fig_name + "_%g" % year_init

        plotdict['plot_stacked'] = True
        plotdict['barwidth'] = 0.8
#        print values
#        print pop_labels.keys()
        _plotBars(new_val, pop_labels.keys(), xlabels=xlabels,
                  colors=colors, y_intercept=y_intercept, yticks=yticks,
                  legendsettings=legendsettings, save_fig=save_fig, **plotdict)


    if plotdict.has_key('legend_off') and plotdict['legend_off']:
        # Do this separately to main iteration so that previous figure are not corrupted
        # Note that colorlist may be different to colors, as it can represent
        # classes of budgets
        # reverse legend order so that it matches top<->bottom of stacked bars
        if use_full_labels:
            if legendsettings is None:
                legendsettings = {}
            separateLegend(labels=full_labels, colors=colors, fig_name=fig_name + "_legend", reverse_order=True, **legendsettings)
        else:
            print cat_labels
            separateLegend(labels=cat_labels, colors=cat_colors, fig_name=fig_name + "_legend", reverse_order=True,)

