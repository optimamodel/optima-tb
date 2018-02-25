# %% Imports
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

"""
Plotting library for Optima TB. 

Examples:

    # setup lines and colors:
    pop_colors = odict() # define colormappings per populations
    pop_colors['ocean'] = ['0-17', '18-59', '60+']
    pop_colors['#41A0BF'] = ['18-59 [D]']
    pop_colors['#BFDFE9'] = ['60+ [D]']  # will be dashed line for PLHIV
    pop_colors['#F1A94E'] = ['Prisoners', 'Prisoners [D]']
    pop_clabels = ['0-17', '18-59', '60+', '18-59 [D]', '60+ [D]', 'Prisoners', 'Prisoners [D]']
    linestyles = odict() # define linestyle per populations
    linestyles['-'] = ['0-17', '18-59', '60+', 'Prisoners']
    linestyles['--'] = ['18-59 [D]', '60+ [D]', 'Prisoners [D]']
    gen_pops = ['0-17', '18-59', '60+']


    ## Plotting line trends
    # plot simple results
    plotResult(proj, results, output_labels=outputids,
               colormappings=pop_colors, linestyles=linestyles, save_fig=save_results, fig_name=filename + "perPop")
    # plot total for all populations
    plotResult(proj, results, output_labels=outputids,
               plot_total=True, save_fig=save_results, fig_name=filename + "perTotal")
    # plot stacked for general populations only
    plotResult(proj, results, output_labels=outputids, plot_type='stacked', pop_labels=gen_pops,
               colormappings=pop_colors, linestyles=linestyles, save_fig=save_results, fig_name=filename + "perGenPopStacked")
    
    # plot without observed data
    plotResult(proj, results, output_labels=outputids, plot_observed_data=False,
               colormappings=pop_colors, linestyles=linestyles, save_fig=save_results, fig_name=filename + "perPopNoData")
    
    # plot relative
    plot_rel = {"year": 2015.} # plot relative to 2015 values
    y_intercepts = [5, 10, 50] # intercepts at 5%, 10%, and 50%
    plotResult(proj, results, output_labels=outputids, plot_observed_data=False, plot_relative=plot_rel, y_intercept=y_intercepts,
               colormappings=pop_colors, linestyles=linestyles, save_fig=save_results, fig_name=filename + "perPopNoDataRelative")
    plotResult(proj, results, output_labels=outputids, plot_observed_data=False, plot_relative=plot_rel, y_intercept=y_intercepts,
               plot_total=True, save_fig=save_results, fig_name=filename + "perTotalNoDataRelative")
    plotResult(proj, results, output_labels=outputids, plot_type='stacked', pop_labels=pops, plot_observed_data=False, plot_relative=plot_rel,
                colormappings=pop_colors, linestyles=linestyles, save_fig=save_results, fig_name=filename + "perPopStackedNoDataRelative")


    See plotting methods for more examples:
        plotResult
        plotCompareResults
        plotYearsBar
        plotCascade
        plotCompsBar
        plotCompareResultsBar
        plotPopulationCrossSection
        plotPopulationCrossSectionBar
"""

CMAPS = ['Blues', 'Purples', 'Reds', 'Oranges', 'Greys', 'Greens', ]
        # Suscep  #Latent  # SP TB    #SN TB    # Dead   # Recovered

PLOTTYPE_LINE = "line"
PLOTTYPE_STACKED = "stacked"
PLOTTYPE_BAR = "bar"
PLOTTYPE_BARSTACKED = "barstacked"

COMPARETYPE_RESULT = 'result'
COMPARETYPE_POP = 'pop'
COMPARETYPE_VALUE = 'output'
COMPARETYPE_YEAR = 'year'
COMPARETYPE_CASCADE = 'cascade'

def gridColorMap(ncolors=10, limits=None, nsteps=10, asarray=False, doplot=False, newwindow=True):
    '''
    Create a qualitative colormap by assigning points according to the maximum pairwise distance in the
    color cube. Basically, the algorithm generates n points that are maximally uniformly spaced in the
    [R, G, B] color cube.
    
    Arguments:
        ncolors: the number of colors to create
        limits: how close to the edges of the cube to make colors (to avoid white and black)
        nsteps: the discretization of the color cube (e.g. 10 = 10 units per side = 1000 points total)
        asarray: whether to return the colors as an array instead of as a list of tuples
        doplot: whether or not to plot the color cube itself
        newwindow: if doplot=True, whether to use a new window

    Usage example:
        from pylab import *
        from colortools import gridcolormap
        ncolors = 10
        piedata = rand(ncolors)
        colors = gridcolormap(ncolors)
        figure()
        pie(piedata, colors=colors)
        gridcolormap(ncolors, doplot=True)
        show()

    Version: 2015dec29 (cliffk)
    '''

    # # Imports
    from numpy import linspace, meshgrid, array, transpose, inf, zeros, argmax, minimum
    from numpy.linalg import norm

    # Steal colorbrewer colors for small numbers of colors
    colorbrewercolors = array([
    [27, 158, 119],
    [217, 95, 2],
    [117, 112, 179],
    [231, 41, 138],
    [255, 127, 0],
    [200, 200, 51],  # Was too bright yellow
    [166, 86, 40],
    [247, 129, 191],
    [153, 153, 153],
    ]) / 255.

    if ncolors <= len(colorbrewercolors):
        colors = colorbrewercolors[:ncolors]

    else:  # Too many colors, calculate instead
        # # Calculate sliding limits if none provided
        if limits is None:
            colorrange = 1 - 1 / float(ncolors ** 0.5)
            limits = [0.5 - colorrange / 2, 0.5 + colorrange / 2]

        # # Calculate primitives and dot locations
        primitive = linspace(limits[0], limits[1], nsteps)  # Define primitive color vector
        x, y, z = meshgrid(primitive, primitive, primitive)  # Create grid of all possible points
        dots = transpose(array([x.flatten(), y.flatten(), z.flatten()]))  # Flatten into an array of dots
        ndots = nsteps ** 3  # Calculate the number of dots
        indices = [0]  # Initialize the array

        # # Calculate the distances
        for pt in range(ncolors - 1):  # Loop over each point
            totaldistances = inf + zeros(ndots)  # Initialize distances
            for ind in indices:  # Loop over each existing point
                rgbdistances = dots - dots[ind]  # Calculate the distance in RGB space
                totaldistances = minimum(totaldistances, norm(rgbdistances, axis=1))  # Calculate the minimum Euclidean distance
            maxindex = argmax(totaldistances)  # Find the point that maximizes the minimum distance
            indices.append(maxindex)  # Append this index

        colors = dots[indices, :]

    # # Wrap up: optionally turn into a list of tuples
    if asarray:
        output = colors
    else:
        output = []
        for i in range(ncolors): output.append(tuple(colors[i, :]))  # Gather output

    # # For plotting
    if doplot:
        from mpl_toolkits.mplot3d import Axes3D  # analysis:ignore
        from pylab import figure, gca
        if newwindow:
            fig = figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = gca(projection='3d')
        ax.scatter(colors[:, 0], colors[:, 1], colors[:, 2], c=output, s=200, depthshade=False)
        ax.set_xlabel('R')
        ax.set_ylabel('G')
        ax.set_zlabel('B')
        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))
        ax.set_zlim((0, 1))
        ax.grid(False)

    return output


def _getColormapRange(cmap, ncols=10, order='alternate'):
    """
    Returns a list of colors values, according to colormap
    
    Params:
        cmap    colormap to be used, specified as name of matplotlib.cmap (String)
        ncols   number of colors to be returned (int)
        order   order in which colors should be chosen. 
                Possible values are 'alternate', 'alternate3','sequential' or 'random'
        
    Usage:
        cmap = 'jet'
        n = 6
        clist = _getColormapRange(cmap,n,'random')
        for i in range(n):
            pylab.plot(xs[i],ys[i],c=clist[i])
    """
    color_norm = colors.Normalize(vmin=-1.5, vmax=ncols + 0.5)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap=cmap)
    order_map = range(ncols)
    if order == 'random':
        shuffle(order_map)
    elif order == 'alternate':
        if ncols < 2:
            return _getColormapRange(cmap, ncols, 'sequential')
        a = order_map[ncols / 2:]
        b = order_map[:ncols / 2]
        order_map[::2] = a
        order_map[1::2] = b
    elif order == 'alternate3':
        if ncols < 3:
            return _getColormapRange(cmap, ncols, 'alternate')
        third = ncols / 3
        a = order_map[2 * third:]
        b = order_map[third:2 * third]
        c = order_map[:third]
        order_map[::3] = a
        order_map[1::3] = b
        order_map[2::3] = c
    else:  # order == 'sequential', so leave in order but add buffer at beginning and end of colormap
        order_map = range(ncols + 2)
        return [scalar_map.to_rgba(index) for index in order_map[1:-1]]

    return [scalar_map.to_rgba(index) for index in order_map]


def getCategoryColors(category_list, order='alternate'):
    """
    For an ordered dictionary of category list, return a mapping of compartment labels to colors (rgba). 
    
    Note that using colormappings will index by compartment keys, so the color maps do not have to be 
    defined in sequence of their appearance within a population's compartment list. 
    The only constraint is that all plottable compartments must be assigned a colour. 
    If a compartment is included in more than one colormap list, its cmap will be set as last list it appeared on. 
    
    Params:
        category_list    odict of list
        order            defines how neighbouring values within a list should be defined. Values are
                            'alternate' : 
                            'alternate3' : 
                            'random'    :
                            'sequential' : 
        
    Returns:
        an odict of (comp_label, rgba) items and a category color (as a representative color for that group)
        
    Usage:
        cat_list = odict()
        cat_list['Blues'] =    ['sus','vac'],
        cat_list['Reds'] =     ['inf','vir'],
        cat_list['Greens'] =   ['treat','rec']
        col_list = self.getCategoryColors(cat_list)
        print col_list['inf'] 
        > (0.96, 0.5,  0.5, 1.)
        print col_list[-1] #rec
        > (0.4, 0.96, 0.4, 1.)
    
    """
    col_list = odict()
    cat_colors = []
    for k, v in category_list.iteritems():

        if isinstance(k, str) and k.startswith('#'):
            # is a hex format
            tmp_list = [k] * len(v)
        elif isinstance(k, str) :
            # must be a colormap
            tmp_list = _getColormapRange(k, len(v), order)
        else:
            # else, unknown
            raise OptimaException('Unknown color format: ' + k)
        for i, label in enumerate(v):
            col_list[label] = tmp_list[i]
        cat_colors.append(tmp_list[0])  # add a representative color
#     print col_list
    return col_list, cat_colors

def getLinemapping(linestyle_dict):
    """
    Generates a dict that maps the selector (population or output_label) to the linestyle, when supplied
    with a dict where the keys are the linestyles.
    
    Params:
        linestyle_dict  dictionary of (key, values) where 
                                keys = linestyle 
                                values = list of selectors (pops or output_label) 
        
    Returns
        linedict        dictionary of (key, values) where 
                                keys = selector
                                values = linestyle
    """
    linedict = odict()
    for k, v in linestyle_dict.iteritems():
        tmp_list = [k] * len(v)
        for i, label in enumerate(v):
            linedict[label] = tmp_list[i]
    return linedict

def getHatchmapping(linestyles, labels):
    """
    Generates a hatch mapping from linestyles. 
    
    Params:
        linestyles    dict of (k,v) = (selector, linestyles)
         
    Returns:
        hatches        dict of (k, v) = (selector, hatchstyles)
    """
    hatches = {}
    if linestyles is not None:
        for (i, lab) in enumerate(labels):
            if linestyles[i] == '-':
                hatches[lab] = None
            elif linestyles[i] == '--':
                hatches[lab] = '///'
            elif linestyles[i] == '.':
                hatches[lab] = '.'
            else:
                logging.debug("Unknown linestyle used --> setting to unmapped hatch content")
                hatches[lab] = '+'
    else:
        for lab in labels:
            hatches[lab] = None
    return hatches


def setupStylings(colormappings, colors, linestyles, series_labels, plotdict):
    """
    
    """
    # generic setup for colors, line and hatches
    if colormappings is not None:
        colors = []
        colors_dict, cat_colors = getCategoryColors(colormappings, 'sequential')
        # extract colors so that they are in the same order as expected for plotting the population
        try:
            for (j, pop_label) in enumerate(series_labels):
                colors.append(colors_dict[pop_label])
        except Exception as e:
            colors = gridColorMap(len(series_labels))
            cat_colors = colors
            logging.warn("Had to set default colormap as labels didn't match colormappings")
            logging.warn(e)
    elif colors is not None and len(colors) >= len(series_labels):
        # colors as defined in the args should be used as is
        cat_colors = colors
    else:
        colors = gridColorMap(len(series_labels))
        cat_colors = colors
        logger.info("Plotting: setting color scheme to be default colormap, as not all lines had color assigned")

    if linestyles is None:
        linestyles = [plotdict['default_linestyle']] * len(series_labels)  #default value
    else:
        # convert from odict (key: style) to odict (key: population)
        linestyles = getLinemapping(linestyles)
        try:
            linestyles = [linestyles[pop] for pop in series_labels] # extract used linestyles from dict into list
        except Exception as e:
            linestyles = [plotdict['default_linestyle']] * len(series_labels)
            logging.warn("Had to use default linestyle as labels didn't match linestyles")

    hatches = getHatchmapping(linestyles, series_labels)
    hatches = [hatches[pop] for pop in series_labels] # extract used hatches from dict into list

    return colors, linestyles, hatches, cat_colors



def separateLegend(labels, colors, fig_name, reverse_order=False, linestyles=None, **legendsettings):
    """
    
    
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    hatches = getHatchmapping(linestyles, labels)

    if reverse_order:
        labels = labels[::-1]
        colors = colors[::-1]


    fig = plt.figure() # figsize=(5, 5))  # silly big
    patches = [  mpatches.Patch(color=color, label=label, ec='white', hatch=hatches[label]) for label, color in zip(labels, colors)]
    legendsettings['loc'] = 'center'
    legendsettings['frameon'] = False
    legendsettings['bbox_to_anchor'] = None
    fig.legend(patches, labels, **legendsettings)
    plt.savefig("%s_legend" % fig_name)  # ,bbox_inches='tight')

def getYearLabels(timeperiods):
    """
    Simple function for mapping years to string equivalent
    """
    return ['%g' % yr for yr in timeperiods]

def _turnOffBorder():
    """
    Turns off top and right borders, leaving only bottom and left borders on.
    """
    pl.gca().spines['right'].set_color('none')
    pl.gca().spines['top'].set_color('none')
    pl.gca().xaxis.set_ticks_position('bottom')
    pl.gca().yaxis.set_ticks_position('left')


def plotResult(proj, result, output_labels=None, pop_labels=None,
               plot_total=False, plot_type=None, plot_relative=None,
               plot_observed_data=True, observed_data_label=None,
               plot_ybounds=None,
               colormappings=None, colors=None, linestyles=None,
               title=None, save_fig=False, fig_name=None, **kwargs):
    """
    Plots either characteristics, compartment size, or flow rate, as line. 
    
    Supports plotting of:
        - multiple populations
        - plotting total across multiple populations
        - plotting observed datapoints (where applicable)
    
    If neither colormappings or colors is specified, default color list is used from gridColorMap()
    
    Params:
        proj            project object, containing plotting settings and observed data points
        result          result object 
        output_labels   list of compartment labels, flow rate labels, characteristics
        pop_labels      populations to be plotted. Default (None) plots all populations.
        plot_total      plot total across populations
        plot_type       plot type. Currently supported = line | stacked   (see PLOTTYPE_ fields)
                        if None, default value = line
        plot_relative   whether to normalize plot as a percentage. Dict with key, value indicating
                            ("year" : year value), calculated for each population or total population if plot_total=True
                            ("value": value), a numeric value
        plot_observed_data    add observed datapoints as scatter plot, if corresponding datapoints exist. See note.
        observed_data_label   dict of mappings of what datapoints should be used. See note.
        colormappings   colormappings that should be used to generate colors for populations. Supercedes colors. 
                        Format: odict with color / colormap as key, value = list of populations for corresponding key 
        colors          list of colors that should be used for population. Superceded by colormappings.
        linestyles      odict of (k,v): (linestypes, list of populations). Uses default linestyle if not supplied.
        title           title for plot
        save_fig        boolean flag, whether to save plot
        fig_name        if plot is saved, filename
        kwargs          Other params acceptable to innerPlotTrend i.e. 
                            y_intercept
        
    Returns:
        fig             fig handle
        
    Replaces:
        plotCharacteristic
        plotPopulationFlows
    
    Example:
        output_labels = ['lt_inf', # characteristic
                         'lt_prev', # characteristic which is a percentage
                         'spdd', # compartment
                         'infac_per100K'] # flow rate
        subsetPop = ['15-64', '15-64 HIV+']
        # plot per population for all populations, using colormappings
        plotResult(proj, results, output_labels=output_labels, save_fig=save_results, fig_name=filename + '_PerPopulation',
                   colormappings=pop_colors, linestyles=linestyles)
        # plot stacked populations
        plotResult(proj, results, output_labels=output_labels, save_fig=save_results, fig_name=filename + 'Stacked',
                   colormappings=pop_colors, plot_type='stacked')
        # plot only a subset of populations
        plotResult(proj, results, output_labels=output_labels, pop_labels=subsetPop, 
                    save_fig=save_results, fig_name=filename + '_15-64Only', colormappings=pop_colors, linestyles=linestyles)
        # plot total across all populations
        plotResult(proj, results, output_labels=output_labels, save_fig=save_results, fig_name=filename + '_Total', plot_total=True)
        # plot total across subset of populations
        plotResult(proj, results, output_labels=output_labels, pop_labels=subsetPop, 
                    save_fig=save_results, fig_name=filename + '_TotalFor15-64', plot_total=True)
           
    """
    if output_labels is None:
        logging.error("No output label specified for plotting")
        return None

    for out_label in output_labels:
        fig = innerPlotTrend(proj, [result], [out_label], compare_type=COMPARETYPE_POP, pop_labels=pop_labels,
                             plot_total=plot_total, plot_type=plot_type, plot_relative=plot_relative,
                             plot_observed_data=plot_observed_data, observed_data_label=observed_data_label,
                             plot_ybounds=plot_ybounds,
                             colormappings=colormappings, colors=colors, linestyles=linestyles,
                             title=title, save_fig=save_fig, fig_name=fig_name, **kwargs)
    return fig # return last fig handle

def plotCompareResults(proj, resultset, output_labels, pop_labels=None,
                       plot_total=False, plot_observed_data=True, observed_data_label=None,
                       colormappings=None, colors=None, linestyles=None,
                       plot_relative=None, y_intercept=None,
                       title=None, save_fig=False, fig_name='', **kwargs):
    """
    Plots either characteristics, compartment size, or flow rate, as across a result set. 
    
    Supports plotting of:
        - plotting total across multiple populations
        - plotting observed datapoints (where applicable)
    
    If neither colormappings or colors is specified, default color list is used from gridColorMap()
    
    Params:
        proj            project object, containing plotting settings and observed data points
        resultset       resultset object 
        output_labels   list of compartment labels, flow rate labels, characteristics
        pop_labels      populations to be plotted. Default (None) plots all populations.
        plot_total      plot total across populations
        plot_observed_data    add observed datapoints as scatter plot, if corresponding datapoints exist. See note.
        observed_data_label_dict   dict of mappings of what 
        colormappings   colormappings that should be used to generate colors for populations. Supercedes colors. 
                        Format: odict with color / colormap as key, value = list of populations for corresponding key 
        colors          list of colors that should be used for population. Superceded by colormappings.
        linestyles      odict of (k,v): (linestypes, list of populations). Uses default linestyle if not supplied.
        title           title for plot
        save_fig        boolean flag, whether to save plot
        fig_name        if plot is saved, filename
        kwargs          Other parameters accepted by innerPlotTrend
                        e.g. y_intercept
        
    Replaces:
        plotScenario
        plotScenarioFlows
    
    Example:
        output_labels = ['lt_inf', # characteristic
                         'lt_prev', # characteristic which is a percentage
                         'spdd', # compartment
                         'infac_per100K'] # flow rate
        subsetPop = ['15-64', '15-64 HIV+']
        resultset = proj.runScenarios()
        # plot per population for all populations, using colormappings
        plotCompareResults(proj, resultset, output_labels=output_labels, save_fig=save_results, fig_name=filename + '_PerPopulation',
                   colormappings=pop_colors, linestyles=linestyles)
        # plot only a subset of populations
        plotCompareResults(proj, resultset, output_labels=output_labels, pop_labels=subsetPop, 
                    save_fig=save_results, fig_name=filename + '_15-64Only', colormappings=pop_colors, linestyles=linestyles)
        # plot total across all populations
        plotCompareResults(proj, resultset, output_labels=output_labels, save_fig=save_results, fig_name=filename + '_Total', plot_total=True)
        # plot total across subset of populations
        plotCompareResults(proj, resultset, output_labels=output_labels, pop_labels=subsetPop, 
                    save_fig=save_results, fig_name=filename + '_TotalFor15-64', plot_total=True)

    """
    for out_label in output_labels:
        if plot_total:
            fig = innerPlotTrend(proj, resultset, [out_label], compare_type=COMPARETYPE_RESULT, pop_labels=pop_labels, plot_total=True,
                   plot_observed_data=plot_observed_data, observed_data_label=observed_data_label,
                   colormappings=colormappings, colors=colors, linestyles=linestyles, plot_relative=plot_relative, y_intercept=y_intercept,
                   title=title, save_fig=save_fig, fig_name=fig_name, **kwargs)
        else:
            logger.info("Plotting result set per population group")
            if pop_labels is None:
                pop_labels = getPops(resultset[0])
            # plot for each population
            for pop in pop_labels:
                fig = innerPlotTrend(proj, resultset, [out_label], compare_type=COMPARETYPE_RESULT, pop_labels=[pop], plot_total=True,
                   plot_observed_data=plot_observed_data, observed_data_label=observed_data_label,
                   plot_relative=plot_relative, y_intercept=y_intercept,
                   colormappings=colormappings, colors=colors, linestyles=linestyles,
                   title=title, save_fig=save_fig, fig_name=fig_name + '_%s' % pop, **kwargs)
            logger.info("Created multiple plots for plotCompareResults for multiple populations. Returning last plot created")
    return fig # return last


def plotYearsBar(proj, result, output_labels, pop_labels=None, year_periods=None,
               plot_total=False, plot_type=None, ylabel=None,
               colormappings=None, colors=None, linestyles=None,
               title=None, save_fig=False, fig_name=None, **kwargs):
    """
    Plot a collection of outlabels per population set, over a range of years
    
    Params:
        proj
        result
        output_labels     list of items compartments against which to plot
        pop_labels
        year_periods
        plot_total
        plot_type
        y_label
        colormappings
        colors
        linestyles
        title
        save_fig
        
            
    Examples:
        # number of new infections vs reinfections every year from 2005 to 2015
        new_inf = ["new_cases", "retreat_cases", "failed_cases", "tb_death_cases"]
        colors = ['#121675', '#f72d00', '#20a322', '#777777']
        tlinestyles = odict()
        tlinestyles['-'] = ['new_cases']
        tlinestyles['--'] = ['retreat_cases']
        tlinestyles['.'] = ['failed_cases']
        tlinestyles['*'] = ['tb_death_cases']
        year_periods = range(2005, 2016)
        # as absolute numbers for each year 
        plotYearsBar(proj, results, year_periods=year_periods, output_labels=new_inf,
                     pop_labels=None, plot_total=False, colors=colors, linestyles=tlinestyles,
                     save_fig=save_results, fig_name=filename + "compareCompsNewInfRetreat")
        # normalized for each year
        plot_rel = ("normalized")
        y_intercepts = [50, 25]
        plotYearsBar(proj, results, year_periods=year_periods, output_labels=new_inf,
                     y_intercept=y_intercepts, plot_relative=plot_rel,
                     pop_labels=None, plot_total=False, colors=colors, linestyles=tlinestyles,
                     save_fig=save_results, fig_name=filename + "compareCompsNewInfRetreatRelative")

    """
    fig = innerPlotBar(proj, [result], year_periods=year_periods, output_labels=output_labels,
                       compare_type=COMPARETYPE_YEAR, pop_labels=pop_labels, ylabel=ylabel,
                       colormappings=colormappings, colors=colors, linestyles=linestyles,
                       title=title, save_fig=save_fig, fig_name=fig_name, **kwargs)
    return fig

def plotCascade(proj, result, output_labels, pop_labels=None, year_periods=None,
                plot_total=False, plot_relative=None,
               colormappings=None, colors=None, linestyles=None, y_intercept=None,
               title=None, save_fig=False, fig_name=None, **kwargs):
    """
    TODO doc
    """
    fig = innerPlotBar(proj, [result], year_periods=year_periods, output_labels=output_labels,
                       compare_type=COMPARETYPE_CASCADE, pop_labels=pop_labels, ylabel="Number of cases",
                       colormappings=colormappings, colors=colors, linestyles=linestyles,
                       plot_relative=plot_relative, y_intercept=y_intercept, xlim=(0, len(output_labels)),
                       title=title, save_fig=save_fig, fig_name=fig_name, **kwargs)


def plotCompsBar(proj, result, output_labels, pop_labels=None, year_periods=None,
                plot_total=False, plot_relative=None,
#                 plot_type=None,
#                plot_observed_data=True, observed_data_label=None,
               colormappings=None, colors=None, linestyles=None, y_intercept=None,
               title=None, save_fig=False, fig_name=None, **kwargs):
    """
    Plots a collection of different compartments in different bars.
    
    Params:
        proj
        result
        output_labels    odict with (key, value) = (string, list), where each value-list is a list of values
        pop_labels
        year_period      a value or tuple
        
    Examples:
        # in and out flows to a compartment
    
    """
    ylabel = 'Number of cases'
    if plot_total:
        fig = innerPlotBar(proj, [result], year_periods=year_periods, output_labels=output_labels,
                       compare_type=COMPARETYPE_VALUE, pop_labels=pop_labels, ylabel=ylabel,
                       colormappings=colormappings, colors=colors, linestyles=linestyles, plot_relative=plot_relative,
                       title=title, save_fig=save_fig, fig_name=fig_name, **kwargs)
    else:
        if pop_labels is None:
            pop_labels = getPops(result)
        for pop in pop_labels:
            fig = innerPlotBar(proj, [result], year_periods=year_periods, output_labels=output_labels,
                       compare_type=COMPARETYPE_VALUE, pop_labels=[pop], ylabel=ylabel,
                       colormappings=colormappings, colors=colors, linestyles=linestyles, plot_relative=plot_relative,
                       title=title, save_fig=save_fig, fig_name=fig_name + "_%s" % pop, **kwargs)

    return fig



def plotCompareResultsBar(proj, resultset, output_labels, pop_labels=None, year_periods=None,
                       plot_total=False, plot_observed_data=True, observed_data_label=None,
                       colormappings=None, colors=None, linestyles=None, y_intercept=None,
                       title=None, save_fig=False, fig_name=None, **kwargs):
    """
    
    """
    for out_label in output_labels:
        fig = innerPlotBar(proj, resultset, year_periods=year_periods, output_labels=output_labels,
                           compare_type=COMPARETYPE_RESULT, plot_total=plot_total,
                           pop_labels=pop_labels,
                           colormappings=colormappings, colors=colors, linestyles=linestyles,
                           title=title, save_fig=save_fig , fig_name=fig_name + "_%s" % (out_label), **kwargs)
    return fig


def plotCompareCascade(proj, resultset, output_labels, pop_labels=None, year_periods=None,
                       plot_total=False, plot_observed_data=True, observed_data_label=None,
                       colormappings=None, colors=None, linestyles=None, y_intercept=None,
                       plot_relative=None, bar_width=0.5, bar_offset=0.25,
                       title=None, save_fig=False, fig_name=None, **kwargs):
    """
    plot groups of compartments i.e. multiple compartments, such as cascade. 
    """
    if colormappings is not None:
        logging.info("colormappings is ignored for plotCompareCascade")
    if colors is None or len(colors) < len(resultset):
        logging.info("colors updated")
        colors = gridColorMap(len(resultset))

    # calculate required barwidth
    barwidth_r = bar_width / (len(resultset))
    bar_offset_r = bar_offset
    xinds = np.arange(len(output_labels)) + bar_offset_r + (barwidth_r * (len(resultset)) / 2.)
    xlim = (0, len(output_labels))
    # for each result, plot cascade passing fig back along
    fig, _ = pl.subplots()
    for i, rname in enumerate(resultset.keys()):
        colorgroup = 10 * [colors[i]]

        fig = innerPlotBar(proj, [resultset[rname]], year_periods=year_periods, output_labels=output_labels,
                       compare_type=COMPARETYPE_CASCADE, pop_labels=pop_labels, ylabel="Number of cases",
                       colors=colorgroup, linestyles=linestyles,
                       plot_relative=plot_relative, y_intercept=y_intercept, xlim=xlim,
                       fig=fig, barwidth=barwidth_r, bar_offset=bar_offset_r + i * barwidth_r, xinds=xinds,
                       title=title, save_fig=save_fig, fig_name=fig_name, **kwargs)
    return fig

def plotPopulationCrossSection(proj, results, output_labels=None, pop_labels=None,
               plot_total=False, plot_type=None,
               plot_observed_data=True, observed_data_label="alive",
               colormappings=None, colors=None, linestyles=None, cat_labels=None,
               title=None, save_fig=False, fig_name=None, **kwargs):
    """

    """
    sim_settings = results.sim_settings
    # setup: determine compartment indices to be plotted
    if output_labels is None:
        output_labels = sorted(results.m_pops[0].comp_ids, key=results.m_pops[0].comp_ids.get)

    # select only compartments that are plottable
    output_labels = [comp_label for comp_label in output_labels if isPlottableComp(comp_label, sim_settings, results.comp_specs)]

    if plot_total:
        fig = innerPlotTrend(proj, [results], output_labels=output_labels,
                   compare_type=COMPARETYPE_VALUE,
                   pop_labels=pop_labels, plot_total=plot_total, plot_type='stacked',
                   plot_observed_data=plot_observed_data, observed_data_label=observed_data_label,
                   colormappings=colormappings, colors=colors, cat_labels=cat_labels,
                   save_fig=save_fig, fig_name=fig_name, **kwargs)
    else:
        if pop_labels is None:
            pop_labels = getPops(results)
            # plot for each population
        for pop in pop_labels:
            fig = innerPlotTrend(proj, [results], output_labels=output_labels,
                   compare_type=COMPARETYPE_VALUE,
                   pop_labels=[pop], plot_total=plot_total, plot_type='stacked',
                   plot_observed_data=True, observed_data_label=observed_data_label,
                   colormappings=colormappings, colors=colors, cat_labels=cat_labels,
                   save_fig=save_fig, fig_name="%s_%s" % (fig_name, pop), **kwargs)
        logger.info("Created multiple plots for plotCompareResults for multiple populations. Returning last plot created")

    return True

def plotPopulationCrossSectionBar(proj, results, output_labels, pop_labels=None, year_periods=None,
                       plot_total=False, plot_observed_data=True, observed_data_label=None,
                       plot_type=None,
                       colormappings=None, colors=None, linestyles=None, y_intercept=None,
                       title=None, save_fig=False, fig_name=None, **kwargs):
    """
    
    """
    if plot_type is None:
        plot_type = PLOTTYPE_BARSTACKED

    for out_label in output_labels:
        fig = innerPlotBar(proj, [results], year_periods=year_periods, output_labels=[out_label],
                           compare_type=COMPARETYPE_POP, plot_total=plot_total,
                           pop_labels=pop_labels, plot_type=plot_type,
                           colormappings=colormappings, colors=colors, linestyles=linestyles, plot_relative=None,
                           title=title, save_fig=save_fig, fig_name=fig_name + "_%s" % out_label, **kwargs)
    return fig


def innerPlotTrend(proj, resultset, output_labels, pop_labels=None,
                   compare_type=None, plot_total=False, plot_type=None,
                   plot_observed_data=True, observed_data_label=None,
                   plot_relative=None, plot_ybounds=None,
                   colormappings=None, colors=None, linestyles=None, cat_labels=None,
                   title=None, save_fig=False, fig_name=None, **kwargs):
    """
    Common functionality, used by plotResult and plotCompareResults
    
    Plots either characteristics, compartment size, or flow rate, as line, across a set of results for a set of populations. 
    
    Supports plotting of:
        - multiple populations
        - plotting total across multiple populations
        - plotting observed datapoints (where applicable)
    
    If neither colormappings or colors is specified, default color list is used from gridColorMap()
    
    Params:
        proj            project object, containing plotting settings and observed data points
        result          result object 
        output_labels   list of compartment labels, flow rate labels, characteristics
        compare_type    compare type. Currently supported = result | pop | value | year
        pop_labels      populations to be plotted. Default (None) plots all populations.
        plot_total      plot total across populations
        plot_type       plot type. Currently supported = line | stacked   (see PLOTTYPE_ fields)
                        if None, default value = line
        plot_observed_data    add observed datapoints as scatter plot, if corresponding datapoints exist. See note.
        observed_data_label_dict   dict of mappings of what 
        colormappings   colormappings that should be used to generate colors for populations. Supercedes colors. 
                        Format: odict with color / colormap as key, value = list of populations for corresponding key 
        colors          list of colors that should be used for population. Superceded by colormappings.
        linestyles      odict of (k,v): (linestypes, list of populations). Uses default linestyle if not supplied.
        title           title for plot
        save_fig        boolean flag, whether to save plot
        fig_name        if plot is saved, filename

        
    """
    # -------------------------------------------------------
    # extract relevant objects
    data = proj.data
    settings = proj.settings
    plotdict = proj.settings.plot_settings
    charac_specs = proj.settings.charac_specs
    plot_over = (proj.settings.tvec_start, proj.settings.tvec_end)
    if 'xlim' in kwargs.keys() and kwargs['xlim'] is not None:
        plot_over = kwargs['xlim']
    try:
        tmp_plotdict = dcp(plotdict)
        tmp_kwargs = dcp(kwargs)
    except NotImplementedError:
        tmp_plotdict = ndcp(plotdict)
        tmp_kwargs = ndcp(kwargs)

    tmp_plotdict.update(kwargs)

    # -------------------------------------------------------
    # generic setup for data
    if pop_labels is None:
        pop_labels = getPops(resultset[0])

    if fig_name is None:
        fig_name = "PlotValue"

    if plot_type is None:
        plot_type = PLOTTYPE_LINE

    logging.debug("compare_type = " + compare_type)
    if compare_type == COMPARETYPE_RESULT:
        series_labels = resultset.keys()
    elif compare_type == COMPARETYPE_VALUE:
        series_labels = output_labels
    elif compare_type == COMPARETYPE_POP:
        if plot_total:
            series_labels = ['Total']
        else:
            series_labels = pop_labels
    else:
        logger.info("Plotting: compare_type not specified, assuming comparing populations")
        raise OptimaException("Unknown compare_type for plotting: %s" % compare_type)
#         series_labels = pop_labels

    if observed_data_label is None:
        observed_data_label = output_labels[0]

    if title is None:
        title = ""

    if len(output_labels) == 1:
        name = output_labels[0]
    else:
        name = "Combined_%s_%s" % (output_labels[0], output_labels[-1])

    # -------------------------------------------------------
    # generic setup for colors, line and hatches
    colors, linestyles, hatches, cat_colors = setupStylings(colormappings, colors, linestyles, series_labels, plotdict)
    if cat_labels is not None:
        legend_labels = cat_labels
        legend_cols = cat_colors
    else:
        legend_labels = series_labels
        legend_cols = None # technically, it is colors, but if legend_cols => None, then the legend plots to

    # -------------------------------------------------------
    # loop over values to be plotted
    # reset reused variables
    ys = []
    ts = []
    dataobs = ([], [])
    unit_tag = ""

    for value_label in output_labels:

        # get values
        if compare_type == COMPARETYPE_RESULT:
            # we loop over results for a single population set
            for resultname in resultset.keys():
                result = resultset[resultname]
                y, t = result.getValuesAt(value_label, year_init=plot_over[0], year_end=plot_over[1], pop_labels=pop_labels, integrated=False)
                y, unit_tag = _convertPercentage(y, value_label, charac_specs)
                ys.append(y)
                ts.append(t)

        elif compare_type == COMPARETYPE_VALUE:
            result = resultset[0]
            y, t = result.getValuesAt(value_label, year_init=plot_over[0], year_end=plot_over[1], pop_labels=pop_labels, integrated=False)
            y, unit_tag = _convertPercentage(y, value_label, charac_specs)
            ys.append(y)
            ts.append(t)

        elif compare_type == COMPARETYPE_POP:
            result = resultset[0]
            if plot_total:

                y, t = result.getValuesAt(value_label, year_init=plot_over[0], year_end=plot_over[1], pop_labels=pop_labels, integrated=False)
                y, unit_tag = _convertPercentage(y, value_label, charac_specs)
                ys = [y]
                ts = [t]
            else:
                for pop in pop_labels:
                    y, t = result.getValuesAt(value_label, year_init=plot_over[0], year_end=plot_over[1], pop_labels=pop, integrated=False)
                    y, unit_tag = _convertPercentage(y, value_label, charac_specs)
                    ys.append(y)
                    ts.append(t)


    # get observed data points
    if plot_observed_data:
        dataobs = _extractDatapoint(result, proj, observed_data_label, pop_labels, charac_specs, plot_total=plot_total)
        dataobs, unit_tag = _convertPercentage(dataobs, value_label, charac_specs)

    fullname = getName(observed_data_label, proj)
    if plotdict.has_key('use_full_labels') and plotdict['use_full_labels']:
        name = fullname

    # if plotting relative, get relative values for 100%
    if plot_relative is not None:
        plot_relative_values, relative_tag = _calcRelativeDatapoint(plot_relative, ys, ts)
        ys = np.array(ys)
        ys = ys / plot_relative_values[:, None]
        ys *= 100
        fullname = fullname + relative_tag
        unit_tag = ' (%)' # TODO get rid of magic variable

    if plot_ybounds is not None:
        # it should be a tuple with either (label_low, label_high) or [(data_years, data_low, data_high), (]
        if len(plot_ybounds) == 3:
            tmp_plotdict['y_bounds'] = plot_ybounds
        elif len(plot_ybounds) == 2:
            dataobs_yrs, dataobs_low = _extractDatapoint(result, proj, plot_ybounds[0], pop_labels, charac_specs, plot_total=plot_total)
            _, dataobs_high = _extractDatapoint(result, proj, plot_ybounds[1], pop_labels, charac_specs, plot_total=plot_total)
            # Cleaning data structures ...
            dataobs_low = dataobs_low
            dataobs_high = dataobs_high
            # and removing nan values
            dataobs_yrs = dataobs_yrs[~np.isnan(dataobs_low)]
            dataobs_low = dataobs_low[~np.isnan(dataobs_low)]
            dataobs_high = dataobs_high[~np.isnan(dataobs_high)]
            tmp_plotdict['y_bounds'] = (dataobs_yrs, dataobs_low, dataobs_high)
        else:
            raise OptimaException("Unknown format for ybounds:" , plot_ybounds)


    # setup for plot:
    final_dict = {
              'xlabel':'Year',
              'ylabel': fullname + unit_tag,
              'title': '%s' % title,
              'save_figname': '%s_%s' % (fig_name, name),
              'y_hat': dataobs[1],
              't_hat': dataobs[0]}

    tmp_plotdict.update(final_dict)
    # plot values
    fig = _plotTrends(ys, ts, legend_labels, plot_type=plot_type,
            save_fig=save_fig, colors=colors, cat_colors=legend_cols,
            linestyles=linestyles, hatches=hatches, **tmp_plotdict)

    # plot separate legend
    if plotdict.has_key('legend_off') and plotdict['legend_off']:
        # Note that color list may be different to colors, as it represents
        # classes of compartments i.e. ['Latent disease states','Active disease states']
        legendsettings = plotdict['legendsettings']
        separateLegend(labels=legend_labels, colors=cat_colors, fig_name=fig_name, linestyles=linestyles, **legendsettings)

    return fig



def innerPlotBar(proj, resultset, output_labels, pop_labels=None,
                 compare_type=None, plot_type=None, year_periods=None,
                 plot_total=False, plot_relative=None,
                 observed_data_label=None, ylabel=None, xlabels=None,
                 colormappings=None, colors=None, linestyles=None, cat_labels=None,
                 title=None, save_fig=False, fig_name=None, **kwargs):
    """
    Params:
        proj
        results
        years
        output_labels    list, or odict of lists, of 
        compare_type
        pop_labels
        colormappings
        colors
        linestyles
        cat_labels
        title
        save_fig
        fig_name
        
    Types of plots:
        compare_type = 'results' : each bar will be the scenario result, representing 
                                    a value for a given population subset for multiple year objects 
                                i.e. scenario results, for the number diagnosed for all high risk pops
                                If plot_total = True: plot as stacked
                                If plot_total = False: plot with bargroups
        compare_type = 'pops': each bar will be the population result, representing a 
                                    value for a multiple year objects, with the same result
                                i.e. for specified populations, what is the number of new infections across 
                                2015 vs 2025. 
                                If plot_total = True: plot as stacked
                                If plot_total = False: plot with bargroups
        compare_type = 'value': each bar is the result for a group of values
                                 plot_total --> False        
        compare_type = 'years': 
                                plot_total --> False
                                
                                
    Examples:
        # Comparison of input to compartment and output from compartment for all populations:
        TODO

        # Comparison of input to compartment and output from compartment for a subset of populations:
        TODO

        # Comparison of number of active TB cases for populations for total pop, for 2017 vs 2025
        TODO
                
        # Comparison of number of active TB cases for populations per population, for 2017 vs 2025
        TODO
        
    """
 # -------------------------------------------------------
    # extract relevant objects
    data = proj.data
    settings = proj.settings
    plotdict = proj.settings.plot_settings
    charac_specs = proj.settings.charac_specs
    plot_over = (proj.settings.tvec_start, proj.settings.tvec_observed_end)
    try:
        tmp_plotdict = dcp(plotdict)
        tmp_kwargs = dcp(kwargs)
    except NotImplementedError:
        tmp_plotdict = ndcp(plotdict)
        tmp_kwargs = ndcp(kwargs)
    tmp_plotdict.update(kwargs)
    # -------------------------------------------------------
    # generic setup for data
    unit_tag = ""

#     if plot_type is None:
#         plot_type = PLOTTYPE_BARSTACKED
#     if plot_type == PLOTTYPE_BAR:
#         plot_stacked = False
#     elif plot_type == PLOTTYPE_BARSTACKED:
#         plot_stacked = True
#     else:
#         raise OptimaException("Plotting: only plot_type=%s or %s allowed; plot_type=%s supplied" % (PLOTTYPE_BAR, PLOTTYPE_BARSTACKED, plot_type))
    plot_stacked = True # TODO implement non-stacked bar plots, update this and uncomment the above code


    if isinstance(output_labels, odict):
        if not isinstance(output_labels[0], list):
            raise OptimaException("Plotting: output_labels structure must have values provided as a list ")
        if not (compare_type == COMPARETYPE_VALUE or compare_type == COMPARETYPE_CASCADE):
            raise OptimaException("Plotting: odict only valid for comparing values. Current attempt for compare_type=%s" % compare_type)


    # should be a list with at least one year
    if year_periods is None:
        year_periods = [plotdict['default_year']]
    elif isinstance(year_periods, numbers.Real):
        year_periods = [year_periods]
    elif isinstance(year_periods, list):
        if len(year_periods) > 2 and (compare_type == COMPARETYPE_RESULT or compare_type == COMPARETYPE_VALUE or compare_type == COMPARETYPE_CASCADE):
#             logging.warn("Multiple years chosen; unable to determine which year set to be used from year=%s" % (str(year_periods)))
            logging.info("Note that will be comparing over period = (%g, %g)" % (year_periods[0], year_periods[-1]))
    else:
        raise OptimaException("Unknown time format for year_periods:" + year_periods)

    if ylabel is None: # string on y-axis
        if observed_data_label is not None:
            ylabel = getName(observed_data_label, proj) + unit_tag
        elif len(output_labels) == 1:
            ylabel = getName(output_labels[0], proj) + unit_tag
        else:
            ylabel = getName(None, proj) + unit_tag

    # components of the overall title
    if pop_labels is None: # detail what populations
        pop = plotdict['default_pops']
    elif pop_labels is dict:
        pop = ' vs '.join(pop_labels.keys())
    else:
        pop = ' vs '.join(pop_labels)
    years = "%g" % year_periods[0] # detail what years
    if len(year_periods) >= 2:
        years += " to %g" % year_periods[-1]
    # update title to include population and year details
    if title is None:
        title = ylabel
        title = "%s\n%s, %s" % (title, pop, years)

    if pop_labels is None: # explicitly update populations
        pop_labels = getPops(resultset[0])

    if fig_name is None:
        fig_name = plotdict['default_figname']

    # set series (x-axis labels) and color labels (colormapping/legend elements),
    # dependent on plot type
    if compare_type == COMPARETYPE_RESULT:
        series_labels = resultset.keys()
        if plot_total:
            color_labels = ['Total']
        else:
            color_labels = pop_labels
    elif compare_type == COMPARETYPE_VALUE:
        if isinstance(output_labels, odict):
            color_labels = []
            series_labels = output_labels.keys()
            for series in series_labels:
                color_labels += output_labels[series]
        else:
            series_labels = output_labels
            color_labels = output_labels
    elif compare_type == COMPARETYPE_CASCADE:
        if isinstance(output_labels, odict):
            series_labels = output_labels.keys()
        else:
            series_labels = output_labels

        if isinstance(pop_labels, dict):
            color_labels = pop_labels.keys()
        else:
            color_labels = pop_labels

    elif compare_type == COMPARETYPE_YEAR:
        series_labels = getYearLabels(year_periods)
        color_labels = output_labels
    elif compare_type == COMPARETYPE_POP:
        if plot_total or plot_stacked :
            series_labels = getYearLabels(year_periods)
        else:
            series_labels = pop_labels
        color_labels = pop_labels
    else:
        logger.info("Plotting: compare_type not specified, assuming comparing populations")
        series_labels = pop_labels
        color_labels = pop_labels

    if xlabels is not None:
        series_labels = xlabels

    # category labels
    if cat_labels is not None:
        legend_labels = cat_labels # TODO revisit
    elif plotdict.has_key('use_full_labels') and plotdict['use_full_labels']:
        legend_labels = [getName(lab, proj) for lab in color_labels]
    else:
        legend_labels = color_labels


    # -------------------------------------------------------
    # generic setup for colors, line and hatches
    colors, linestyles, hatches, cat_colors = setupStylings(colormappings, colors, linestyles, color_labels, plotdict)

    # -------------------------------------------------------
    # extract values from results for plotting
    values = []
    if compare_type == COMPARETYPE_RESULT:
        value_label = output_labels[0]
        time_period = year_periods
        for resultname in resultset.keys():
            result = resultset[resultname]
            ys = []
            if plot_total:
                y = getValueHandler(proj, result, value_label, time_period, pop_labels)
                ys.append(y)
            else:
                for pop in pop_labels:
                    y = getValueHandler(proj, result, value_label, time_period, [pop])
                    ys.append(y)
            values.append(ys)

    elif compare_type == COMPARETYPE_VALUE:
        time_period = year_periods # TODO confirm that only one time period required
        result = resultset[0] # TODO confirm that only one time period required
        for output_key in output_labels.keys():
            output_group = output_labels[output_key]
            ys = []
            for value_label in color_labels:
                if value_label in output_group:
                    y = getValueHandler(proj, result, value_label, time_period, pop_labels)
                else:
                    y = 0
                ys.append(y)
            values.append(ys)

    elif compare_type == COMPARETYPE_CASCADE:
        # similar to COMPARETYPE_VALUE, but we aggregate over populations
        time_period = year_periods # TODO confirm that only one time period required
        result = resultset[0] # TODO confirm that only one time period required
        for output_key in output_labels.keys():
            output_group = output_labels[output_key]
            ys = []
            if isinstance(pop_labels, dict):
                for pop_lab, pop_set in pop_labels.iteritems():
                    count = 0
                    for lab in output_group:
                        y = getValueHandler(proj, result, lab, time_period, pop_set)
                        count += y
                    ys.append(count)
            else:
                for pop in pop_labels:
                    count = 0
                    for lab in output_group:
                        y = getValueHandler(proj, result, lab, time_period, [pop])
                        count += y
                    ys.append(count)
            values.append(ys)

    elif compare_type == COMPARETYPE_YEAR:
        # loop over values, over year periods
        result = resultset[0]
        for time_period in year_periods:
            ys = []
            for value_label in output_labels: # TODO extend so that it returns per pop_labels
                y, _ = result.getValuesAt(value_label, year_init=time_period, pop_labels=pop_labels, integrated=False)
                y = y[0]
                y, unit_tag = _convertPercentage(y, value_label, charac_specs)
                ys.append(y)
            values.append(ys)

    elif compare_type == COMPARETYPE_POP:
        value_label = output_labels[0] # TODO check that output_labels len == 1
        result = resultset[0]
        for time_period in year_periods:
            ys = []
            for pop in pop_labels:
                y, _ = result.getValuesAt(value_label, year_init=time_period, pop_labels=[pop], integrated=False)
                y = y[0]
                y, unit_tag = _convertPercentage(y, value_label, charac_specs)
                ys.append(y)
            values.append(ys)
    else:
        logger.error("Non valid type encountered")
        raise OptimaException("Unknown compare_type: %s" % compare_type)

    if plot_relative is not None:
        plot_relative_values, relative_tag = _calcRelativeDatapoint(plot_relative, values)
        values = np.array(values)
        values = values / plot_relative_values[:, None]
        values *= 100
        ylabel = ylabel + relative_tag # ", relative to <AMOUNT> (%)" # TODO get rid of magic variable


    # update (and selective update) the plotting dict
    final_dict = {'ylabel': ylabel,
                  'plot_stacked': plot_stacked,
                  'save_figname': fig_name}
    print "\n"*5
    print tmp_plotdict
    print "\n"*5

    tmp_plotdict.update(final_dict)
    # separately update the title, to allow if the plot settings would allow for plots
    if tmp_plotdict.has_key('title') and tmp_plotdict['title'] is None:
        tmp_plotdict['title'] = title


    fig = _plotBars(values, labels=legend_labels, colors=colors,
            linestyles=linestyles, hatches=hatches, xlabels=series_labels,
            save_fig=save_fig, **tmp_plotdict)

    if plotdict.has_key('legend_off') and plotdict['legend_off']:
        # Note that color list may be different to colors, as it represents
        # classes of compartments i.e. ['Latent disease states','Active disease states']
        legendsettings = plotdict['legendsettings']
        # TODO: fix usage when legend should use colors rather than cat_colors
        separateLegend(labels=legend_labels, colors=cat_colors, fig_name=fig_name, linestyles=linestyles, **legendsettings)

    return fig


def getValueHandler(proj, result, value_label, year_period, pop_labels):
    """
    Simple wrapper around getValues for usage within innerPlotBar
    """
    charac_specs = proj.settings.charac_specs
    if len(year_period) >= 2:
        y, _ = result.getValuesAt(value_label, year_init=year_period[0], year_end=year_period[-1], pop_labels=pop_labels, integrated=True)
    else:
        y, _ = result.getValuesAt(value_label, year_init=year_period[0], pop_labels=pop_labels, integrated=False)
        y = y[0]
    y, unit_tag = _convertPercentage(y, value_label, charac_specs)
    return y


def plotBudgets(budgets, settings, title="", labels=None, xlabels=None, xlabel="Budget", currency="USD",
                colormappings=None, cat_labels=None, use_full_labels=False, full_labels=None,
                save_fig=False, fig_name=None, legendsettings=None, linestyles=None):
    """
    
    Params:
        budgets     list of dicts, with key:val of program:budget
        title       string, with plot title
        labels      list of programs
    """
    xlim = 3
    if len(xlabels) > 3:
        xlim = len(xlabels)

    plotdict = settings.plot_settings
    if plotdict is None:
        plotdict = {}

    # setup: determine colors to be used
    if labels is None:
        # create super set of all programs. We could use itertools, but we'll use maps
        progkeys = [b.keys() for b in budgets]
        labels = list(set.union(*map(set, progkeys)))
        labels.sort()

    colors = []
    print linestyles
    colors, linestyles, hatches, cat_colors = setupStylings(colormappings, colors, linestyles, labels, plotdict)

#     if colormappings is not None:
#         colors_dict, cat_colors = getCategoryColors(colormappings, 'sequential')
#         # reorder so that colors are same as expected for plotting the population
#         for (j, prog_label) in enumerate(labels):
#             colors.append(colors_dict[prog_label])
#
#     if linestyles is None:
#         linestyles = getLinemapping(linestyle_dict)

    if legendsettings is None:
        legendsettings = {}


    # unfortunately we have to do it this way to ensure that the programs are all extracted in the same order
    values = [[b[k] if b.has_key(k) else 0 for k in labels ] for b in budgets]


    final_dict = {'xlim': (0, xlim),
                  'title': "", # 'Budgets for %s' % (title),
                  'ylabel': "%s (%s)" % (xlabel, currency),
                  'save_figname': '%s_budget' % fig_name}
    plotdict.update(final_dict)

    if use_full_labels:
        clabels = full_labels
    else:
        clabels = labels

    _plotBars(values, clabels, colors=colors, linestyles=linestyles, hatches=hatches,
              xlabels=xlabels, # legendsettings=legendsettings,
              save_fig=save_fig, **plotdict)

    if plotdict.has_key('legend_off') and plotdict['legend_off']:
        # Do this separately to main iteration so that previous figure are not corrupted
        # Note that colorlist may be different to colors, as it can represent
        # classes of budgets
        # reverse legend order so that it matches top<->bottom of stacked bars
        if use_full_labels:
            # legendsettings = {'ncol':2}
            separateLegend(labels=full_labels, colors=colors, linestyles=linestyles,
                           fig_name=fig_name, reverse_order=True, **legendsettings)
        else:
            separateLegend(labels=cat_labels, colors=cat_colors, linestyles=linestyles,
                           fig_name=fig_name, reverse_order=True,)

def getName(output_id, proj):
    """
    For a given output_id, returns the user-friendly version of the name. 
    
    This works for:
    - characteristic labels
    - compartment labels
    - parameter labels
    - parameter tags
    
    Specifying output_id = None will return the default ylabel as specified in settings.plot_settings dict 
    """
    settings = proj.settings
    if isinstance(output_id, list):
        name = "List_(%s,%s)" % (output_id[0], output_id[-1])
    elif output_id in settings.charac_specs: # characteristic
        name = settings.charac_specs[output_id]['name']
    elif output_id in settings.linkpar_specs: # parameter
        name = settings.linkpar_specs[output_id]['name']
    elif output_id in settings.node_specs: # compartment
        name = settings.node_specs[output_id]['name']
    elif output_id in proj.data['pops']['label_names'].keys(): # population label
        name = output_id
    elif output_id is None: # default label, as user hasn't specified ylabel for aggregate value
        name = settings.plot_settings["default_ylabel"]
    else:
        try:
            # this is a parameter specified by tag, rather than by parameter label
            # TODO: propose that this could more effective by moving the tag->label dict into settings, and generated upon init
            tags = {settings.linkpar_specs[label]['tag']: label for label in settings.linkpar_specs.keys() if settings.linkpar_specs[label].has_key('tag')}
            tmp_id = tags[output_id]
            name = settings.linkpar_specs[tmp_id]['name'] + " (%s)" % settings.plot_settings['effective_rate']
        except:
            name = "Unknown_%s" % output_id
            logging.warn('ERROR: Attempting to plot characteristic "%s" but cannot locate it in either characteristic or parameter specs.' % output_id)
    return name

def getPops(result):
    """
    Returns the full list of populations
    """
    return [pop.label for pop in result.m_pops]

def _calcRelativeDatapoint(plot_relative, ys, ts=None):
    """
    Calculates datapoints corresponding to 100%, so that trends and values can then 
    be normalized. The structure is flexible to allow for normalized either based on 
    an explicit value, a year, or for intended usage for bar plots, to be normalized to 100%
    for each bar, or for the first bar only.
    
    Params:
        plot_relative     tuple pair (key, value), where key can be:
                                year
                                value
                                normalized*
                                normalizedFirst*    will normalize to the first height
                          * value is ignored for these keys
        ys
        ts                timepoints. Only required when plotting relative for year, otherwise can be set to None
        
    Returns:
        vals             an 1D numpy array, of length matching the first dimension of ys     
    
    Examples:
        plot_rel_year = ('year', 2015.) # will return a list where each element is the value of 
                                        # each ys at the timepoint for t=2015. 
        plot_rel_val  = ('value', 0.5)  # will return a list where each element is 0.5 
        plot_rel_norm = ('normalized')  # will return a list where each element is the corresponding total for each bar
        plot_rel_normF = ('normalizedYear') # will return a list of len(ys) where each element is the sum of the first bar 
        
    """
    relative_tag = ", normalized"
    if plot_relative[0] == 'year':
        idx = np.nonzero(np.in1d(ts[0], [plot_relative[1]]))[0]
        vals = [ float(yy[idx]) for yy in ys]
        relative_tag = ", relative to %g" % plot_relative[1]
    elif plot_relative[0] == 'value':
        values = plot_relative[1]
        if type(values) == 'list':
            if len(values) >= len(ys):
                vals = values[:len(ys)]
            else:
                vals = values * len(ys)
        else: # else, a number # TODO check and confirm
            vals = [values] * len(ys)
        relative_tag += " to y=%g" % plot_relative[1]
    elif plot_relative == 'normalized' or plot_relative[0] == 'normalized':
        ys = np.array(ys)
        vals = ys.sum(axis=1)
    elif plot_relative == 'normalizedFirst' or plot_relative[0] == 'normalizedFirst' :
        ys = np.array(ys)
        vs = ys.sum(axis=1)
        vals = [vs[0]] * len(vs)
    else:
        raise OptimaException("Unknown indicator for plot_relative: %s" % plot_relative[0])

    return np.array(vals), relative_tag


def _extractDatapoint(results, proj, value_label, pop_labels, charac_specs, plot_total=False):
    """
    Extract the characteristic datapoints for populations
    
    Param:
        data     proj.data
        value_label     string for characteristic label
        pop_labels    list of populations
        charac_specs    
        plot_total    sums over observed datapoints
        
    Returns:
        tuple with (timepoints, y-values) for each population
        
    TODO: make proj parameters, so that charac_specs and data are encapsulated
    """
    dataobs = None
    data = proj.data

    data_locations = ['characs', 'linkpars']
    for data_loc in data_locations:
        if value_label in data[data_loc].keys():
            break # we've likely identifed where it is

    try:
        ys = [data[data_loc][value_label][poplabel]['y'] for poplabel in pop_labels]
        ts = [data[data_loc][value_label][poplabel]['t'] for poplabel in pop_labels]
    except:
        """
        ys = [list(proj.parsets[0].getPar(value_label).interpolate(tvec=t, pop_label=poplabel, extrapolate_nan=True)) for poplabel in pop_labels]
        ts = [t for poplabel in pop_labels]
        """
        logging.info("Could not find datapoints with label '%s'" % value_label)
        ys = [[]]
        ts = [[]]

    try:
        if 'plot_percentage' in charac_specs[value_label].keys():
            ys *= 100
    except:
        pass

    if plot_total:
        # we have to be a bit more careful, as there could be missing elements for the different timesteps
        try:
            ys = np.array(ys)
            ys = ys.sum(axis=0)
            ts = ts[0]
        except ValueError:
            # could raise a ValueError if ys elements are different lengths
            ys, ts = [], []

    dataobs = (ts, ys)
    print dataobs
    return dataobs

def _convertPercentage(datapoints, output_label, charac_specs):
    """
    Checks whether output characteristic should be plottable as a 
    percentage, and if so, converts data points to percentage i.e. 0.02 to 2 (%). 
    
    Params:
        datapoints
        output_label
        charac_specs
        
    Returns:
        datapoint
        unittags
    """
    if output_label in charac_specs and 'plot_percentage' in charac_specs[output_label].keys():
        datapoints *= 100
        unit_tags = ' (%)' # TODO remove magic string
    else:
        unit_tags = ''

    return datapoints, unit_tags


def isPlottableComp(comp_label, sim_settings, comp_specs):
    """
    Returns bool indicating whether a population label should be included in metrics
    for population reporting when plotting cascade
    """
    # TODO move this method to sim_settings / comp_specs superobject
    if comp_label not in comp_specs.keys():
        logging.info("Attempting to plot a non-compartment: %s" % comp_label)
        return True # TODO revisit this assumption
    if comp_label in sim_settings['tag_no_plot']:
        return False
    if comp_specs[comp_label].has_key('junction'):
        return False
    if comp_specs[comp_label].has_key('tag_birth'):
        return False
    if comp_specs[comp_label].has_key('tag_dead'):
        return False
    return True

def _plotYIntercept(ax, y_intercept, xlim, **kwargs):
    yintercept_plotdict = {}
    if kwargs.has_key('y_intercept_line'):
        yintercept_plotdict = kwargs['y_intercept_line']
    ax.hlines([y_intercept], xmin=xlim[0], xmax=xlim[1], zorder=1, **yintercept_plotdict)



def _plotTrends(ys, ts, labels, colors=None, y_hat=[], t_hat=[], plot_type=None, cat_colors=None,
             legendsettings=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, y_ticks=None, x_ticks=None,
             y_intercept=None, reverse_order=False, y_bounds=None, linestyles=None, hatches=None,
             smooth=False, symmetric=False, repeats=5,
             alpha=0.3, marker='o', s=40, facecolors='none', linewidth=3, zorder=10,
             save_fig=False, save_figname=None, legend_off=False,
             box_width=0.9, box_offset=0.0, formatter=None, **kwargs):
    """
    Plots multiple lines, with additional option of overlaying observed datapoints
    
    Params:
        ys        list of values for ys, with each entry corresponding to a line
        ts        list of values for xs, with each entry corresponding to a line
        labels    list of labels for each line
        colors    list of colors for each line
        y_hat
        t_hat
        plot_type
        cat_colors
        legendsettings
        title
        xlabel
        ylabel
        xlim
        ylim
        x_ticks
        y_ticks
        y_intercept
        reverse_order
        y_bounds    list of array for each ys entry, with format of (tbound, ybound_min, ybound_ymax), thus can be specified independently of ts
        linestyles
        hatches
        smooth
        symmetric
        repeats
        alpha
        marker
        s
        facecolors
        linewidth
        zorder
        save_fig
        save_figname
        legend_off
        box_width
        box_offset
        formatter
        **kwargs    further keyword arguments, such as ylims, legend_off, edgecolors, etc.
        
    Returns
        fig         figure handle
        
    Replaces
         _plotLine, _plotStackedCompartments
         
        
    """
    logging.debug("Xlim = ", xlim)
    if len(ys) == 0:
        logging.error("No values supplied; cannot plot")
        raise OptimaException("No values supplied; cannot plot")

    if plot_type is None:
        plot_type = PLOTTYPE_LINE

    # Some pre-processing and setting first pass values for ymin / ymax
    if xlim is None: xlim = (ts[0][0], ts[0][-1])
    ymin_val = np.min(ys[0])
    indices = (ts[0] >= xlim[0]) * (ts[0] <= xlim[1])
    ymax_val = np.max(ys[0][indices])

    # Setup of plotting bling
    if legendsettings is None: legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':1}

    if colors is None or len(colors) < len(ys):
        colors = gridColorMap(len(ys))
        logger.info("Plotting: setting color scheme to be default colormap, as not all lines had color assigned")

    if linestyles is None or len(linestyles) < len(ys):
        try: linestyles = [kwargs['default_linestyle']] * len(ys)
        except: linestyles = ['-'] * len(ys)
        logger.info("Plotting: setting linestyles to be default value, as not all lines had styles assigned")

    if hatches is None or len(hatches) < len(ys):
        try: hatches = [kwargs['default_hatch']] * len(ys)
        except: linestyles = [None] * len(ys)
        logger.info("Plotting: setting hatches to be default value, as not all hatches had styles assigned")

    # Plotting proper
    fig, ax = pl.subplots()
    bottom = 0  # setup if we are using stacked plots

    if y_intercept is not None:
        _plotYIntercept(ax, y_intercept, xlim, **kwargs)

    # plot ys, but reversed - and also reverse the labels (useful for scenarios, and optimizations):
    order_ys = range(len(ys))
    if reverse_order:
        logger.info("Reversing order of plot lines")
        order_ys = order_ys[::-1]  # surely there are more elegant ways to do this ...
        labels = labels[::-1]

    for k in order_ys:

        yval = ys[k]

        # if there are confidence bounds provided, plot using fill_between
        if y_bounds is not None:
            try:
                t_bound, y_min_bound, y_max_bound = y_bounds[0], y_bounds[1], y_bounds[2]
            except:
                t_bound, y_min_bound, y_max_bound = zip(*y_bounds[k])[0] , zip(*y_bounds[k])[1] , zip(*y_bounds[k])[2]
#                 t_bound, y_min_bound, y_max_bound = zip(y_bounds[k])[0] , zip(y_bounds[k])[1] , zip(y_bounds[k])[2]
            ax.fill_between(t_bound, y_min_bound, y_max_bound, facecolor=colors[k], alpha=alpha, linewidth=0.1, edgecolor=colors[k])

        # smooth the value
        if smooth:
            yval = smoothfunc(yval, symmetric, repeats)

        # plot the thing
        if plot_type == PLOTTYPE_LINE:
            ax.plot(ts[k], yval, c=colors[k], ls=linestyles[k])
            # keep track of max y-value seen, as sometimes plotting can rescale
            if np.max(yval[indices]) > ymax_val:
                ymax_val = np.max(yval[indices])
            scatter_color = colors[k]
        elif plot_type == PLOTTYPE_STACKED:
            top = bottom + yval
            lw = 1.
            if hatches[k] is None:
                ec = colors[k]
            else:
                ec = kwargs['hatch_bg']
            ax.fill_between(ts[0], bottom, top, facecolor=colors[k], alpha=1, lw=lw, hatch=hatches[k], edgecolor=ec)  # for some reason, lw=0 leads to no plot if we then use fig.savefig()
            reg, = ax.plot((0, 0), (0, 0), color=colors[k], linewidth=10)  # TODO fix this by using xlims and ylims appropriately
            bottom = dcp(top)
            scatter_color = kwargs['marker_color']
        else:
            logger.error("Unknown plot_type = %s. Aborting plotting" % plot_type)

        # scatter points for observed data points
        if len(y_hat) > 0:
            try:
                if len(y_hat[k]) > 0:  # i.e. we've seen observable data
                    ax.scatter(t_hat[k], y_hat[k], marker=marker, edgecolors=scatter_color, facecolors=facecolors, s=s, zorder=zorder, linewidth=linewidth)
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

    # set position
    box = ax.get_position()
    ax.set_position([box.x0 + box.width * box_offset, box.y0, box.width * box_width, box.height])

    # set title, labels and legends
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if formatter is not None:
        ax.yaxis.set_major_formatter(formatter)
    else:
        ax.get_yaxis().get_major_formatter().set_scientific(False)

    if not legend_off:
        if cat_colors is not None:
            import matplotlib.patches as mpatches
            patches = [  mpatches.Patch(color=color, label=label, ec='white', hatch=hatch) for label, color, hatch in zip(labels, cat_colors, hatches)]
            ax.legend(patches, labels, **legendsettings)
        else:
            ax.legend(labels, **legendsettings)

    # automatic choices for setting ylim bounds
    ax.set_ylim(ymin=0)
    ax.set_ylim(ymax=ax.get_ylim()[1] * 1.05)


    # overwrite with specified choice
    if ylim is not None:
        ax.set_ylim(ylim)
    if xlim is not None:
        ax.set_xlim(xlim)
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

def _plotBars(values, labels=None, colors=None, title="", orientation='v', legendsettings=None, y_intercept=None,
              xlabel="", ylabel="", xlabels=None, yticks=None, barwidth=0.5, bar_offset=0.25, xlim=None, ylim=None,
              linewidth=3, inds=None, xinds=None, alphas=None, hatches=None, plot_stacked=True,
              save_fig=False, save_figname=None, fig=None,
              box_width=0.9, box_offset=0.0, legend_off=False, formatter=None, reverse_order=True, **kwargs):
    """
    Plots bar graphs. Intended for budgets, characteristics, and other model output
    
    Params:
        values    list of lists, containing values
        labels    
        colors
        title
        orientation
        legendsettings
        y_intercept
        xlabel
        ylabel
        xlabels
        yticks
        barwidth
        bar_offset
        xlim
        ylim
        linewidth
        inds
        xinds        indices of xlabels
        alphas
        hatches
        plot_stacked
        save_fig
        save_figname
        box_width
        box_offset
        legend_off
        formatter
        reverse_order
        **kwargs
        
    Returns
        fig    figure handle
    
    TODO:
        implement orientation = 'h'
        implement non-stacked plots
    """
    print kwargs.keys()
    print kwargs
    print barwidth
    print bar_offset
    print xlim, "<----------- xlim"
    # ---------------------
    # setup
    num_bars = len(values)
    num_cats = len(values[0])  # label categories
    if inds is None:
        inds = np.arange(num_bars) + bar_offset

    if xinds is not None:
        pass
    elif plot_stacked:
        xinds = np.arange(num_bars) + bar_offset + barwidth / 2.
    else: # TODO modify if not stacked
        xinds = np.arange(num_cats) + bar_offset + barwidth / 2.

    if alphas is None:
        alphas = [1.] * num_bars * num_cats

    if hatches is None:
        hatches = [None] * num_bars * num_cats

    if xlabels is None:
        x_ticks = (xinds, range(num_bars))
    else:
        x_ticks = (xinds, xlabels)

    if colors is None:
        colors = gridColorMap(num_cats)
        logger.info("Plotting: setting color scheme to be default colormap, as not all lines had color assigned")

    if legendsettings is None: legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':1}

    # preprocessing to make our lives easier:
    cat_values = map(list, zip(*values))
    cumulative = np.zeros(num_bars)

    if not plot_stacked:
        # I hate myself for writing this. I'm so, so, sorry.
        cat_values = values
        num_cats = num_bars

    # ---------------------
    # and plot:
    if fig is None:
        fig, ax = pl.subplots()
    else:
        ax = fig.gca()

    for k in range(num_cats):

        if plot_stacked:
            indices = inds
        else:
            indices = inds[k]

        lw = 0.
        if hatches[k] is None:
            ec = colors[k]
        else:
            ec = kwargs['hatch_bg']

        # as we go through each category, update where we're plotting our bar
        if k == 0: # ... except for the first category, which is in the same location regardless
            ax.bar(indices, cat_values[k], color=colors[k], hatch=hatches[k], edgecolor=ec, width=barwidth, alpha=alphas[k], lw=lw, label=labels[k])
        elif plot_stacked:
            ax.bar(indices, cat_values[k], color=colors[k], hatch=hatches[k], edgecolor=ec, width=barwidth, bottom=cumulative, alpha=alphas[k], lw=lw, label=labels[k])
        else:
            ax.bar(indices, cat_values[k], color=colors[k], hatch=hatches[k], edgecolor=ec, width=barwidth, alpha=alphas[k], lw=lw, label=labels[k])

        if plot_stacked:
            cumulative += cat_values[k]

    # set our ylim and xlim so that if we draw a y_intercept, it will be included
    if ylim is not None:
        ax.set_ylim(ylim)

    if True: # xlim is None:
        xlim = (indices[0] - barwidth / 2, indices[-1] + 0.25 + barwidth)
    ax.set_xlim(xlim)

    if y_intercept is not None:
        _plotYIntercept(ax, y_intercept, xlim, **kwargs)
    # ---------------------
    # post-plotting: set position and formatting
    box = ax.get_position()
    ax.set_position([box.x0 + box.width * box_offset, box.y0, box.width * box_width, box.height])

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if formatter is not None:
        ax.yaxis.set_major_formatter(formatter)
    else:
        ax.get_yaxis().get_major_formatter().set_scientific(False)

    print x_ticks
    print xlabels

    if x_ticks is not None:
        ax.set_xticks(x_ticks[0])
        ax.set_xticklabels(x_ticks[1])

    if yticks is not None:
        ax.set_yticks(yticks[0])
        ax.set_yticklabels(yticks[1])

    if not legend_off:
        handles, labels_leg = ax.get_legend_handles_labels()
        # TODO: this doesn't actually display as reversed - requires fixing
        ax.legend(handles=handles[::-1], labels=labels_leg[::-1], **legendsettings)

    _turnOffBorder()

    if save_fig:
        fig.savefig('%s' % (save_figname))
        logger.info("Saved figure: '%s'" % save_figname)

    return fig

def smoothfunc(ys, symmetric=False, repeats=3):
    """
    Params:
        ys = y values to smooth 
        symmetric = whether or not to use a symmetric convolution kernel
        repeats = the number of times to apply the kernel
        
    Returns:
        smoothed version of ys (as an array)
    
    Notes:
        A symmetric kernel produces less "distortion" for the plot, but
        could exacerbate unphysical effects (e.g., an intervention having
        an impact before it begins). The function pads the ends of the
        time series, so the first and last point shouldn't be affected
        by the smoothing.
        
        The smoothness is directly proportional to the number of repeats.
        In general, 1 repeat will smooth out time point pairs, 2 repeats
        will smooth out time point triplets, etc.
    
    Example:
        import pylab as pl
        y = pl.rand(50)
        ys = smoothfunc(y, symmetric=True, repeats=10)
        pl.plot(y)
        pl.plot(ys)
    """
    ys = np.array(ys)  # Convert to an array

    # Choose the kernel
    if symmetric:
        kernel = np.array([0.25, 0.5, 0.25])  # The tiniest imaginable Gaussian
    else:
        kernel = np.array([0.125, 0.25, 0.5, 0.125])  # The tiniest imaginable asymmetric Gaussian

    npad = repeats * len(kernel)  # Figure out how big the padding on each end needs to be
    yspad = np.concatenate([np.ones(npad) * ys[0], ys, np.ones(npad) * ys[-1]])  # Pad the ends
    for repeat in range(repeats):  # Do the convolution
        yspad = np.convolve(yspad, kernel, 'same')

    ys = yspad[npad:-npad]  # Trim off the padding we added

    return ys

"""
/////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//////////////  All code below this line is slated for removal after code clean up of plotting \\\\\\\\\\\\\\\\
/////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
"""


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
    logging.warn("plotPopulation() is to be deprecated. Please move to new plotting library function")
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
    mpops = results.m_pops
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
    logging.warn("_calculateDatapoints() is to be deprecated. Please move to new plotting library function")

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
    logging.warn("extractCompartment() is to be deprecated. Please move to new plotting library function")

    datapoints, pop_labels, comp_labels = results.getCompartmentSizes(pop_labels=pop_labels, comp_label=comp_labels, use_observed_times=False)
    yhat, that = [], []

    if use_full_labels:
        labels = [comp.label for comp in results.m_pops[0].comps]
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
    logging.warn("extractCharacteristic() is to be deprecated. Please move to new plotting library function")

    datapoints, _, _ = results.getCharacteristicDatapoints(pop_label=pop_labels, char_label=charac_labels, use_observed_times=False)

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
    logging.warn("extractFlows() is to be deprecated. Please move to new plotting library function")

    all_rates = []
    all_tvecs = []

    for (j, pid) in enumerate(pop_labels):

        comp = results.m_pops[pid].getComp(comp_label)

        all_labels = []
        pop_rates = []
        pop_tvecs = []
        for in_out in xrange(2):
            if (in_out == 0 and plot_inflows) or (in_out == 1 and plot_outflows):
                comp_link_ids = [comp.inlink_ids, comp.outlink_ids][in_out]
                label_tag = ['In: ', 'Out: '][in_out]
                for link_tuple in comp_link_ids:
                    link = results.m_pops[link_tuple[0]].links[link_tuple[1]]
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
                            comp_source = results.m_pops[link.index_from[0]].comps[link.index_from[1]]
                        else:
                            comp_source = comp
                        was_proportion = False
                        if link.val_format == 'proportion':
                            denom_val = sum(results.m_pops[lid_tuple[0]].links[lid_tuple[-1]].vals for lid_tuple in comp_source.outlink_ids)
                            num_flow /= denom_val
                            was_proportion = True
                        if link.val_format == 'fraction' or was_proportion is True:
                            if was_proportion is True:
                                num_flow *= comp_source.popsize_old
                            else:
                                num_flow[num_flow > 1.] = 1.
                                num_flow = 1 - (1 - num_flow) ** results.dt  # Fractions must be converted to effective timestep rates.
                                num_flow *= comp_source.popsize
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
    logging.warn("plotSingleCompartmentFlow() is to be deprecated. Please move to new plotting library function")

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
        plot_pids = range(len(results.m_pops))
        plot_pops = [pop.label for pop in results.m_pops]

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

            comp = results.m_pops[pid].getComp(comp_label)

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
    logging.warn("plotPopulationFlows() is to be deprecated. Please move to new plotting library function")

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
        pop_labels = [pop.label for pop in results.m_pops]
        plot_pids = range(len(results.m_pops))

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
    logging.warn("DEPRECATED (plotCharacterstic): please replace usage with plotResult()")
    charac_specs = settings.charac_specs
    # setup:
    if outputIDs is None:
        outputIDs = results.outputs.keys()
        # now only select characteristics which are plottable
        outputIDs = [output_id for output_id in outputIDs if isPlottableCharac(output_id, charac_specs)]
    # else we already know which characteristics to plot
#     print "-------", outputIDs


    if pop_labels is None:
        pop_labels = [pop.label for pop in results.m_pops]
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
    logging.warn("plotScenarioFlows() is to be deprecated. Please move to new plotting library function")

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
    logging.warn("isPlottableCharac() is to be deprecated. Please move to new plotting library function")

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
    logging.warn("getPIDs() is to be deprecated. Please move to new plotting library function")

    pids = []
    for poplabel in poplabels:
        for i, pop in enumerate(results.m_pops):
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
    logging.warn("plotProjectResults() is to be deprecated. Please move to new plotting library function")

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
    logging.warn("plotScenarios() is to be deprecated. Please move to new plotting library function")

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
    logging.warn("_plotLine() is to be deprecated. Please move to new plotting library function")


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
    logging.warn("plotStackedCompartments() is to be deprecated. Please move to new plotting library function")

    if colors is None or len(colors) != len(comps):
        if len(colors) != len(comps):
            logger.info("Plotting: setting color scheme to be default colormap, as not all compartments had color assigned")
        colors = gridColorMap(len(comps))

    # setup
    fig, ax = pl.subplots()
    bottom = 0 * tvec
    max_val = 0

    for (k, comp) in enumerate(comps):
        top = bottom + comp.popsize
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
    logging.warn("plotAllOutflows() is to be deprecated. Please move to new plotting library function")

    mpops = results.m_pops
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
    logging.warn("_plotLine() is to be deprecated. Please move to new plotting library function")

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

    values = [ [scen_results[rname].getValueAt(output, year) for output in output_list] for rname in scen_results.keys()]
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
    logging.warn("plotStackedBarOutputs() is to be deprecated. Please move to new plotting library function")

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
    values = [[results.getValueAt(output_label, year, pop_labels=pop_labels)  for output_label in output_list ]  for year in year_list]

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
    print plotdict
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
    logging.warn("plotValueBars() is to be deprecated. Please move to new plotting library function")

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
    logging.warn("plotCareCascade() is to be deprecated. Please move to new plotting library function")

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


        print values
        print plotdict
        print colors
        print "--> convert"
        new_val = np.array(values)
        new_val.transpose()

#         plotdict['bar_offset'] = i * plotdict['bar_width']
        plotdict['xlim'] = xlim
        plotdict['save_figname'] = fig_name + "_%g" % year_init

        plotdict['plot_stacked'] = True
        plotdict['barwidth'] = 0.8
        print values
        print pop_labels.keys()
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

