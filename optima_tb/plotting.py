#%% Imports
import logging
logger = logging.getLogger(__name__)

from optima_tb.utils import odict

import numpy as np

import pylab as pl
from copy import deepcopy as dcp
import matplotlib.cm as cmx
import matplotlib.colors as colors
from random import shuffle


CMAPS = ['Blues', 'Purples', 'Reds','Oranges',  'Greys', 'Greens', ]
        # Suscep  #Latent  # SP TB    #SN TB    # Dead   # Recovered

 
#%% Function to generate useful colormap for plots

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

    ## Imports
    from numpy import linspace, meshgrid, array, transpose, inf, zeros, argmax, minimum
    from numpy.linalg import norm
    
    # Steal colorbrewer colors for small numbers of colors
    colorbrewercolors = array([
    [27,  158, 119],
    [217, 95,  2],
    [117, 112, 179],
    [231, 41,  138],
    [255, 127, 0],
    [200, 200, 51], # Was too bright yellow
    [166, 86,  40],
    [247, 129, 191],
    [153, 153, 153],
    ])/255.
    
    if ncolors<=len(colorbrewercolors):
        colors = colorbrewercolors[:ncolors]
        
    else: # Too many colors, calculate instead
        ## Calculate sliding limits if none provided
        if limits is None:
            colorrange = 1-1/float(ncolors**0.5)
            limits = [0.5-colorrange/2, 0.5+colorrange/2]
        
        ## Calculate primitives and dot locations
        primitive = linspace(limits[0], limits[1], nsteps) # Define primitive color vector
        x, y, z = meshgrid(primitive, primitive, primitive) # Create grid of all possible points
        dots = transpose(array([x.flatten(), y.flatten(), z.flatten()])) # Flatten into an array of dots
        ndots = nsteps**3 # Calculate the number of dots
        indices = [0] # Initialize the array
        
        ## Calculate the distances
        for pt in range(ncolors-1): # Loop over each point
            totaldistances = inf+zeros(ndots) # Initialize distances
            for ind in indices: # Loop over each existing point
                rgbdistances = dots - dots[ind] # Calculate the distance in RGB space
                totaldistances = minimum(totaldistances, norm(rgbdistances,axis=1)) # Calculate the minimum Euclidean distance
            maxindex = argmax(totaldistances) # Find the point that maximizes the minimum distance
            indices.append(maxindex) # Append this index
        
        colors = dots[indices,:]
    
    ## Wrap up: optionally turn into a list of tuples
    if asarray:
        output = colors
    else:
        output = []
        for i in range(ncolors): output.append(tuple(colors[i,:])) # Gather output
    
    ## For plotting
    if doplot:
        from mpl_toolkits.mplot3d import Axes3D # analysis:ignore
        from pylab import figure, gca
        if newwindow:
            fig = figure()
            ax = fig.add_subplot(111, projection='3d')
        else: 
            ax = gca(projection='3d')
        ax.scatter(colors[:,0], colors[:,1], colors[:,2], c=output, s=200, depthshade=False)
        ax.set_xlabel('R')
        ax.set_ylabel('G')
        ax.set_zlabel('B')
        ax.set_xlim((0,1))
        ax.set_ylim((0,1))
        ax.set_zlim((0,1))
        ax.grid(False)
    
    return output


def _getColormapRange(cmap,ncols=10,order='alternate'):
    """
    Returns a list of colors values, according to colormaps
    
    Params:
        cmap    name of matplotlib.cmap (String)
        ncols   number of colors to be returned (int)
        order   order in which colors should be chosen. Possible values are 'alternate','sequential' or 'random'
        
    Usage:
        cmap = 'jet'
        n = 6
        clist = _getColormapRange(cmap,n,'random')
        for i in range(n):
            pylab.plot(xs[i],ys[i],c=clist[i])
    """
    color_norm  = colors.Normalize(vmin=-1.5, vmax=ncols+0.5)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap=cmap) 
    order_map = range(ncols)
    if order == 'random':
        shuffle(order_map)
    elif order == 'alternate':
        if ncols < 2:
            return _getColormapRange(cmap, ncols, 'sequential')
        a = order_map[ncols/2:]
        b = order_map[:ncols/2]
        order_map[::2] = a 
        order_map[1::2] = b 
    elif order == 'alternate3':
        if ncols < 3:
            return _getColormapRange(cmap, ncols, 'alternate')
        third = ncols/3
        a = order_map[2*third:]
        b = order_map[third:2*third]
        c = order_map[:third]
        order_map[::3] = a 
        order_map[1::3] = b 
        order_map[2::3] = c
    else: # order == 'sequential', so leave in order 
        pass
    
    return [scalar_map.to_rgba(index) for index in order_map]


def getCategoryColors(category_list,order='alternate'):
    """
    For an ordered dictionary of category list, return a mapping of compartment labels to colors (rgba). 
    
    Note that using colormappings will index by compartment keys, so the color maps do not have to be 
    defined in sequence of their appearance within a population's compartment list. 
    The only constraint is that all plottable compartments must be assigned a colour. 
    If a compartment is included in more than one colormap list, its cmap will be set as last list it appeared on. 
    
    Params:
        category_list    odict of list
        order            defines how neighbouring values within a list should be defined. Values are
                            'alternate' 
                            'random'    :
                            'sequential' : 
        
    Returns:
        an odict of (comp_label, rgba) items
        
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
    for k,v in category_list.iteritems():
        
        if isinstance(k, str) and k.startswith('#'):
            # is a hex format 
            tmp_list = [k]*len(v)
        elif isinstance(k, str) : 
            # must be a colormap
            tmp_list = _getColormapRange(k,len(v),order)
        else:
            # else, unknown
            raise OptimaException('Unknown color format: '+k)
        for i,label in enumerate(v):
            col_list[label] = tmp_list[i]
        cat_colors.append(tmp_list[0])
    return col_list,cat_colors
    
    
def separateLegend(labels,colors,fig_name):
    
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    print labels, colors
    fig = plt.figure()
    patches = [
        mpatches.Patch(color=color, label=label)
        for label, color in zip(labels, colors)]
    fig.legend(patches, labels, loc='center', frameon=False)
    plt.savefig("%s_legend"%fig_name)


def _turnOffBorder():
    """
    Turns off top and right borders, leaving only bottom and left borders on.
    """
    pl.gca().spines['right'].set_color('none')
    pl.gca().spines['top'].set_color('none')
    pl.gca().xaxis.set_ticks_position('bottom')
    pl.gca().yaxis.set_ticks_position('left')


def isPlottableComp(comp_label,sim_settings,comp_specs):
    """ 
    Returns bool indicating whether a population label should be included in metrics
    for population reporting when plotting cascade
    """
    # TODO make unplottable comp_specs a setting
    #print comp_label, comp_specs[comp_label].keys()
    if comp_label in sim_settings['tag_no_plot']:
        return False
    if comp_specs[comp_label].has_key('junction'):
        return False
    if comp_specs[comp_label].has_key('tag_birth'):
        return False
    if comp_specs[comp_label].has_key('tag_dead'):
        return False
    return True

def isPlottableCharac(output_id,charac_specs):
    """
    Returns False if specified in cascade spreadsheet as 'n' or 'N' for plot characteristic.
    Else returns True
    """
    try: # Better to try and ask for forgiveness than for permission ... 
        if charac_specs[output_id]['plot_characteristic'].lower() == 'n':
            return False
    except: 
        return True
    
def getPIDs(results,poplabels):
    """
    Takes in a list of poplabels and returns the corresponding PIDs
    
    TODO: this can be improved and made more efficient by either a better implementation of this
    look up, OR (better yet) by improving the data structures of mpops.
    """
    pids = []
    for poplabel in poplabels:
        for i,pop in enumerate(results.m_pops):
            if pop.label == poplabel:
                pids.append(i)
    return pids
    

def plotProjectResults(results,settings, data, title='', colormappings=None, colorlabels=None, pop_labels=None, plot_comp_labels=None, debug=False, plot_observed_data=True, save_fig=False, fig_name=None):
            
    
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
    #plotPopulation(results=results, data=data, title=title, colormappings=colormappings, cat_labels=colorlabels, pop_labels=pop_labels, plot_observed_data=plot_observed_data, \
    #               save_fig=save_fig, fig_name=fig_name, plotdict=plotdict, plot_comp_labels=plot_comp_labels)
    
    # plot characteristics
    plotCharacteristic(results=results, charac_specs=charac_specs, data=data, title=title, plot_observed_data=plot_observed_data, save_fig=save_fig, fig_name=fig_name, plotdict=plotdict)
    # internal plotting
    if debug:
        plotAllOutflows(results)
    
def plotScenarios(scen_results,scen_labels,settings,data,plot_charac=None,plot_pops=None,
                  colormappings=None,colors=None,plot_observed_data=True,save_fig=False,fig_name=None):
    """
    Line plots for scenarios. Should be used for characteristics only.
    
    Params
        scen_results        list of results
        scen_labels         list of scenario labels, to be displayed
        plot_characs        list of characeristics to be plotted. If None, then the default list from databook is used
        
    Notes: 
        TODO: replace the default list of plot_characs as a setting value. 
    
    """
    
    
    
    # close all remaining windows
    pl.close("all") 
    # setup    
    charac_specs = settings.charac_specs
    plotdict = settings.plot_settings
    year_inc = 5.  # TODO: move this to setting
    tvec = scen_results[0].sim_settings['tvec'] # TODO
    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]
    yr_range = np.arange(start_year,end_year+0.1,year_inc,dtype=int)    
    
    
    if plot_charac is None:
        pass
    
    if plot_pops is None:
        plot_pids = getPIDs(scen_results[0],plot_pops) #####
    else:
        plot_pids = range(len(scen_results[0].m_pops))
        plot_pops = [pop.label for pop in scen_results[0].m_pops]
    
    # generate plots
    for (i, pid) in enumerate(plot_pids):
        plot_label = plot_pops[i]
        
        for charac in plot_charac:
        
            yvals = []
            tvals = []
            labels= []
            observed_data = []
            yhat = []
            that = []
            
            for (i,result_name) in enumerate(scen_results.keys()):
                result = scen_results[result_name] ############### GET VALUES JUST FOR THIS POPULATION
                y_values_cur, t_values_cur, final_dict_cur, _ = extractCharacteristic(results=result, charac_label=charac, charac_specs=charac_specs, data=data)
                yvals.append(y_values_cur[pid])
                tvals.append(t_values_cur[pid])
                labels.append(scen_labels[i])
            
            if plot_observed_data:
                pass # TODO: include 
            
            unit_tag = ''
            if 'plot_percentage' in charac_specs[charac].keys():
                ### we don't vals *= 100, as this is already done in extractCharacteristic()
                unit_tag = ' (%)'
            
            final_dict = {'y_hat': yhat,
                  't_hat': that,
                  'ylim' : 0,
                  'unit_tag': unit_tag,
                  'xlabel':'Year',
                  'ylabel': charac_specs[charac]['name'] + unit_tag,
                  'x_ticks' : (yr_range,yr_range),
                  'title': 'Scenario comparison:\n%s [%s]' % (charac_specs[charac]['name'],plot_label),
                  'save_figname': '%s_ScenarioComparision_%s_%s'%(fig_name, plot_label, charac_specs[charac]['name'])}
            final_dict.update(plotdict)
            
            figure = _plotLine(ys = yvals, ts = tvals, labels = labels, legendsettings=None, save_fig=save_fig, fig_name=fig_name, colors=colors, **final_dict)#, y_hat=[final_dict_cur['y_hat'][pid],final_dict_com['y_hat'][pid]], t_hat=[final_dict_cur['t_hat'][pid],final_dict_com['t_hat'][pid]])
            
        
def plotScenarioFlows(scen_results,scen_labels,settings,data,plot_charac=None,plot_pops=None,
                      comp_labels = None, comp_titles = None, pop_labels = None, pop_titles = None, 
                      link_labels = None, include_link_not_exclude = True, link_legend = None, 
                      plot_inflows = True, plot_outflows = True, exclude_transfers = False, sum_total=False,
                      colormappings=None, plot_observed_data=True, save_fig=False, fig_name=None, colors=None):
    """
    Line plots for scenarios. Should be used for flows (instantaneous rates) only.
    
    Params
        scen_results        list of results
        scen_labels         list of scenario labels, to be displayed
        plot_characs        list of characeristics to be plotted. If None, then the default list from databook is used
        
    Notes: 
        TODO: replace the default list of plot_characs as a setting value. 
    
    """
    
    print comp_labels
    
    # close all remaining windows
    pl.close("all") 
    # setup    
    charac_specs = settings.charac_specs
    plotdict = settings.plot_settings
    year_inc = 5.  # TODO: move this to setting
    tvec = scen_results[0].sim_settings['tvec'] # TODO
    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]
    yr_range = np.arange(start_year,end_year+0.1,year_inc,dtype=int)    
    
    
    if plot_charac is None:
        pass
    
    if plot_pops is None:
        plot_pids = getPIDs(scen_results[0],plot_pops)
    else:
        plot_pids = range(len(scen_results[0].m_pops))
        plot_pops = [pop.label for pop in scen_results[0].m_pops]
    
    if pop_labels is None: pop_labels = scen_results[0].pop_labels
    
    if link_legend is None: link_legend = dict()
    
    if comp_labels is None:
        logger.info("No compartments have been selected for flow-plots.")
        comp_labels = []
        
    if comp_titles is not None and len(comp_titles) != len(comp_labels):
        logger.error("Flow-plot failure due to the number of compartment plot titles not matching the number of compartments to analyse.")
    if pop_titles is not None and len(pop_titles) != len(pop_labels):
        logger.error("Flow-plot failure due to the number of population plot titles not matching the number of populations to analyse.")
    
    
    # generate plots
    """
    for (i, pid) in enumerate(plot_pids):
        plot_label = plot_pops[i]
        
        for charac in plot_charac:
        
            yvals = []
            tvals = []
            labels= []
            observed_data = []
            yhat = []
            that = []
    """       
            
    for (i,comp_label) in enumerate(comp_labels):
        
        for (j, pid) in enumerate(plot_pids):
            
            plot_label = plot_pops[j]
            
            yvals = []
            tvals = []
            labels= []
            observed_data = []
            yhat = []
            that = []
            
            for (k,result_name) in enumerate(scen_results.keys()):
                
                result = scen_results[result_name]  
                
                comp = result.m_pops[pid].getComp(comp_label)
                
                all_rates, all_tvecs, all_labels = extractFlows(comp=comp,
                                                                    results=result, 
                                                                    settings=settings,
                                                                    tvec=tvec,
                                                                    link_labels=link_labels,
                                                                    include_link_not_exclude=include_link_not_exclude,
                                                                    link_legend=link_legend,
                                                                    plot_inflows=plot_inflows,
                                                                    plot_outflows=plot_outflows,
                                                                    sum_total=sum_total,
                                                                    exclude_transfers=exclude_transfers)
                
                yvals.append(all_rates[0])
                tvals.append(all_tvecs[0])
                labels.append(scen_labels[k])
        
            if comp_titles is not None:
                title_comp = comp_titles[i]
            else:
                title_comp = 'Compartment: "%s"' % settings.node_specs[comp_label]['name']
            if pop_titles is not None:
                title_pop = plot_pops[j]
            else:
                title_pop = '\nPopulation: "%s"' % pop.label
            
            final_dict = {'ylim' : 0,
              'xlabel':'Year',
              'ylabel': 'Number of People',
              'x_ticks' : (yr_range,yr_range),
              'title': 'Scenario : %s' % (title_comp+title_pop),
              'save_figname': '%s_ScenarioFlowComparision_%s_%s'%(fig_name,comp_label,plot_label)
              }
            final_dict.update(plotdict)
        
            _plotLine(ys=yvals, ts=tvals, labels=labels, colors=colors, save_fig=save_fig, **final_dict)
    
        
               
def plotPopulation(results, data, pop_labels, title='',colormappings=None, 
                   plot_observed_data=True, plot_observed_label="alive", save_fig=False, fig_name=None, 
                   use_full_labels=True, plot_comp_labels=None, plotdict=None,cat_labels=None):
    """ 
    
    Plot all compartments for all populations
    
    Params:
        results         Results object
        title           String, for title of plots
        colormapping    odict with comp.labels as keys and corresponding rgba colors as values
        plotObservedData    boolean flag: whether observed data should be plotted
        save_fig        boolean flag, whether to save figure as "Full_compartment_<pop_label>"
        use_full_labels     boolean flag: whether to use full compartment label (i.e. 'Susceptible') or short label version (i.e. 'sus')
                            The latter is easier for debugging purposes
        
    """
    # setup data structures
    tvec = results.sim_settings['tvec']
    year_inc = 5.  # TODO: move this to setting
    if plotdict is not None and 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]
    
    yr_range = np.arange(start_year,end_year+0.1,year_inc,dtype=int)
    mpops = results.m_pops
    sim_settings = results.sim_settings
    save_figname=None
    dataobs = None # placeholder for observed data
    
    if use_full_labels:
        comp_labels = results.comp_label_names
        ncol = 1
    else:
        ncol = 2
    
    # setup: determine compartment indices to be plotted
    if plot_comp_labels is None:
        plot_comp_labels = sorted(mpops[0].comp_ids, key=mpops[0].comp_ids.get)
        
    comp_indices = [mpops[0].comp_ids[i] for i in plot_comp_labels]
    
    
    # setup: determine colors to be used
    colors = []
    if colormappings is not None:
        colors_dict, cat_colors = getCategoryColors(colormappings,plotdict['colormapping_order'])
        # reorder so that colors are same as expected for plotting the population
        for (j,comp_label) in enumerate(plot_comp_labels):
            if isPlottableComp(comp_label,sim_settings,results.comp_specs):
                colors.append(colors_dict[comp_label])
    
    # setup: plotting dict structures  
    if plotdict is None:
        plotdict = {}
        
        
    # iterate for each key population group
    for (i,pop) in enumerate(mpops):
        
        if pop.label not in pop_labels:
            continue
        
        yvals, labels, dataobs = extractCompartment(results, sim_settings, pop, comp_indices, data, comp_labels=comp_labels, plot_observed_data=plot_observed_data, 
                                           plot_observed_label=plot_observed_label, use_full_labels=use_full_labels)
        
        xlim = ()
        
        pl_title = title+' Population: %s' % (pop.label)
        if save_fig:
            save_figname = fig_name + "_compartment_%s"%pop.label
        dict = {  'ymin': 0,
                  'xlabel': 'Year',
                  'year_inc' :  5.,
                  'ylabel': 'People',
                  'mec' : 'k',
                  'x_ticks' : (yr_range,yr_range),
                  'colors': colors
                  }
        dict.update(plotdict)
    
        legendsettings =  {'loc':'center left', 
                           'bbox_to_anchor':(1.05, 0.5), 
                           'ncol':ncol}
   
        _plotStackedCompartments(tvec, yvals, labels,datapoints=dataobs,title=pl_title,legendsettings=legendsettings, catlabels=cat_labels,catcolors=cat_colors,
                                     save_fig=save_fig,save_figname=save_figname,**dict)
        
        
    if dict.has_key('legend_off') and dict['legend_off']:
        # Do this separately to main iteration so that figure not corrupted
        # Note that colorlist may be different to colors, as it represents 
        # classes of compartments
        separateLegend(labels=cat_labels,colors=cat_colors,fig_name=fig_name)
        
         

def plotCharacteristic(results, charac_specs, data, title='', outputIDs=None, 
                       pop_labels = None, plot_total = False,
                       plot_observed_data=True, save_fig=False, fig_name=None, 
                       colors=None, plotdict=None):
    """
    Plot a characteristic across all populations
    
    Params:
        results
        charac_specs
        title
        outputIDs        list of compartment labels (characs.keys()) which will be selectively be plotted
                         Default: None, which causes all labels to be plotted
        pop_ids          list of labels for subset of populations to be plotted. If None, all populations are plotted
        plot_total       sum and plot the total. 
        plotObservedData plot observed data points on top of simulated data. Useful for calibration
        save_fig         bool to indicate whether to save the figure to file
        colors           list of colors for populations
        
        
    Example dict:
        dict = {'y_hat': yhat,
                't_hat': that,
                'unit_tag': '(%)', # default = ''
                'xlabel':'Year',
                'ylabel': charac_specs[output_id]['name'] + '(%)',
                'title': '%s Outputs - %s' % (title, charac_specs[output_id]['name']),
                'marker': 'o',
                'x_ticks' : ([2000,2030],[2000,2030]),
                'save_figname': 'MyPlot'}
                
    
        
    """
    # setup
    outputs = results.outputs
    
#    legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':1}
    
    if outputIDs is None:
        outputIDs = outputs.keys()
           
        
    for output_id in outputIDs:
        if isPlottableCharac(output_id, charac_specs):
            y_values, t_values, final_dict, pop_labels = extractCharacteristic(results=results, charac_label=output_id, 
                                                                   charac_specs=charac_specs, data=data, plot_observed_data=plot_observed_data,
                                                                   pop_labels = pop_labels, plot_total = plot_total,
                                                                   fig_name=fig_name, plotdict=plotdict)
            
            _plotLine(y_values, t_values, pop_labels, legendsettings=None, save_fig=save_fig, colors=colors, **final_dict)




def plotSingleCompartmentFlow(results, settings, comp_labels = None, comp_titles = None, plot_pops = None, pop_labels = None, pop_titles = None, 
              link_labels = None, include_link_not_exclude = True, link_legend = None, sum_total=False,
              plot_inflows = True, plot_outflows = True, exclude_transfers = False, observed_data = None,
              save_fig=False, fig_name=None, colors=None, suppress_plot=False):
    """
    Plot flows rates in and out of a compartment.
    """
    plotdict = settings.plot_settings
    year_inc = 5. # TODO remove hardcoded ref
    tvec = results.sim_settings['tvec']
    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]
    yr_range = np.arange(start_year,end_year+0.1,year_inc,dtype=int)    
    
    print pop_labels, plot_pops
    
    
    if pop_labels is None:
        pop_labels = results.pop_labels
    
    
    if link_legend is None: link_legend = dict()
    
    if plot_pops is None:
        plot_pids = getPIDs(results,pop_labels)
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
    
    
    for (i,comp_label) in enumerate(comp_labels):
        
        
        for (j, pid) in enumerate(plot_pids):
            
            plot_label = plot_pops[j]
            
            comp = results.m_pops[pid].getComp(comp_label)
    #            print comp_label
    #            print 'out'
    #            print comp.outlink_ids
    #            print 'in'
    #            print comp.inlink_ids

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
                title_pop = '\nPopulation: "%s"' % pop.label
        
            
            final_dict = {
              'ylim' : 0,
              'xlabel':'Year',
              'ylabel': 'Number of People',
              'x_ticks' : (yr_range,yr_range),
              'title': title_comp+title_pop,
              'save_figname': '%s_FlowComparision_%s_%s'%(fig_name,comp_label,plot_label)
              }
            
            if observed_data is not None:
                final_dict['y_hat'] = [observed_data[0]]
                final_dict['t_hat'] = [observed_data[1]]
            final_dict.update(plotdict)
            
            
            if len(all_rates) > 0:
                _plotLine(ys=all_rates, ts=all_tvecs, labels=all_labels,colors=colors,save_fig=save_fig, **final_dict)
            else:
                logger.warn("No flows selected for plotting")
            

def plotPopulationFlows(results, settings, comp_labels = None, comp_titles = None, plot_pops = None, pop_labels = None, pop_titles = None, 
              link_labels = None, include_link_not_exclude = True, link_legend = None, sum_total=False, sum_population=False,
              plot_inflows = True, plot_outflows = True, exclude_transfers = False, observed_data = None,
              save_fig=False, fig_name=None, colors=None, suppress_plot=False):
    """
    Plot flows rates in and out of a compartment, across populations. 
    
    Suitable for total new infections, deaths, etc. where the net flow is required. 
    """
    plotdict = settings.plot_settings
    year_inc = 5. # TODO remove hardcoded ref
    tvec = results.sim_settings['tvec']
    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[1]
    yr_range = np.arange(start_year,end_year+0.1,year_inc,dtype=int)    
    
    print pop_labels, plot_pops
    
    
    if pop_labels is None:
        pop_labels = results.pop_labels
    
    
    if link_legend is None: link_legend = dict()
    
    if plot_pops is None:
        plot_pids = getPIDs(results,pop_labels)
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
    
    
    for (i,comp_label) in enumerate(comp_labels):
        
        all_rates = []
        all_tvecs = []
        
        
        for (j, pid) in enumerate(plot_pids):
            
            plot_label = plot_pops[j]
            
            comp = results.m_pops[pid].getComp(comp_label)
    #            print comp_label
    #            print 'out'
    #            print comp.outlink_ids
    #            print 'in'
    #            print comp.inlink_ids

            rates, tvecs, all_labels = extractFlows(comp=comp,
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
            
            print len(rates), len(tvecs)
        
            
            all_rates.append(rates)
            all_tvecs.append(tvecs)
            
        if comp_titles is not None:
            title_comp = comp_titles[i]
        else:
            title_comp = 'Compartment: "%s"' % settings.node_specs[comp_label]['name']
        if pop_titles is not None:
            title_pop = plot_pops[j]
        else:
            title_pop = '\nPopulation: "%s"' % pop.label
    
        
        final_dict = {
          'ylim' : 0,
          'xlabel':'Year',
          'ylabel': 'Number of People',
          'x_ticks' : (yr_range,yr_range),
          'title': title_comp,
          'save_figname': '%s_FlowComparision_%s'%(fig_name,comp_label)
          }
            
        if observed_data is not None:
            final_dict['y_hat'] = [observed_data[0]]
            final_dict['t_hat'] = [observed_data[1]]
        final_dict.update(plotdict)
        
        print all_rates
        print all_tvecs
        
        if len(all_rates) > 0:
            _plotLine(ys=all_rates, ts=all_tvecs, labels=all_labels,colors=colors,save_fig=save_fig, **final_dict)
        else:
            logger.warn("No flows selected for plotting")
            



def extractCompartment(results, sim_settings, pop, comp_indices, data, comp_labels=None, plot_observed_data=True, plot_observed_label="alive", 
                   use_full_labels=False):
    
    comps = []
    labels = []
    dataobs = None
        
    for comp_id in comp_indices:
        
        comp_label = pop.comps[comp_id].label
        
        if isPlottableComp(comp_label,sim_settings,results.comp_specs):
        
            comps.append(pop.comps[comp_id])
            
            if use_full_labels:
                c_label = comp_labels[comp_label]
            else:
                c_label = comp_label
            labels.append(c_label)
            
    if plot_observed_data:
        ys = data['characs'][plot_observed_label][pop.label]['y']
        ts = data['characs'][plot_observed_label][pop.label]['t']
        dataobs = (ts,ys)

    return comps, labels, dataobs

    

def extractCharacteristic(results, charac_label, charac_specs, data, title='', 
                          pop_labels = None, plot_total = False, plot_observed_data=True, fig_name=None, plotdict=None):
    
    if plotdict is None: plotdict = {}    
    tvec = results.sim_settings['tvec']
    year_inc = 5.  # TODO: move this to setting
    if 'xlim' in plotdict.keys():
        xlim = plotdict['xlim']
        start_year, end_year = xlim[0], xlim[1]
    else:
        start_year, end_year = tvec[0], tvec[-1]
    yr_range = np.arange(start_year,end_year+0.1,year_inc,dtype=int)    
      
    if pop_labels is None:
        pop_labels = [pop.label for pop in results.m_pops]
        
    outputs = results.outputs
    sim_settings = results.sim_settings
    unit_tag = ''
    
    output_id = charac_label
    
    y_values = []
    t_values = []
    yhat = []
    that = []
    
    for k,poplabel in enumerate(pop_labels):
        
        vals = dcp(outputs[output_id][poplabel])
        
        if 'plot_percentage' in charac_specs[output_id].keys():
            vals *= 100
            unit_tag = ' (%)'
        y_values.append(vals)
        t_values.append(sim_settings['tvec'])
        
        if plot_observed_data:
            if output_id in data['characs'].keys():
                ys = data['characs'][output_id][poplabel]['y']
                ts = data['characs'][output_id][poplabel]['t']
            else:   # For the case when plottable characteristics were not in the databook and thus not converted to data.
                ys = []
                ts = []
            if 'plot_percentage' in charac_specs[output_id].keys():
                ys *= 100
#            if len(ys)==0:
#                # add an empty list to preserve label colours
#                ys,ts = [],[]
            yhat.append(ys)
            that.append(ts)
            
    if plot_total:
        y_values = np.array(y_values)
        y_values = [y_values.sum(axis=0)]
        t_values = [t_values[0]]
        # enforce so that no observed data is included (as we can't enforce that all points are supplied).
        # This can be circumvented i.e. plotting against totals, by updating the values in 
        # the returned plotting dictionary
        yhat = [[]]
        that = [[]]
        pop_labels = ['Total']
            
    final_dict = {'y_hat': yhat,
                  't_hat': that,
                  'unit_tag': unit_tag,
                  'xlabel':'Year',
                  'ylabel': charac_specs[output_id]['name'] + unit_tag,
                  'x_ticks' : (yr_range,yr_range),
                  'title': '%s Outputs: %s' % (title, charac_specs[output_id]['name']),
                  'save_figname': '%s_characteristic_%s'%(fig_name, charac_specs[output_id]['name'])}
    final_dict.update(plotdict)
    
    return y_values, t_values, final_dict, pop_labels
        
      
def extractFlows(comp, results, settings, tvec, link_labels = None, include_link_not_exclude = True, link_legend = None, 
                  plot_inflows = True, plot_outflows = True, exclude_transfers = False, sum_total=False):
    all_labels = []
    all_rates = []
    all_tvecs = []
    for in_out in xrange(2):
        if (in_out == 0 and plot_inflows) or (in_out == 1 and plot_outflows):
            comp_link_ids = [comp.inlink_ids, comp.outlink_ids][in_out]
            label_tag = ['In: ','Out: '][in_out]
            for link_tuple in comp_link_ids:
                link = results.m_pops[link_tuple[0]].links[link_tuple[1]]
                #print link.label, link_labels, include_link_not_exclude
                if link_labels is None or (include_link_not_exclude and link.label in link_labels) or (not include_link_not_exclude and link.label not in link_labels):
                    try: 
                        legend_label = label_tag + settings.linkpar_specs[link.label]['name']
                    except: 
                        if exclude_transfers: continue
                        else: legend_label = label_tag + link.label
                    if link.label in link_legend:
                        legend_label = link_legend[link.label]    # Overwrite legend for a label.
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
                            num_flow = 1 - (1 - num_flow) ** results.dt     # Fractions must be converted to effective timestep rates.
                            num_flow *= comp_source.popsize
                        num_flow /= results.dt      # All timestep-based effective fractional rates must be annualised.
                        
                    all_labels.append(legend_label)
                    all_rates.append(num_flow)
                    all_tvecs.append(tvec)
                    
    if sum_total:
        all_tvecs = all_tvecs[:1]
        all_labels = ['Total summed movement']
        all_rates_tmp = np.array(all_rates)
        all_rates = [all_rates_tmp.sum(axis=0)]
        
    return all_rates, all_tvecs, all_labels
    
     
def _plotStackedCompartments(tvec,comps,labels=None,datapoints=None,title='',ylabel=None,xlabel=None,xlim=None,ymin=None,ylim=None,save_figname=None,legend_off=False,
                            save_fig=False,colors=None,catlabels=None,catcolors=None,
                            marker='o',edgecolors='k',facecolors='none',s=40,zorder=10,linewidth=3,x_ticks=None,legendsettings={},**kwargs):  
    """ 
    Plot compartments in time. 
    This creates a stacked plot of several compartments, over a given period.
    Observed data points can also be additionally plotted, overlaying the data
        
    Params:
        tvec        time period
        comps       compartment sizes
        labels      list of labels
        datapoints  observed datapoints, specified as a list of tuples 
        **kwargs    further keyword arguments, such as ylims
    """
    if colors is None or len(colors) != len(comps):
        if len(colors) != len(comps):
            logger.info("Plotting: setting color scheme to be default colormap, as not all compartments had color assigned")
        colors = gridColorMap(len(comps))
    
    # setup
    fig, ax = pl.subplots()
    bottom = 0*tvec
    max_val = 0
    
    for (k,comp) in enumerate(comps):
        top = bottom + comp.popsize
        lw = 0
        if save_fig:
            lw = 0.1
        ax.fill_between(tvec, bottom, top, facecolor=colors[k], alpha=1,lw=lw) # for some reason, lw=0 leads to no plot if we then use fig.savefig()
        reg, = ax.plot((0, 0), (0, 0), color=colors[k], linewidth=10)
        bottom = dcp(top)
        
    max_val = max(top)
        
    if datapoints is not None:
        ts,ys= datapoints[0],datapoints[1]
        ax.scatter(ts,ys,marker=marker,edgecolors=edgecolors,facecolors=facecolors,s=s,zorder=zorder,linewidth=linewidth)
        if max(ys) > max_val:
            max_val = max(ys) 
            
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])   
    
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.set_ylim(ymax=max_val*1.05)
    if ylim is not None:
        ax.set_ylim(ylim)
    elif ymin is not None:
        ax.set_ylim(ymin=ymin)
    
    if x_ticks is not None:
        ax.set_xticks(x_ticks[0])
        ax.set_xticklabels(x_ticks[1])
    
    if not legend_off:
        ax.legend(labels, **legendsettings)
    
    _turnOffBorder()
    pl.suptitle('')
    
    if save_fig:
        fig.savefig('%s' % (save_figname))
        logger.info("Saved figure: '%s'"%save_figname)
        
    
 
    
def _plotLine(ys,ts,labels,colors=None,y_hat=[],t_hat=[],
             legendsettings=None,title=None,xlabel=None,ylabel=None,xlim=None,ylim=None,y_ticks=None,x_ticks=None,
             marker='o',s=40,facecolors='none',linewidth=3,zorder=10,save_fig=False,save_figname=None,legend_off=False,**kwargs):
    """
    
    """
    
    if legendsettings is None: legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':1}    
    
    ymin_val = np.min(ys[0])
    
    if colors is None:        
        colors = gridColorMap(len(ys))  
    elif len(colors) < len(ys):
        colors = gridColorMap(len(ys))    
        
    fig, ax = pl.subplots()
    print len(ys), len(ts), len(colors)
    for k,yval in enumerate(ys):
        
        ax.plot(ts[k], yval, c=colors[k])
        if np.min(yval) < ymin_val:
            ymin_val = np.min(yval)
            
        if len(y_hat) > 0 and len(y_hat[k]) > 0: # i.e. we've seen observable data
            ax.scatter(t_hat[k],y_hat[k],marker=marker,edgecolors=colors[k],facecolors=facecolors,s=s,zorder=zorder,linewidth=linewidth)
            if np.min(y_hat[k]) < ymin_val:
                ymin_val = np.min(y_hat[k])
        
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])   
    
          
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    
    if not legend_off:
        ax.legend(labels, **legendsettings)
    
    ymin = ax.get_ylim()[0]
    # Set the ymin to be halfway between the ymin_val and current ymin. 
    # This seems to get rid of the worst of bad choices for ylabels[0] = -5 when the real ymin=0
    tmp_val = (ymin+ymin_val)/2.
    ax.set_ylim(ymin=tmp_val)
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
        logger.info("Saved figure: '%s'"%save_figname)
        
    return fig
    

        
              
def plotAllOutflows(results, num_subplots = 5):
    """ 
    Visualise outflows for each compartment in each population as fractions of compartment size
    """
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
                fig, ax = pl.subplots(nrows = num_subplots, sharex = True)
                pl.suptitle('Population (%i): %s' % (pid, pop.label.title()))
                ax[num_subplots-1].set_xlabel('Year')
            bottom = 0*sim_settings['tvec']
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
            ax[plot_id].set_position([box.x0, box.y0, box.width*0.75, box.height])
            ax[plot_id].legend(patch_labels,**legendsettings)

            cid += 1
        
        while cid % num_subplots != 0:
            plot_id = cid % num_subplots
            box = ax[plot_id].get_position()
            ax[plot_id].set_position([box.x0, box.y0, box.width*0.75, box.height])
            cid += 1
            
            
        pid += 1
        
        
