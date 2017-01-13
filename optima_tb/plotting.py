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
   
   
class Plotter():
        
    
    def __init__(self,data={}):
        self.updateData(data)
        
    def updateData(self,data):
        self.data = data
    
    #%% Function to generate useful colormap for plots
    
    def gridColorMap(self,ncolors=10, limits=None, nsteps=10, asarray=False, doplot=False, newwindow=True):
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
    
    
    def _getColormapRange(self,cmap,ncols=10,order='alternate'):
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
            a = order_map[ncols/2:]
            b = order_map[:ncols/2]
            order_map[::2] = a 
            order_map[1::2] = b 
        else: # order == 'sequential', so leave in order 
            pass
        
        return [scalar_map.to_rgba(index) for index in order_map]
    
    
    def getCategoryColors(self,category_list,order='alternate'):
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
        for k,v in category_list.iteritems():
            tmp_list = self._getColormapRange(k,len(v))
            for i,label in enumerate(v):
                col_list[label] = tmp_list[i]
        return col_list
        
        
    
    def turnOffBorder(self):
        
        pl.gca().spines['right'].set_color('none')
        pl.gca().spines['top'].set_color('none')
        pl.gca().xaxis.set_ticks_position('bottom')
        pl.gca().yaxis.set_ticks_position('left')
    
    
    def isPlottable(self,comp_label,sim_settings):
        """ 
        Returns bool indicating whether a population label should be included in metrics
        for population reporting 
        """
        if comp_label in sim_settings['tag_no_plot']:
            return False
        return True
    
    def isPlottableCharac(self,output_id,charac_specs):
        """
        Returns False if specified in cascade spreadsheet as 'n' or 'N' for plot characteristic.
        Else returns True
        """
        try: # Better to try and ask for forgiveness than for permission ... 
            if charac_specs[output_id]['plot_characteristic'].lower() == 'n':
                return False
        except: 
            return True
        
        
    
    def plotProjectResults(self,results,charac_specs,title='',colormappings=None,debug=False):
        """
        Placeholder plotting function, originally located in project.py
        """
        pl.close("all") 
        
        
        self.plotPopulation(results, title, colormappings=colormappings)
        
        # List charac labels for those that were present in databook (i.e. could have cascade-provided defaults overwritten).
        label_list = [x for x in charac_specs.keys() if charac_specs[x]['databook_order'] >= 0]
        self.plotCompartment(results, charac_specs, title, outputIDs = label_list)
        
        if debug:
            self.plotOutflows(m_pops, sim_settings)
        
        
        
    def plotPopulation(self,results,title='',colormappings=None,plotObservedData=True,save_fig=False,use_full_labels=True,**kwargs):
        """ 
        
        Plot all compartments for all populations
        
        Params:
            results         Results objet
            title           String, for title of plots
            colormapping    odict with comp.labels as keys and corresponding rgba colors as values
            plotObservedData    boolean flag: whether observed data should be plotted
            save_fig        boolean flag, whether to save figure as "Full_compartment_<pop_label>"
            use_full_labels     boolean flag: whether to use full compartment label (i.e. 'Susceptible') or short label version (i.e. 'sus')
                                The latter is easier for debugging purposes
            
        """
        # setup
        tvec = results.sim_settings['tvec']
        mpops = results.m_pops
        sim_settings = results.sim_settings
        pop_labels = results.pop_labels
        save_figname=None
        
        if use_full_labels:
            comp_labels = results.comp_label_names
            ncol = 1
        else:
            ncol = 2
        
        colors = []
        if colormappings is not None:
            colors_dict = self.getCategoryColors(colormappings)
            # reorder so that colors are same as expected for plotting the population
            for (j,comp) in enumerate(mpops[0].comps):
                if self.isPlottable(comp.label,sim_settings):
                    colors.append(colors_dict[comp.label])
        
            
        # iterate for each key population group
        for (i,pop) in enumerate(mpops):
            comps = []
            labels = []
            for (j,comp) in enumerate(pop.comps):
                if self.isPlottable(comp.label,sim_settings):
                    comps.append(comp)
                    if use_full_labels:
                        c_label = comp_labels[comp.label]
                    else:
                        c_label = comp.label
                    labels.append(c_label)
            if plotObservedData:
                ys = self.data['characs']['alive'][pop.label]['y']
                ts = self.data['characs']['alive'][pop.label]['t']
                dataobs = (ts,ys)
            xlim = ()
            
            pl_title = 'Compartments for Population: %s' % (pop_labels[i])
            if save_fig:
                save_figname = 'Full_compartment_%s'%pop_labels[i]
            kwargs = {'xlim': (tvec[0],tvec[-1]),
                      'ymin': 0,
                      'xlabel': 'Year',
                      'ylabel': 'People',
                      'mec' : 'k',
                      'marker':'o',
                      's':40,
                      'linewidth':3,
                      'facecolors':'none',
                      'colors': colors
                      }
            legendsettings =  {'loc':'center left', 
                               'bbox_to_anchor':(1.05, 0.5), 
                               'ncol':ncol}
       
            self._plotStackedCompartments(tvec, comps, labels,datapoints=dataobs,title=pl_title,legendsettings=legendsettings, 
                                         save_fig=save_fig,save_figname=save_figname,**kwargs)
            
          
    def _plotStackedCompartments(self,tvec,comps,labels=None,datapoints=None,title='',ylabel=None,xlabel=None,xlim=None,ymin=None,ylim=None,save_figname=None,
                                save_fig=False,colors=None,marker='o',edgecolors='k',facecolors='none',s=40,zorder=10,linewidth=3,legendsettings={},**kwargs):  
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
            colors = self.gridColorMap(len(comps))
        
        # setup
        fig, ax = pl.subplots()
        bottom = 0*tvec
        max_val = 0
        
        for (k,comp) in enumerate(comps):
            top = bottom + comp.popsize
            ax.fill_between(tvec, bottom, top, facecolor=colors[k], alpha=1,lw=0) # for some reason, lw=0 leads to no plot
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
        ax.set_xlim(xlim)
        ax.set_ylim(ymax=max_val*1.05)
        if ylim is not None:
            ax.set_ylim(ylim)
        elif ymin is not None:
            ax.set_ylim(ymin=ymin)
        
        ax.legend(labels,**legendsettings)
        self.turnOffBorder()
        pl.suptitle('')
        
        if save_fig:
            fig.savefig('%s_stacked.png' % (save_figname))
            logger.info("Saved figure: '%s_stacked.png'"%save_figname)
            
        
    
    def plotCompartment(self,results,charac_specs,title='',outputIDs=None,plotObservedData=True,save_fig=False,colors=None):
        """
        Plot a compartment across all populations
        
        Params:
            results
            charac_specs
            title
            outputIDs        list of compartment labels (characs.keys()) which will be selectively be plotted
                             Default: None, which causes all labels to be plotted
            plotObservedData plot observed data points on top of simulated data. Useful for calibration
            save_fig         bool to indicate whether to save the figure to file
            colors           list of colors for populations
            
            
        Example dict:
            dict = {'y_hat': yhat,
                    't_hat': that,
                    'unit_tag': unit_tag,
                    'xlabel':'Year',
                    'ylabel': charac_specs[output_id]['name'] + unit_tag,
                    'title': '%s Outputs - %s' % (title, charac_specs[output_id]['name']),
                    'marker': 'o',
                    'x_ticks' : ([2000,2030],[2000,2030]),
                    'save_figname': 'MyPlot'}
            
        """
        # setup
        tvec = results.sim_settings['tvec']
        year_inc = 5.  # TODO: move this to setting
        yr_range = np.arange(tvec[0],tvec[-1]+0.1,year_inc,dtype=int)
        mpops = results.m_pops
        sim_settings = results.sim_settings
        pop_labels = results.pop_labels
        outputs = results.outputs
        save_figname=None    
        unit_tag = ''
        
        if outputIDs is None:
            outputIDs = outputs.keys()
            
        for output_id in outputIDs:
            if not self.isPlottableCharac(output_id,charac_specs):
                continue
            
            y_values = []
            t_values = []
            yhat = []
            that = []
            
            for k,pop in enumerate(mpops):
                
                vals = dcp(outputs[output_id][pop.label])
                if 'plot_percentage' in charac_specs[output_id].keys():
                    vals *= 100
                    unit_tag = ' (%)'
                y_values.append(vals)
                t_values.append(sim_settings['tvec'])
                
                if plotObservedData:
                    ys = self.data['characs'][output_id][pop.label]['y']
                    ts = self.data['characs'][output_id][pop.label]['t']
                    if 'plot_percentage' in charac_specs[output_id].keys():
                        ys *= 100
                    if len(ys)==0 and ys[0]==0:
                        # add an empty list to preserve label colours
                        ys,ts = [],[]
                    yhat.append(ys)
                    that.append(ts)
                
                    
            dict = {'y_hat': yhat,
                    't_hat': that,
                    'unit_tag': unit_tag,
                    'xlabel':'Year',
                    'ylabel': charac_specs[output_id]['name'] + unit_tag,
                    'title': '%s Outputs - %s' % (title, charac_specs[output_id]['name']),
                    'marker': 'o',
                    'facecolors':'none',
                    's': 40,
                    'linewidth':3,
                    'x_ticks' : (yr_range,yr_range),
                    'save_figname': 'PlotCompartment_%s_%s'%(title, charac_specs[output_id]['name'])}
            legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':1}         
             
            self._plotLine(y_values, t_values, pop_labels, legendsettings=legendsettings,save_fig=save_fig,**dict)
            
        
    def _plotLine(self,ys,ts,labels,colors=None,y_hat=[],t_hat=[],
                 legendsettings=None,title=None,xlabel=None,ylabel=None,xlim=None,ylim=None,y_ticks=None,x_ticks=None,
                 marker='o',s=40,facecolors='none',linewidth=3,zorder=10,save_fig=False,save_figname=None,**kwargs):
        """
        
        """
        if colors is None:        
            colors = self.gridColorMap(len(ys))       
            
        fig, ax = pl.subplots()
        
        for k,yval in enumerate(ys):
            
            ax.plot(ts[k], yval, c=colors[k])
            ax.scatter(t_hat[k],y_hat[k],marker=marker,edgecolors=colors[k],facecolors=facecolors,s=s,zorder=zorder,linewidth=linewidth)
  

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])   
        
              
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(labels, **legendsettings)
        
        if x_ticks is not None:
            ax.set_xticks(x_ticks[0])
            ax.set_xticklabels(x_ticks[1])
        if y_ticks is not None:
            ax.set_yticks(y_ticks[0])
            ax.set_yticklabels(y_ticks[1])
            
        self.turnOffBorder()
        if save_fig:
            fig.savefig('%s.png' % (save_figname))                    
            logger.info("Saved figure: '%s.png'"%save_figname)
        
        
        
                
    def plotOutflows(self, results, sim_settings, num_subplots = 5):
        """ 
        Visualise outflows for each compartment in each population as fractions of compartment size
        """
        
        # NOTE: Hard coding is bad, but results/plots need a lot more refining before propagating stable reference.
        num_links = len(results[0].links)
        
        colors = self.gridColorMap(num_links)          
        legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':2}        

        
        
        pid = 0
        for pop in results:
            cid = 0
            for comp in pop.comps:
                plot_id = cid % num_subplots
                patch_labels = []
                if plot_id == 0:
                    fig, ax = pl.subplots(nrows = num_subplots, sharex = True)
                    pl.suptitle('Population (%i): %s' % (pid, pop.label.title()))
                    ax[num_subplots-1].set_xlabel('Year')
                bottom = 0*sim_settings['tvec']
                for link_id in comp.outlink_ids:
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