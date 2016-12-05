#%% Imports
import logging
logger = logging.getLogger(__name__)

import pylab as pl
from copy import deepcopy as dcp


   
class Plotter():
        
    
    def __init__(self,data):
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
        
        
    
    def plotProjectResults(self,results,outputs,sim_settings,charac_specs,title='',debug=False):
        """
        Placeholder plotting function, originally located in project.py
        """
        pl.close("all") 
        self.plotPopulation(results, outputs, sim_settings, charac_specs, title)
        
        # List charac labels for those that were present in databook (i.e. could have cascade-provided defaults overwritten).
        label_list = [x for x in charac_specs.keys() if charac_specs[x]['databook_order'] >= 0]
        self.plotCompartment(results, outputs, sim_settings, charac_specs, title, outputIDs = label_list)
        
        if debug:
            self.plotOutflows(results, sim_settings)
        
    def plotPopulation(self,results,outputs,sim_settings,charac_specs,title='',plotObservedData=True,saveFig=False):
        """ 
        
        Plot all compartments for a population
        """
        
        for pop in results:
               
            fig, ax = pl.subplots()
            colors = self.gridColorMap(len(pop.comps))
            bottom = 0*sim_settings['tvec']
            
            k = 0
            labels = []
            for comp in pop.comps:
                
                if not self.isPlottable(comp.label,sim_settings):
                    continue
                
                top = bottom + comp.popsize
                labels.append(comp.label)  #explicitly gather comp.label names as not all cascade names are going to be plotted
                ax.fill_between(sim_settings['tvec'], bottom, top, facecolor=colors[k], alpha=1,lw=0) # for some reason, lw=0 leads to no plot
                ax.plot((0, 0), (0, 0), color=colors[k], linewidth=10)
                bottom = dcp(top)
                k += 1
            
            if plotObservedData:
                # TODO confirm that alive will always be tag. It probably won't be, so better to have this named in the settings? userdefined?
                try:
                    ys = self.data['characs']['alive'][pop.label]['y']
                    ts = self.data['characs']['alive'][pop.label]['t']
                    ax.scatter(ts,ys,marker='o',edgecolors='k',facecolors='none',s=40,zorder=10,linewidth=3)
                except:
                    logger.info("No observed data with label='alive' for plotting in plotPopulations() [plotting.py]")
                
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])   
            
            legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':2}
            ax.set_title('%s Cascade - %s' % (title, pop.label.title()))
            ax.set_xlabel('Year')
            ax.set_ylabel('People')
            ax.set_xlim((sim_settings['tvec'][0], sim_settings['tvec'][-1]))
            ax.set_ylim((0, max(top)))
            ax.legend(labels,**legendsettings)
            self.turnOffBorder()
            if saveFig:
                fig.savefig('%sCascade-%s.png' % (title, pop.label.title()))
            
            
    
    def plotCompartment(self,results,outputs,sim_settings,charac_specs,title='',outputIDs=None,plotObservedData=True,saveFig=False):
        """
        Plot a compartment across all populations
        
        Params:
            results
            outputs
            sim_settings
            charac_specs
            title
            outputIDs        list of compartment labels (characs.keys()) which will be selectively be plotted
                             Default: None, which causes all labels to be plotted
            plotObservedData plot observed data points on top of simulated data. Useful for calibration
        """
        if outputIDs is None:
                outputIDs = outputs.keys()
                
        colors = self.gridColorMap(len(results))        
        for output_id in outputIDs:
            
            if not self.isPlottableCharac(output_id,charac_specs):
                continue
            
            
            unit_tag = ''
            fig, ax = pl.subplots()
            for k,pop in enumerate(results):
                vals = dcp(outputs[output_id][pop.label])
                if 'plot_percentage' in charac_specs[output_id].keys():
                    vals *= 100
                    unit_tag = ' (%)'
                ax.plot(sim_settings['tvec'], vals,c=colors[k])
                
                if plotObservedData:
                    ys = self.data['characs'][output_id][pop.label]['y']
                    if len(ys)==0 and ys[0]==0:
                        # don't bother plotting if there's nothing to plot
                        # TODO better methods of indicating absent data?
                        continue
                    if 'plot_percentage' in charac_specs[output_id].keys():
                        ys *= 100
                    ts = self.data['characs'][output_id][pop.label]['t']
                    ax.scatter(ts,ys,marker='o',edgecolors=colors[k],facecolors='none',s=40,zorder=10,linewidth=3)
                    
    
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])   
            
            legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':1}                
            ax.set_title('%s Outputs - %s' % (title, charac_specs[output_id]['name']))
            ax.set_xlabel('Year')
            ax.set_ylabel(charac_specs[output_id]['name'] + unit_tag)
            ax.legend([pop.label for pop in results], **legendsettings)
            self.turnOffBorder()
            if saveFig:
                fig.savefig('%sOutputs-%s.png' % (title, charac_specs[output_id]['name']))                    

                
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