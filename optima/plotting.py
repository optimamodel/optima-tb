#%% Imports




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