
from analysis import evaluateDiseaseProgression




import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

"""
Run model simulation and plot output:

Loads and runs model data

"""

from project import Project
import pylab

proj = Project(name = 'Belarus', cascade_path = './cascade-belarus.xlsx', validation_level = 'avert')

#set the year range we simulate over as starting in 1995:
proj.setYear([2000,2030],False)


proj.loadSpreadsheet(databook_path = './databook-belarus-template.xlsx')
proj.makeParset(name = 'default')

# ================================================
# Define populations and progressions:
# ------------------------------------------------

specified_populations = proj.parsets[0].pop_labels # this outputs all populations

# To evaluate given transitions between compartments, update 'specified_progressions':
specified_progressions = {'sus': ['lteu','ltlu','spdu','spmu'],
                          'lteu': ['ltet','ltlu','spdu','spmu'],
                          'spdu': ['spdd']}

year_track = [1.,2.,3.] # times that we report on. Note that these can only be multiple of dt
outputfile = "BelarusDiseaseProgression.csv"
# ================================================

evaluateDiseaseProgression(proj,
                           specified_progressions=specified_progressions,
                           specified_populations=specified_populations,
                           year_track=year_track,
                           birth_transit='b_rate',
                           output_file=outputfile)

