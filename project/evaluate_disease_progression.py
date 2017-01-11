
import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from analysis import evaluateDiseaseProgression


"""
Determines disease progression from certain compartments to other compartments, 
by (artifically) setting the population of a compartment to a predefined size
(default size=1e6) and then letting the disease run. Note that birth rate is set
to zero, so the total population size (including deaths) should be constant.

This process outputs: 
 1) population plots showing the progress of the disease, for all compartments
 2) a csv file, which for each population outputs certain transitions from one
     compartment to another compartment, expressed as a percentage. The times that
     are returned can be defined as required (as year_track), and can include partial 
     years (i.e. 0.5, provided that they are a multiple of dt)
     

This is done for all populations, and for a list of specific mappings between compartments
i.e. from sus to lteu. Multiple mappings can also be included too i.e. sus to lteu and ltlu.

@author: sjjarvis
@date: 11/01/2017

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

