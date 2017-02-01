from optima_tb.project import Project
from optima_tb.calibration import performSensitivityAnalysis

"""
Aim: Test and validate Sensitivity analysis
"""

num_pop = 2

databook = '../data/databook-simple-cascade-autocalibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = '../data/cascade-simple-calibration.xlsx')
proj.setYear([2000.,2030.]) 
proj.loadSpreadsheet(databook_path = databook)
proj.makeParset()
FitScore = performSensitivityAnalysis(proj=proj, steps=2, sigma = 0.50, savePlot = False)