import sys
try:
    sys.path.remove('d:\\work projects\\optima\\optima 2.0')
    sys.path.remove('d:\\work projects\\optima\\optima 2.0\\optima')
except: pass
sys.path.append('../optima')

from optima_tb.project import Project
import pylab

proj= Project(name = 'Belarus', cascade_path = './cascade-belarus.xlsx')
#proj.makeSpreadsheet(databook_path = './databook_belarus_template.xlsx', num_pops = 6, num_migrations = 2)

#set the year range we simulate over as starting in 1995:
#proj.setYear([1995,2030],False)



proj.loadSpreadsheet(databook_path = './databook-belarus.xlsx')
proj.makeParset()

# Run the autofit calibration process. ----------------------------
# The default amount of runs is 500 iterations. If you find it needs to run for more iterations, then 
# uncomment out the next line:
#proj.settings.autofit_params['MaxIter'] = 1000
# This line runs the autofit calibration process, and fits to only two characteristics ('num_vac', 'ac'_inf'):
proj.runAutofitCalibration(new_parset_name='autofit',target_characs=['num_vac','ac_inf']) 
# End of the autofit calibration process. -------------------------

# Run our simulation using the new process
results2 = proj.runSim(parset_name='autofit',plot=True)


pylab.show()