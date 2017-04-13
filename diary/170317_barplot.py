from optima_tb.plotting import _plotBars, plotBudgets
from optima_tb.settings import Settings
from optima_tb.project import Project
from optima_tb.utils import odict
import pylab



values = [ [5,6,4], [4,10, 1] ] #,[10,3, 2]]

plotting_level = 'dev'

proj = Project(cascade_path='data/cascade-simple-calibration.xlsx', plotting_level=plotting_level)
settings = Settings(cascade_path='data/cascade-simple-calibration.xlsx', plotting_level=plotting_level)

#_plotBars(values, **settings.plot_settings)


b1 = odict()
b1['Prog1'] = values[0][0]
b1['Prog2'] = values[0][1]
b1['Prog3'] = values[0][2]
b2 = odict()
b2['Prog1'] = values[1][0]
b2['Prog2'] = values[1][1]
b2['Prog3'] = values[1][2]
budgets = [b1,b2]
b_labels = ['Current conditions', 'Optimized']

plotBudgets(budgets, plotdict=settings.plot_settings, xlabels=b_labels, title='Example budgets',
            save_fig=True, fig_name="example-budget.png")


pylab.show()