from optima_tb.project import Project
from optima_tb.utils import odict
from optima_tb.analysis import calculateCumulativeDerivatives

import pylab

path = '../tb-ucl-analyses/belarus/'
dt = 1.0/4
plot_over = (2000,2030)
proj = Project(name = 'Belarus', cascade_path = '%s/cascade-belarus.xlsx'%path, validation_level = 'avert', plotting_level = 'presentation')
proj.settings.tvec_dt = dt
#setup_belarusPlotting(proj,dt,plot_over)

#set the year range we simulate over as starting in 1995:
proj.setYear([1999,2030],False)


proj.loadSpreadsheet(databook_path = '%s/databook-belarus.xlsx'%path)
proj.makeParset(name = 'default')
settings = proj.settings


results = proj.runSim()



comp_labels = ['sndd']
comp_titles = ['Smear Negative Notifications']
pop_labels = ['5-14','15-64']
pop_titles = ['5-14','15-64']

link_labels = ['snddiag_rate']
include_link_not_exclude = True
plot_inflows = True
plot_outflows = False
exclude_transfers = True
link_legend = {'snddiag_rate':'Notifications'}
sum_total = True

from_year = 2010
to_year = 2015

snnot_all =     calculateCumulativeDerivatives(results=results, settings=proj.settings, 
                                                from_year = from_year, to_year=to_year, 
                                                comp_labels=comp_labels, comp_titles=comp_titles, 
                                                pop_labels=pop_labels, pop_titles=pop_titles,
                                                link_labels=link_labels, include_link_not_exclude=include_link_not_exclude,
                                                link_legend=link_legend,
                                                plot_inflows=plot_inflows, plot_outflows=plot_outflows,
                                                sum_total = sum_total,
                                                exclude_transfers=exclude_transfers)


