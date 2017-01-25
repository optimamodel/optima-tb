

from optima_tb.scenarios import ParameterScenario
from optima_tb.utils import odict
from optima_tb.settings import Settings, DO_NOT_SCALE
from optima_tb.project import Project


scvalues = odict()
param = 'birth_transit'
scvalues[param] = odict()
scvalues[param]['Pop1'] = odict()
scvalues[param]['Pop1']['y'] = [3e6, 2e6]
scvalues[param]['Pop1']['t'] = [2003.,2005.]
scvalues[param]['Pop1']['y_format'] = 'number'
scvalues[param]['Pop1']['y_factor'] = DO_NOT_SCALE
param = 'ltfu_a'
scvalues[param] = odict()
scvalues[param]['Pop2'] = odict()
scvalues[param]['Pop2']['y'] = [0.9999, 0.9]
scvalues[param]['Pop2']['t'] = [2000.,2005.]
scvalues[param]['Pop2']['y_format'] = 'fraction'
scvalues[param]['Pop2']['y_factor'] = DO_NOT_SCALE

pops = {'Pop1':'Pop1','Pop2':'Pop2'}

pscenario = ParameterScenario(name="testPS",scenario_values=scvalues,pop_labels=pops)

databook = 'data/databook-simple-cascade-autocalibration.xlsx'
  
proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple-calibration.xlsx')
proj.setYear([2000.,2030.],False) 
proj.loadSpreadsheet(databook_path = databook)
proj.makeParset()
parset = proj.parsets[0]

parset1 = parset + (pscenario.scenario)
print parset1.pars['cascade'][0]
print parset1.pars['cascade'][12]
parset1.name = "charlie"
parset.name = "bob"
print parset1.name # should be charlie


parset2 = parset << pscenario.scenario
print parset2.pars['cascade'][0]
print parset2.pars['cascade'][12]
parset2.name = "eve"
parset.name = "alice"
print parset2.name # should be eve
