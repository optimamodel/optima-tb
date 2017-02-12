
import pylab

from optima_tb.scenarios import ParameterScenario
from optima_tb.utils import odict
from optima_tb.settings import Settings, DO_NOT_SCALE
from optima_tb.project import Project


"""

Example usage for parameter scenario: 

This demo shows the underlying functionality of how to 
create and run a parameter scenario, and illustrates the different 
options

"""


scvalues = odict()
param = 'birth_transit'
scvalues[param] = odict()
scvalues[param]['Pop1'] = odict()
scvalues[param]['Pop1']['y'] = [3e6, 1e4, 1e4, 2e6]
scvalues[param]['Pop1']['t'] = [2003.,2004.,2014.,2015.]
scvalues[param]['Pop1']['y_format'] = 'number'
scvalues[param]['Pop1']['y_factor'] = DO_NOT_SCALE
param = 'ltfu_a'
scvalues[param] = odict()
scvalues[param]['Pop2'] = odict()
scvalues[param]['Pop2']['y'] = [0.999, 0.0, 0.9]
scvalues[param]['Pop2']['t'] = [2000.,2010.,2025.]
scvalues[param]['Pop2']['y_format'] = 'fraction'
scvalues[param]['Pop2']['y_factor'] = DO_NOT_SCALE

pops = {'Pop1':'Pop1','Pop2':'Pop2'}


"""

Creating and running parameter scenarios directly.

Note that it's not advised to do this in practice, but this example
is included 

"""
def testDirect():
    """
    
    
    """
    pscenario = ParameterScenario(name="testPS",scenario_values=scvalues,pop_labels=pops)
    
    databook = 'data/databook-simple-cascade-autocalibration.xlsx'
      
    proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple-calibration.xlsx')
    proj.setYear([2000.,2030.],False) 
    proj.loadSpreadsheet(databook_path = databook)
    proj.makeParset()
    parset = proj.parsets[0]
    
    """ Adds the values of pscenario to those of the first parset """
    parset1 = parset + (pscenario.scenario_parset)
    print parset1.pars['cascade'][0] # should be 301000 for t=2003
    print parset1.pars['cascade'][12] # should be capped at 1 for t=2000
    parset1.name = "charlie"
    parset.name = "bob"
    print parset1.name # should be charlie
    
    """ Overwrites overlapping values of the first parset with those of pscenario """
    parset2 = parset << pscenario.scenario_parset
    print parset2.pars['cascade'][0] # should be 3e6 for t=2003
    print parset2.pars['cascade'][12] # should be 0.999 for t=2000
    parset2.name = "eve"
    parset.name = "alice"
    print parset2.name # should be eve

#testDirect()



"""

Using the same functionality via project interface:

See the docstring documentation in project.py / createScenarios

"""
def testProject():
    databook = 'data/databook-simple-cascade-autocalibration.xlsx'
    proj= Project(name = 'test-Belarus-simple', cascade_path = 'data/cascade-simple-calibration.xlsx')
    proj.setYear([2000.,2030.],False) 
    proj.loadSpreadsheet(databook_path = databook)
    proj.makeParset()
    parset = proj.parsets[0]
    scen_values = { 'test_scenario': {'type': 'Parameter',
                                      'overwrite' : True, # it will overwrite scenario to the parset
                                      'run_scenario' : True,
                                      'scenario_values': scvalues}
                   }
    proj.createScenarios(scen_values)
    resultset = proj.runScenarios(parset.name,include_bau=True)
    for resultName in resultset.keys():
        proj.plotResults(resultset[resultName],savePlot=True,figName=resultName)
        pylab.show()
        
    return resultset, proj

results, proj = testProject()



