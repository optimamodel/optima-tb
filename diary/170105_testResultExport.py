from optima_tb.project import Project

"""
Example to illustrate how to export a set of results for a simple model run

Required: 
    The following files are required and placed at 'Homefolder'\data
        1. databook-simple-cascadeBirths.xlsx
        2. cascade-simple-births.xlsx
          
Generates: sample-output.xlsx, which will be the same format as databook-simple-cascade.xlsx


@date:   05-Jan-2017
@author: hussainazfar

"""

databook = './data/databook-simple-cascadeBirths.xlsx'
cascade =  './data/cascade-simple-births.xlsx'


proj= Project(name = 'test-simple-birth', cascade_path = cascade)
proj.loadSpreadsheet(databook_path = databook)
# make initial parameter set 
proj.makeParset()
#run model on parset
r1 = proj.runSim(plot = False)
#export characteristics (i.e. save results)
r1.export()

