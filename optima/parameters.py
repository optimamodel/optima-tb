#%% Imports





#%% Parset class that contains one set of parameters converted from raw project data

class ParameterSet(object):
    ''' Class to hold all parameters and information on how they were generated, and perform operations on them'''
    
    def __init__(self, name='default'):
        self.name = name 
        self.pop_names = []       # List of population names.
    
    def makePars(self, data):
        self.pop_names = data['pops']['name_labels'].keys()