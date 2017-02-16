from optima_tb.utils import odict, OptimaException

import logging
logger = logging.getLogger(__name__)

from uuid import uuid4 as uuid

class ProgramSet:
    
    def __init__(self, name='default'):
        self.name = name
        self.uid = uuid()
        
        self.progs = list()
        self.prog_ids = dict()
        
        logging.info("Created ProgramSet: %s"%self.name)
        
    def makeProgs(self, data):
        for l, prog_label in enumerate(data['progs']):
            prog_name = data['progs'][prog_label]['name']
            prog_type = data['progs'][prog_label]['prog_type']
            t = data['progs'][prog_label]['t']
            cost = data['progs'][prog_label]['cost']
            cov = data['progs'][prog_label]['cov']
            new_prog = Program(name = prog_name, label = prog_label, prog_type = prog_type, t = t, cost = cost, cov = cov)
            self.progs.append(new_prog)
            self.prog_ids[prog_label] = l 
        
    def getProg(self, label):
        if label in self.prog_ids.keys():
            return self.progs[self.prog_ids[label]]
        raise OptimaException('ERROR: Label "%s" cannot be found in program set "%s".' % (label, self.name))
        
        

class Program:
    
    def __init__(self, name, label, prog_type, t = None, cost = None, cov = None):#, duration, category=None):
        """
        
        """
        self.name = name 
        self.label = label
        self.prog_type = prog_type
        self.uid = uuid()
        
        if t is None: t = []
        if cost is None: cost = []
        if cov is None: cov = []
        self.t = t                      # Time data.
        self.cost = cost                # Spending data.
        self.cov = cov                  # Coverage data.
        
#        self.duration = duration
#        self.category = category
        
#        self.cost_coverage = []
        
        
#class TreatmentProgram(Program):
#    
#    def __init__(self,efficacy,adherence,*args):
#        
#        super(TreatmentProgram,self).init(*args)
#        self.efficacy = efficacy
#        self.adherence = adherence
#        
#class TestingProgram(Program):
#    
#    def __init__(self,specificity,sensitivity,*args):
#        
#        super(TreatmentProgram,self).init(*args)
#        self.specificity = specificity
#        self.sensitivity = sensitivity
        
