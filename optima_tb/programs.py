from uuid import uuid4 as uuid

class ProgramSet:
    
    def __init__(self):
        
        self.programs = odict()
        
        

class Program:
    
    def __init__(self,name,label,duration,category=None):
        """
        
        """
        self.name = name 
        self.label = label
        self.uid = uuid()
        
        self.duration = duration
        self.category = category
        
        self.cost_coverage = []
        
        
class TreatmentProgram(Program):
    
    def __init__(self,efficacy,adherence,*args):
        
        super(TreatmentProgram,self).init(*args)
        self.efficacy = efficacy
        self.adherence = adherence
        
class TestingProgram(Program):
    
    def __init__(self,specificity,sensitivity,*args):
        
        super(TreatmentProgram,self).init(*args)
        self.specificity = specificity
        self.sensitivity = sensitivity
        
