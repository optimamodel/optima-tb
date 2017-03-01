from optima_tb.utils import OptimaException
from optima_tb.interpolation import interpolateFunc

import logging
logger = logging.getLogger(__name__)

from copy import deepcopy as dcp
import numpy as np
from uuid import uuid4 as uuid

class ProgramSet:
    
    def __init__(self, name='default'):
        self.name = name
        self.uid = uuid()
        
        self.progs = list()
        self.prog_ids = dict()
        
        self.impacts = dict()
        
        logging.info("Created ProgramSet: %s"%self.name)
        
    def makeProgs(self, data, settings):
        for l, prog_label in enumerate(data['progs']):
            prog_name = data['progs'][prog_label]['name']
            prog_type = data['progs'][prog_label]['prog_type']
            t = data['progs'][prog_label]['t']
            cost = data['progs'][prog_label]['cost']
            cov = data['progs'][prog_label]['cov']
            cost_format = data['progs'][prog_label]['cost_format']
            cov_format = data['progs'][prog_label]['cov_format']
            target_pops = data['progs'][prog_label]['target_pops']
            prog_type_name = data['progs'][prog_label]['prog_type']
            prog_type_label = settings.progtype_name_labels[prog_type_name]
            target_pars = settings.progtype_specs[prog_type_label]['impact_pars']
            for target_par in target_pars:
                if target_par not in self.impacts: self.impacts[target_par] = []
                self.impacts[target_par].append(prog_label)
            new_prog = Program(name = prog_name, label = prog_label, prog_type = prog_type,
                               t = t, cost = cost, cov = cov, 
                               cost_format = cost_format, cov_format = cov_format,
                               target_pops = target_pops, target_pars = target_pars)
            
            func_pars = dict()
            func_pars['unit_cost'] = data['progs'][prog_label]['unit_cost']
            new_prog.genFunctionSpecs(func_pars = func_pars)
            
            self.progs.append(new_prog)
            self.prog_ids[prog_label] = l
        
    def getProg(self, label):
        if label in self.prog_ids.keys():
            return self.progs[self.prog_ids[label]]
        raise OptimaException('ERROR: Label "%s" cannot be found in program set "%s".' % (label, self.name))
        
        

class Program:
    
    def __init__(self, name, label, prog_type, t = None, cost = None, cov = None, cost_format = None, cov_format = None, target_pops = None, target_pars = None):
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
        self.cost_format = cost_format
        self.cov_format = cov_format
        
        if target_pops is None: target_pops = []
        self.target_pops = target_pops
        
        if target_pars is None: target_pars = []
        self.target_pars = target_pars
        
        self.func_specs = dict()
        
    def interpolate(self, tvec = None, attribute = None):
        ''' Takes attribute values and constructs a dictionary of interpolated array corresponding to the input time vector. Ignores np.nan. '''
        
        # Validate input.
        if tvec is None: raise OptimaException('ERROR: Cannot interpolate attributes of program "%s" without providing a time vector.' % self.label)
        if not len(self.t) > 0: raise OptimaException('ERROR: There are no timepoint values for program "%s".' % self.label)
        if not len(self.t) == len(self.cost): raise OptimaException('ERROR: Program "%s" does not have corresponding time and cost values.' % self.label)
        if not len(self.t) == len(self.cov): raise OptimaException('ERROR: Program "%s" does not have corresponding time and coverage values.' % self.label)

        output = dict()

        input_cost = dcp(self.cost)
        input_cov = dcp(self.cov)
        val_types = ['cost','cov']
        val_arrays = [input_cost,input_cov]
        
        for k in xrange(len(val_types)):
            val_type = val_types[k]
            val_array = val_arrays[k]
            if attribute is not None and attribute != val_type: continue
        
            if len(self.t) == 1:    # If there is only one timepoint, corresponding cost and cov values should be real valued after loading databook. But can double-validate later.
                output[val_type] = np.ones(len(tvec))*(val_array)[0]   # Don't bother running interpolation loops if constant. Good for performance.
            else:
                t_array = dcp(self.t)   # Need to refresh this during each attribute interpolation loop.
                
                # Eliminate np.nan from value array before interpolation. Makes sure timepoints are appropriately constrained.
                t_temp = dcp(t_array)
                val_temp = dcp(val_array)
                t_array = dcp(t_temp[~np.isnan(val_temp)])
                val_array = dcp(val_temp[~np.isnan(val_temp)])
                             
                # Pad the input vectors for interpolation with minimum and maximum timepoint values, to avoid extrapolated values blowing up.
                ind_min, t_min = min(enumerate(t_array), key = lambda p: p[1])
                ind_max, t_max = max(enumerate(t_array), key = lambda p: p[1])
                val_at_t_min = val_array[ind_min]
                val_at_t_max = val_array[ind_max]
                
                # This padding effectively keeps edge values constant for desired time ranges larger than data-provided time ranges.
                if tvec[0] < t_min:
                    t_array = np.append(tvec[0], t_array)
                    val_array = np.append(val_at_t_min, val_array)
                if tvec[-1] > t_max:
                    t_array = np.append(t_array, tvec[-1])
                    val_array = np.append(val_array, val_at_t_max)
                
                output[val_type] = interpolateFunc(t_array, val_array, tvec)
        
        return output
        
    def genFunctionSpecs(self, func_pars, func_type = 'linear'):
        
        self.func_specs['type'] = func_type
        
#        # WARNING: HARD-CODED AND SIMPLISTIC UNIT-COST GENERATION METHOD. IMAGINE IF THERE IS ZERO SPENDING FOR A PROGRAM IN THE LAST YEAR. AMEND ASAP.
#        output = self.interpolate(tvec=[max(self.t)])    # Use the latest year stored in the program to inform unit costs.
        self.func_specs['pars'] = dict()
#        self.func_specs['pars']['unit_cost'] = output['cov'][-1]/output['cost'][-1]
        self.func_specs['pars']['unit_cost'] = func_pars['unit_cost']
        
    def getDefaultBudget(self, year = None):
        '''
        Returns program cost, interpolated either for a given year or the last year that data exists for, as a budget.
        Note that users may enter a unit cost into assumptions, which will provide a flawed budget value in the absence of other data.
        '''
        
        if year is None: year = max(self.t)
        output = self.interpolate(tvec=[year], attribute='cost')
        budget = output['cost'][-1]
        return budget
        
    def getImpact(self, budget):
        
        # WARNING: ASSUMING COVERAGE IS IMPACT.
        if self.cov_format.lower() == 'fraction':
            return self.func_specs['pars']['unit_cost']*budget*0.01     # Unit cost is per percentage when format is a fraction.
        else:
            return self.func_specs['pars']['unit_cost']*budget
        

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
        
