# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 11:32:07 2017

@author: Azfar
"""

from uuid import uuid4 as uuid
from optima_tb.utils import odict

class ProgramSet:
    
    def __init__(self):
        
        self.programs = odict()
        
        

class Program:
    def __init__(self,name,label,category=None):
        """
        
        """
        self.name = name 
        self.label = label
        self.uid = uuid()
        
        self.category = category
    
        self.cost_cov = CostCovOut()
    
    
    
    def getCoverage(self,budget):
        '''Returns coverage for a time/budget array'''
        pass
    
    def getBudget(self,coverage):
        '''Returns budget for a time/coverage array'''
        pass
    
    def getImpacts(self,budget):
        ''' Returns impacts for a time/budget array'''
        pass
    
    def plotCostCov(self):
        ''' Plot the cost-coverage curve for a single program'''
        # should pass data to plotting.py
        pass
        
class TreatmentProgram(Program):
    
    def __init__(self,efficacy,adherence,duration,*args):
        
        super(TreatmentProgram,self).init(*args)
        self.efficacy = efficacy
        self.adherence = adherence
        
        self.duration = duration
        
class TestingProgram(Program):
    
    def __init__(self,specificity,sensitivity,frequency,*args):
        
        super(TreatmentProgram,self).init(*args)
        self.specificity = specificity
        self.sensitivity = sensitivity
        self.frequency = frequency

class CostCovOut:
    '''
    Cost-coverage, coverage-outcome and cost-outcome objects
    '''
    def __init__(self):
        self.cost_data = {'t': [], 'spending': [], 'coverage': []}
        self.cost_cov_param = {'t': [], 'saturation': [], 'unitcost': []}
        self.curve = CostCoverageCurve()
    
    def addCostData(self, cost_data, overwrite=False):
        '''add or replace spending-coverage data'''
        pass
    
    def removeCostData(self, t):
        '''remove spending-coverage data against a given year'''
        pass
    
    def getCostData(self, t):
        '''retrieve/return odict of cost-coverage data'''
        pass
    
    def addCostCoverageParam(self, ccopars, overwrite=False):
        '''add or replace cost-coverage parameters'''
        pass
    
    def removeCostCoverageParam(self, t):
        '''remove cost-coverage parameters'''
        pass
    
    def getCostCoverageParam(self, t):
        '''retrieve/return odict of cost-coverage parameters'''
        pass
    
    def calcCostCoverage(self, budget, costCovInfo, popsize, eps=None):
        '''Function that calculates and returns coverage in a given year for a given spending amount.'''
        pass
        
    def calcInverseCostCoverage(self, current_coverage, costCovInfo, popsize, eps=None):
        '''Function that calculates and returns budget in a given year for a given coverage amount.'''
        pass
    
    
    
    
class CostCoverageCurve:
    
    import math
    """
    
    Curves currently supported: 
    1) exponentially saturating, from (0,0)
    2) exponentially saturating with offset, i.e. curve only defined for >= (offset, 0); 0, otherwise
    3) sigmoidal
    4) fixed cost
    
    """
    
    def __init__(self,curve_shape="exponential",params={'exponent':1.,'saturation':1e6}):
        
        try:
            self.curve_shape = curve_shape
            self.curve_params = params
            self.curve = getattr(self, '_curve_%s'%curve_shape)(params)
        except:
            raise NotImplementedError("No curve associated with _curve_%s (programs.py)"%curve_shape)
            
    def getValue(self,budget):
        return self.curve(budget)
            
    def getInverseValue(self,current_coverage,**params):
        try:
            inverse_function = getattr(self, '_inverse_curve_%s'%self.curve_shape)(params)
            return inverse_function(current_coverage)
        except:
            raise NotImplementedError("No curve associated with _curve_%s (programs.py)"%curve_shape)
       
    """ Internal representations of curves """
    def _curve_exponential( exponent, saturation):
        return saturation*(1-math.exp(-exponent))
    
    def _inverse_curve_exponential(self,exponent,saturation):
        pass
    
    def _curve_exponential_offset(self, exponent, saturation, offset):
        pass
    
    def _inverse_exponential_offset(self,exponent,saturation, offset):
        pass
    
    def _curve_sigmoidal(self,sigma, saturation):
        pass

    def _inverse_curve_sigmoidal(self,sigma,saturation):
        pass
    
    def _curve_fixed(self,fixed):
        pass
    
    def _inverse_curve_fixed(self,fixed):
        pass
    
 
    