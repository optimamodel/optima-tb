# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 11:32:07 2017

@author: Azfar
"""

from uuid import uuid4 as uuid
from utils import odict

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
        
        self.cost_cov = CostCovOut()
    
    def getCoverage():
        '''Returns coverage for a time/spending array'''
        pass
    
    def getBudget():
        '''Returns budget for a time/coverage array'''
        pass
    
    def plotCostCov():
        ''' Plot the cost-coverage curve for a single program'''
        pass
        
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

class CostCovOut(object):
    '''
    Cost-coverage, coverage-outcome and cost-outcome objects
    '''
    def __init__(self):
        self.cost_data = {'t': [], 'spending': [], 'coverage': []}
        self.ccopar = {'t': [], 'saturation': [], 'unitcost': []}
    
    def addcostdata(self, cost_data, overwrite=False):
        '''add or replace spending-coverage data'''
        pass
    
    def rmcostdata(self, t):
        '''remove spending-coverage data against a given year'''
        pass
    
    def getcostdata(self, t):
        '''retrieve/return odict of cost-coverage data'''
        pass
    
    def addccopar(self, ccopars, overwrite=False):
        '''add or replace cost-coverage parameters'''
        pass
    
    def rmccopar(self, t):
        '''remove cost-coverage parameters'''
        pass
    
    def getccopar(self, t):
        '''retrieve/return odict of cost-coverage parameters'''
        pass
    
    def costcov(budget, costCovInfo, popsize, eps=None):
        '''Function that calculates and returns coverage in a given year for a given spending amount.'''
        pass
        
    def inversecostcov(current_coverage, costCovInfo, popsize, eps=None):
        '''Function that calculates and returns budget in a given year for a given coverage amount.'''
        pass