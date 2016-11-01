#%% Imports.

from pyparsing import Word, Literal, Optional, alphanums, nums, ZeroOrMore, Group, Forward
import operator



#%% Parser for functions written as strings

class functionParser(object):
    '''
    A parser that can decompose a string of numbers and variables into tokens, then evaluate the expression.
    Incorporates BODMAS operations, as well as exponentiation and negative numbers.
    Any modifications to this should be made with care, given the left-to-right parsing grammar.
    
    Script written by Paul McGuire ~2006 used as a basis.
    Accordingly, implementing further tokens (e.g. exponential functions) should search fourFn.py online.
    
    Version: 2016nov01
    '''

    def pushFirst(self, s, loc, toks):
        self.expr_stack.append(toks[0])
        
    def pushUnaryMinus(self, s, loc, toks):
        if toks and toks[0]=='-': 
            self.expr_stack.append('u-')

    def __init__(self, debug = False):
        self.debug = debug        
        
        plus    = Literal("+")
        minus   = Literal("-")
        mult    = Literal("*")
        div     = Literal("/")
        neg     = Optional("-")
        addop  = plus | minus
        multop = mult | div
        expop = Literal("^")
        lpar  = Literal('(').suppress()
        rpar  = Literal(')').suppress()
        num = Word(nums + ".")
        var = Word(alphanums + "_")
        
        self.op_dict = {"+": operator.add,
                        "-": operator.sub,
                        "*": operator.mul,
                        "/": operator.truediv,
                        "^": operator.pow}

        self.blah = {"spd_infxness": 2,
                     "spd_infx": 1.5,
                     "alive": 8.2}
        
        self.grammar = Forward()
        primary = (neg + ((num | var).setParseAction(self.pushFirst) | Group(lpar + self.grammar + rpar))).setParseAction(self.pushUnaryMinus)        
        factor = Forward()      # Making sure chain-exponentiation tokens are evaluated from right to left.
        factor << primary + ZeroOrMore((expop + factor).setParseAction(self.pushFirst))
        term = factor + ZeroOrMore((multop + factor).setParseAction(self.pushFirst))
        self.grammar << term + ZeroOrMore((addop + term).setParseAction(self.pushFirst))
        
    def evaluateStack(self, stack, level = None):
        if level == None:
            if self.debug: print('Progressing through stack evaluation...')
            level = 0
        op = stack.pop()
        if op == 'u-':
            return -self.evaluateStack(stack, level = level + 1)
        elif op in "+-*/^":
            op2 = self.evaluateStack(stack, level = level + 1)
            op1 = self.evaluateStack(stack, level = level + 1)
            if self.debug: print('Level %i: %s %s %s = %s' % (level, op1, op, op2, self.op_dict[op](op1, op2)))
            return self.op_dict[op](op1, op2)
        elif op[0].isalpha():
            return self.blah[op]
        else:
            return float(op)
            
    def parse(self, string):
        self.expr_stack = []
        val = self.grammar.parseString(test)
        if self.debug:
            print('Token decomposition...')
            print val
            print('Token stack for evaluation...')
            print self.expr_stack
        return self.evaluateStack(self.expr_stack)

test = "-6+(-spd_infxness*spd_infx+(-3/-4.0)^(-2+6)^0.5)/(alive-6.2)"

fp = functionParser(debug = True)
print fp.parse(test)
