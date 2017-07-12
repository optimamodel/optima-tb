from optima_tb.utils import OptimaException, odict
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

        logging.info("Created ProgramSet: %s" % self.name)

    def makeProgs(self, data, settings):
        for l, prog_label in enumerate(data['progs']):
            prog_name = data['progs'][prog_label]['name']
            prog_type = data['progs'][prog_label]['prog_type']

            special = None
            if 'special' in settings.progtype_specs[prog_type]:
                special = settings.progtype_specs[prog_type]['special']

            t = data['progs'][prog_label]['t']
            cost = data['progs'][prog_label]['cost']
            cost_format = data['progs'][prog_label]['cost_format']
            cov = data['progs'][prog_label]['cov']
            cov_format = data['progs'][prog_label]['cov_format']
            attributes = data['progs'][prog_label]['attributes']
            target_pops = data['progs'][prog_label]['target_pops']
            target_pars = settings.progtype_specs[prog_type]['impact_pars']
            for target_par in target_pars.keys():
                if target_par not in self.impacts: self.impacts[target_par] = []
                self.impacts[target_par].append(prog_label)

            new_prog = Program(name=prog_name, label=prog_label, prog_type=prog_type,
                               t=t, cost=cost, cov=cov,
                               cost_format=cost_format, cov_format=cov_format,
                               attributes=attributes,
                               target_pops=target_pops, target_pars=target_pars)


            func_pars = dict()
            if special == 'cost_only':
                func_type = 'cost_only'
            else:
                func_type = 'linear'
                func_pars['unit_cost'] = data['progs'][prog_label]['unit_cost']
            new_prog.genFunctionSpecs(func_pars=func_pars, func_type=func_type)

            self.progs.append(new_prog)
            self.prog_ids[prog_label] = l

    def getProg(self, label):
        if label in self.prog_ids.keys():
            return self.progs[self.prog_ids[label]]
        raise OptimaException('ERROR: Label "%s" cannot be found in program set "%s".' % (label, self.name))

    def getBudgets(self, coverages):
        """
        Returns the budgets for each of the programs in this program set
        
        Params:
            coverages   dict of prog.label,coverage pairs
            
        Returns
                        dict of prog.label, budget pairs
        """
        budgets = {}
        for prog_label in coverages.keys():
            prog = self.progs[self.prog_ids[prog_label]]
            budgets[prog.label] = prog.getBudget(coverages[prog_label])
        return budgets

    def getCoverages(self, budgets):
        """
        Returns the coverage for each of the programs in this program set
        
        Params:
            budgets     dict of prog.label, budget pairs
            
        Returns
                        dict of prog.label, coverage pairs
        """
        coverage = {}
        for prog_label in budgets.keys():
            prog = self.progs[self.prog_ids[prog_label]]
            coverage[prog.label] = prog.getCoverage(budgets[prog_label])
        return coverage

    def copy(self):
        pass


class Program:

    def __init__(self, name, label, prog_type, t=None, cost=None, cov=None, cost_format=None, cov_format=None, attributes={}, target_pops=None, target_pars=None):
        """
        
        """
        self.name = name
        self.label = label
        self.prog_type = prog_type
        self.uid = uuid()

        if t is None: t = []
        if cost is None: cost = []
        if cov is None: cov = []
        if attributes is None: attributes = odict()
        self.t = t  # Time data.
        self.cost = cost  # Spending data.
        self.cov = cov  # Coverage data.
        self.cost_format = cost_format
        self.cov_format = cov_format
        self.attributes = attributes

        if target_pops is None: target_pops = []
        self.target_pops = target_pops

        if target_pars is None: target_pars = []
        self.target_pars = target_pars

        self.func_specs = dict()

    def interpolate(self, tvec=None, attributes=None):
        ''' Takes attribute values and constructs a dictionary of interpolated array corresponding to the input time vector. Ignores np.nan. '''

        # Validate input.
        if tvec is None: raise OptimaException('ERROR: Cannot interpolate attributes of program "%s" without providing a time vector.' % self.label)
        if not len(self.t) > 0: raise OptimaException('ERROR: There are no timepoint values for program "%s".' % self.label)
        if not len(self.t) == len(self.cost): raise OptimaException('ERROR: Program "%s" does not have corresponding time and cost values.' % self.label)
        if not len(self.t) == len(self.cov): raise OptimaException('ERROR: Program "%s" does not have corresponding time and coverage values.' % self.label)

        output = dict()

        input_cost = dcp(self.cost)
        input_cov = dcp(self.cov)
        val_types = ['cost', 'cov'] + [x for x in self.attributes.keys()]
        val_arrays = [input_cost, input_cov] + [y for y in self.attributes.values()]

        for k in xrange(len(val_types)):
            val_type = val_types[k]
            val_array = val_arrays[k]
            if attributes is not None and val_type not in attributes: continue

            t_array = dcp(self.t)  # Need to refresh this during each attribute interpolation loop.

            # Eliminate np.nan from value array before interpolation. Makes sure timepoints are appropriately constrained.
            t_temp = dcp(t_array)
            val_temp = dcp(val_array)
            t_array = dcp(t_temp[~np.isnan(val_temp)])
            val_array = dcp(val_temp[~np.isnan(val_temp)])

            if len(t_array) == 1:  # If there is only one timepoint, corresponding cost and cov values should be real valued after loading databook. But can double-validate later.
                output[val_type] = np.ones(len(tvec)) * (val_array)[0]  # Don't bother running interpolation loops if constant. Good for performance.
            else:
                # Pad the input vectors for interpolation with minimum and maximum timepoint values, to avoid extrapolated values blowing up.
                ind_min, t_min = min(enumerate(t_array), key=lambda p: p[1])
                ind_max, t_max = max(enumerate(t_array), key=lambda p: p[1])
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

    def genFunctionSpecs(self, func_pars, func_type='linear'):

        self.func_specs['type'] = func_type

#        # WARNING: HARD-CODED AND SIMPLISTIC UNIT-COST GENERATION METHOD. IMAGINE IF THERE IS ZERO SPENDING FOR A PROGRAM IN THE LAST YEAR. AMEND ASAP.
#        output = self.interpolate(tvec=[max(self.t)])    # Use the latest year stored in the program to inform unit costs.
        self.func_specs['pars'] = dcp(func_pars)
#        self.func_specs['pars']['unit_cost'] = output['cov'][-1]/output['cost'][-1]
#        self.func_specs['pars']['unit_cost'] = func_pars['unit_cost']

    def getDefaultBudget(self, year=None):
        '''
        Returns program cost, interpolated either for a given year or the last year that data exists for, as a budget.
        Note that users may enter a unit cost into assumptions, which will provide a flawed budget value in the absence of other data.
        '''

        if year is None: year = max(self.t)
        output = self.interpolate(tvec=[year], attributes=['cost'])
        budget = output['cost'][-1]
        return budget

    def getBudget(self, coverage):
        '''
        Returns a budget that corresponds to a program coverage. In simplest form, this is coverage times unit cost.
        Note that this method is not typically used during model processing.
        '''

        # If coverage is a scalar, make it a float. If it is a list or array, make it an array of floats.
        try: coverage = float(coverage)
        except: coverage = dcp(np.array(coverage, 'float'))

        if self.func_specs['type'] == 'cost_only':
            bud = np.nan    # A fixed-cost program has no coverage, so a coverage-budget conversion should return a NaN.
                            # This will not bother value-updating in the model, as fixed-cost programs should not have impact parameters anyway.
        elif self.cov_format is None:
            raise OptimaException('ERROR: Attempted to convert coverage to budget for a program (%s) that does not have coverage.' % self.label)
        else:
            if self.cov_format.lower() == 'fraction':
                bud = coverage * self.func_specs['pars']['unit_cost'] / 0.01     # Unit cost is per percentage when format is a fraction.
            else:
                bud = coverage * self.func_specs['pars']['unit_cost']
        return bud

    def getCoverage(self, budget):
        '''
        Returns prospective coverage for a program. In simplest form, this is budget divided by unit cost.
        Note that this coverage will be dynamically divided amongst populations/compartments in the model prior to scaling into an impact.
        Excess coverage will also be ignored in the model.
        '''

        # If budget is a scalar, make it a float. If it is a list or array, make it an array of floats.
        try: budget = float(budget)
        except: budget = dcp(np.array(budget, 'float'))

        if self.cov_format is None:
            cov = np.nan  # Fixed cost programs have no coverage.
        else:
            if self.cov_format.lower() == 'fraction':
                cov = budget * 0.01 / self.func_specs['pars']['unit_cost']     # Unit cost is per percentage when format is a fraction.
            else:
                cov = budget / self.func_specs['pars']['unit_cost']
        return cov

    def getImpact(self, budget, impact_label=None, parser=None, years=None, budget_is_coverage=False):

        if self.func_specs['type'] == 'cost_only':
            return 0.0

        budget = dcp(budget)    # Just in case.

        # Baseline impact is just coverage.
        if budget_is_coverage:
            imp = budget
        else:
            imp = self.getCoverage(budget)

        # If impact parameter has an impact function, this is what coverage is scaled by.
        if not impact_label is None:
            if 'f_stack' in self.target_pars[impact_label].keys():
                if parser is None:
                    raise OptimaException('ERROR: Cannot calculate "%s" impact for "%s" without a parser, due to the existence of an impact function.' % (self.label, impact_label))
                f_stack = dcp(self.target_pars[impact_label]['f_stack'])
                attribs = dcp(self.target_pars[impact_label]['attribs'])
                for attrib_label in attribs.keys():
                    if attrib_label in self.attributes.keys():
                        if years is None: years = [max(self.t)]
                        output = self.interpolate(tvec=years, attributes=[attrib_label])
                    else:
                        raise OptimaException('ERROR: Cannot calculate "%s" impact for "%s" due to a missing reference "%s".' % (self.label, impact_label, attrib_label))
                    attribs[attrib_label] = output[attrib_label]# [-1]
#                    if len(output[attrib_label]) > 1:
#                        print attribs
                new_val = parser.evaluateStack(stack=f_stack, deps=attribs)
                imp *= new_val      # Scale coverage.

        return imp
