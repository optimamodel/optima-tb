from optima_tb.utils import OptimaException, odict
from optima_tb.interpolation import interpolateFunc
from functools import partial

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
            flag = settings.progtype_specs[prog_type]['special']
            sat = data['progs'][prog_label]['sat']
            for target_par in target_pars.keys():
                if target_par not in self.impacts: self.impacts[target_par] = []
                self.impacts[target_par].append(prog_label)

            # TODO remove hard-coded tag and use ModelPopulation.getAlive() instead
            # sum over all population sizes associated with the program
            target_pop_size = np.sum([data['characs']['h_alive'][pop]['y'][0] for pop in target_pops])

            new_prog = Program(name=prog_name, label=prog_label, prog_type=prog_type, t=t, cost=cost, cov=cov,
                               cost_format=cost_format, cov_format=cov_format, attributes=attributes,
                               target_pops=target_pops, target_pop_size=target_pop_size, target_pars=target_pars,
                               flag=flag)

            # func_pars = dict()
            # TODO: This is where new function types can be specified, such as sigmoid cost-curves.
            # The func_type and func_pars dict will need to reflect this, e.g. with func_pars['unit_cost'] and func_pars['asymptote'].
            if special == 'cost_only':
                cost_only = True
                # required for defaults.py which requires func_specs['type'] to be set
                new_prog.func_specs['type'] = 'cost_only'
            else:
                cost_only = False
                unit_cost = data['progs'][prog_label]['unit_cost']
                if None in sat: # if there is one None in the sat-list, then all of the entries are None
                    sat = None
                # required for defaults.py which requires func_specs['type'] to be set
                new_prog.func_specs['type'] = 'other'
            new_prog.genFunctionSpecs(unit_cost, sat, cost_only, settings)

            self.progs.append(new_prog)
            self.prog_ids[prog_label] = l

        # set global and second-order programs
        self._setSpecialPrograms(gp_flags=['scale_prop'], sop_flags=['supp'])
        self._findDependentLinks(settings)


    def getProg(self, label):
        if label in self.prog_ids.keys():
            return self.progs[self.prog_ids[label]]
        raise OptimaException('ERROR: Label "%s" cannot be found in program set "%s".' % (label, self.name))


    def getBudgets(self, coverages, model):
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


    def getCoverages(self, budgets, model):
        """
        Returns the coverage for each of the programs in this program set
        
        Params:
            budgets     dict of prog.label, budget pairs
            
        Returns
                        dict of prog.label, coverage pairs
        """
        coverage = {}
        if budgets == None:
            for prog in self.progs:
                coverage[prog.label] = dcp(prog.cov)
        else:
            for prog_label in budgets.keys():
                prog = self.progs[self.prog_ids[prog_label]]
                coverage[prog.label] = prog.getCoverage(budgets[prog_label])
        return coverage


    def copy(self):
        pass


    # extract list of global and second order programs
    def _setSpecialPrograms(self, gp_flags, sop_flags):
        self.GPs = []
        self.GP_pars = []
        self.SOPs = []
        for p in self.progs:
            if p.flag[0] in gp_flags:
                self.GPs.append(p.label)
                self.GP_pars.extend(p.flag[1])
            elif p.flag[0] in sop_flags:
                self.SOPs.append(p.label)
            else:
                pass

        # uniquify
        self.GPs =list(set(self.GPs))
        self.SOPs = list(set(self.SOPs))
        self.GP_pars = list(set(self.GP_pars))


    def _findDependentLinks(self, settings):
        for prog in filter(lambda x: x.flag[0] == 'scale' and x.flag[1], self.progs):
            prog.deps = {}
            for par in prog.flag[1]:
                if not 'tag' in settings.linkpar_specs[par]: # meaning: par is not a link
                    # find function/formula where par is used and trace that back to a link
                    for p in filter(lambda x: 'deps' in settings.linkpar_specs[x] and par in settings.linkpar_specs[x]['deps'], settings.linkpar_specs):
                        if 'tag' in settings.linkpar_specs[p]:
                            # ambiguity check:
                            if not par in prog.deps:
                                prog.deps[par] = settings.linkpar_specs[p]['tag']
                            else:
                                raise OptimaException('Dependency on parameter \'%s\' is not unigue: \'%s\' and \'%s\' require \'%s\'' % (par, self.deps[par], p, par))
                else:
                    prog.deps[par] = settings.linkpar_specs[par]['tag']


class Program:
    def __init__(self, name, label, prog_type, t=None, cost=None, cov=None, cost_format=None, cov_format=None,
                 attributes={}, target_pops=None, target_pop_size=None, target_pars=None, flag='replace'):
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

        self.pop_size = target_pop_size

        if target_pops is None: target_pops = []
        self.target_pops = target_pops

        if target_pars is None: target_pars = []
        self.target_pars = target_pars

        self.func_specs = dict()

        # currently sanity check is disabled
        self.flag = self.parseSpecialTag(flag, None)
        if self.flag[0] == 'supp':
            self._setSOP()


    def _setSOP(self):
        # extract all attributes which refer to another program
        refs = filter(lambda x: x.startswith('$ref_'), self.attributes)
        # list of all attributes other than programs
        var = list(set(self.attributes.keys()).difference(set(refs)))

        self.ref = {}

        for p in refs:
            # suffix of the program, everything after the last '_' in program label (incl. '_')
            suff = p[p.rfind('_'):]
            # find all corresponding label of self.attributes
            tmp = filter(lambda x: x.endswith(suff), var)
            # refine labels: keep only those that reference values which are not all zero
            self.ref[self.attributes[p][0]] = filter(lambda x: not np.any(self.attributes[x]), tmp)


    # expects a 'special tag' of a program. The passed string is split at whitespaces. The first word defines the
    # special behaviour of the program, the other words are put in a list because they are specific for the special
    # behaviour of the program. At the end a sanity check is performed if prog != None
    def parseSpecialTag(self, tag, prog=None):
        # decompose special tag word by word (string separated by whitespaces)
        words = tag.split()

        if prog == None:
            return words[0], words[1:]

        # check if all the provided parameters are indeed impact parameters
        if words[0] != '' and len(words[1:]) > 0:
            for w in words[1:]:  # check if impact parameters and passed words coincide
                if not w in prog.target_pars:
                    raise OptimaException(
                        'ERROR: Parameters passed in the \'special tag\' field after \'%s\' must be impact parameters. \'%s\' is not in list [%s].'
                        % (words[0], w, ' '.join(prog.target_pars.keys())))

        return (words[0], words[1:])


    def insertValuePair(self, t, y, attribute, rescale_after_year=False):
        '''
        Check if the inserted t value already exists for the attribute type.
        If not, append y value. If so, overwrite y value.
        If optional argument rescale_after_year is True, relevant attribute values for timepoints greater than t will be proportionally scaled by y divided by the value it is replacing.
        If the scaling factor is infinite or nan, the later values will just be directly replaced by y.
        '''
        
        val_dict = {'cost':self.cost, 'cov':self.cov}
        for val_type in self.attributes.keys():
            val_dict[val_type] = self.attributes[val_type]
        
        orig_y = self.interpolate(tvec=[t], attributes=[attribute])[attribute][-1]
        try: scaling_factor = float(y)/orig_y
        except: raise OptimaException('ERROR: Unable to insert value "%f" at year "%f" for attribute "%s" of program "%s"' % (y, t, attribute, self.label))
            
        value_replacement = False
        k = 0
        for t_val in self.t:
            if t_val == t:
                value_replacement = True
            if t_val == t or (rescale_after_year is True and t_val > t):
                if np.isnan(scaling_factor) or np.isinf(scaling_factor):
                    new_val = float(y)
                else:
                    new_val = scaling_factor*val_dict[attribute][k]
                val_dict[attribute][k] = new_val
                logging.info('Inserted/replaced value "%f" at year "%f" for attribute "%s" of program "%s"' % (new_val, t_val, attribute, self.label))
            k += 1
        
        if value_replacement: return    # No need to continue and append new values if the target year exists.

        try:
            self.t = np.append(self.t, float(t))
            self.cost = np.append(self.cost, np.nan)
            self.cov = np.append(self.cov, np.nan)
            for val_type in self.attributes.keys():
                self.attributes[val_type] = np.append(self.attributes[val_type], np.nan)

            if attribute == 'cost':
                self.cost[-1] = float(y)
            elif attribute == 'cov':
                self.cov[-1] = float(y)
            else:
                self.attributes[attribute][-1] = float(y)
        except:
            raise OptimaException(
                'ERROR: Unable to insert value "%f" at year "%f" for attribute "%s" of program "%s"' % (
                y, t, attribute, self.label))


    def interpolate(self, tvec=None, attributes=None):
        ''' Takes attribute values and constructs a dictionary of interpolated array corresponding to the input time vector. Ignores np.nan. '''

        # Validate input.
        if tvec is None: raise OptimaException(
            'ERROR: Cannot interpolate attributes of program "%s" without providing a time vector.' % self.label)
        if not len(self.t) > 0: raise OptimaException(
            'ERROR: There are no timepoint values for program "%s".' % self.label)
        if not len(self.t) == len(self.cost): raise OptimaException(
            'ERROR: Program "%s" does not have corresponding time and cost values.' % self.label)
        if not len(self.t) == len(self.cov): raise OptimaException(
            'ERROR: Program "%s" does not have corresponding time and coverage values.' % self.label)

        output = dict()

        input_cost = dcp(self.cost)
        input_cov = dcp(self.cov)
        val_types = ['cost', 'cov'] + [x for x in self.attributes.keys()]
        val_arrays = [input_cost, input_cov] + [y for y in self.attributes.values()]

        for k in xrange(len(val_types)):
            val_type = val_types[k]
            val_array = val_arrays[k]
            if attributes is not None and val_type not in attributes: continue
            # in this case val_array contains strings which raise an exception on isnan()
            if val_type.startswith('$ref'): continue

            # Need to refresh this during each attribute interpolation loop.
            t_array = dcp(self.t)

            # Eliminate np.nan from value array before interpolation. Makes sure timepoints are appropriately
            # constrained.
            t_temp = dcp(t_array)
            val_temp = dcp(val_array)
            t_array = dcp(t_temp[~np.isnan(val_temp)])
            val_array = dcp(val_temp[~np.isnan(val_temp)])

            # If there is only one timepoint, corresponding cost and cov values should be real valued after loading
            # databook. But can double-validate later.
            if len(t_array) == 1:
                # Don't bother running interpolation loops if constant. Good for performance.
                output[val_type] = np.ones(len(tvec)) * (val_array)[0]
            else:
                # Pad the input vectors for interpolation with minimum and maximum timepoint values, to avoid
                # extrapolated values blowing up.
                ind_min, t_min = min(enumerate(t_array), key=lambda p: p[1])
                ind_max, t_max = max(enumerate(t_array), key=lambda p: p[1])
                val_at_t_min = val_array[ind_min]
                val_at_t_max = val_array[ind_max]

                # This padding effectively keeps edge values constant for desired time ranges larger than data-provided
                # time ranges.
                if tvec[0] < t_min:
                    t_array = np.append(tvec[0], t_array)
                    val_array = np.append(val_at_t_min, val_array)
                if tvec[-1] > t_max:
                    t_array = np.append(t_array, tvec[-1])
                    val_array = np.append(val_array, val_at_t_max)

                output[val_type] = interpolateFunc(t_array, val_array, tvec)

        return output

    def genFunctionSpecs(self, unit_cost, sat=None, cost_only=False, settings=None):
        # getCoverage/getBudget requires 2 parameters, but not all cost curve types need 2 parameters. So lambda
        # function takes to arguments so that getCoverage/getBudget can be uniformly called and simply discards
        # unnecessary parameters
        if sat is None and cost_only:
            # cost only case
            self.__bud_func = constCostCov
            self.__cov_func = constCostCov

        # Always return the same type (fraction/number) for coverage as the given input!

        elif sat is not None and not cost_only:
            # sigmoidal case
            # raises exception if sat is not iterable; but that is what it is meant to be (a single float), so don't do
            # anything in exception handling
            try: sat = sat[0]
            except: pass

            frac = (self.cov_format.lower() == 'fraction')

            # sat is assumed to be a fraction
            self.__bud_func = partial(logCovCost, saturation=sat, unit_cost=unit_cost,
                                      pop_size=self.pop_size, is_fraction=frac, dt_year=settings.tvec_dt)
            self.__cov_func = partial(logCostCov, saturation=sat, unit_cost=unit_cost,
                                      pop_size=self.pop_size, is_fraction=frac, dt_year=settings.tvec_dt)

            # Use truncated linear functions as fallback plan, BUT DO get rid of the lambda expressions or parallelism
            # will fail
            # if self.cov_format.lower() == 'fraction':
            #     self.__bud_func = lambda cov, pop: np.minimum(cov, sat) * unit_cost / 0.01
            #     self.__cov_func = lambda bud, pop: np.minimum(bud * 0.01 / unit_cost, sat)
            # else:
            #     self.__bud_func = lambda cov, pop: np.minimum(cov * unit_cost, sat * pop)
            #     self.__cov_func = lambda bud, pop: np.minimum(bud / unit_cost, sat * pop)

        else:
            # linear case
            frac = (self.cov_format.lower() == 'fraction')
            # Unit cost is per percentage when format is a fraction.
            self.__bud_func = partial(linCovCost, unit_cost=unit_cost, is_fraction=frac)
            self.__cov_func = partial(linCostCov, unit_cost=unit_cost, is_fraction=frac)

    def getDefaultBudget(self, year=None):
        '''
        Returns program cost, interpolated either for a given year or the last year that data exists for, as a budget.
        Note that users may enter a unit cost into assumptions, which will provide a flawed budget value in the absence
        of other data.
        '''

        if year is None: year = max(self.t)

        # if all values are NaN, the no costs have been provided and must be computed from the coverate
        if all(np.isnan(v) for v in self.cost):
            output = self.interpolate(tvec=[year], attributes=['cov'])
            budget = self.getBudget(output['cov'][-1])
        else:
            output = self.interpolate(tvec=[year], attributes=['cost'])
            budget = output['cost'][-1]

        return budget


    # Note that inverse functions may be more challenging than the functions employed in Program.getCoverage().
    def getBudget(self, coverage, pop_size=None):
        '''
        Returns a budget that corresponds to a program coverage. In simplest form, this is coverage times unit cost.
        Note that this method is not typically used during model processing.
        '''

        # if pop_size is provided, update current pop_size
        if pop_size is not None:
            self.pop_size = pop_size

        # check if budget function makes use of this parameter, if it does, update the value
        if 'pop_size' in self.__bud_func.keywords:
            self.__bud_func.keywords['pop_size'] = self.pop_size

        return self.__bud_func(coverage)


    def getCoverage(self, budget, pop_size=None):
        '''
        Returns prospective coverage for a program. In simplest form, this is budget divided by unit cost.
        Note that this coverage will be dynamically divided amongst populations/compartments in the model prior to
        scaling into an impact. Excess coverage will also be ignored in the model.
        '''

        # if pop_size is provided, update current pop_size
        if pop_size is not None:
            self.pop_size = pop_size

        # check if coverage function makes use of this parameter, if it does, update the value
        if 'pop_size' in self.__bud_func.keywords:
            self.__bud_func.keywords['pop_size'] = self.pop_size

        return self.__cov_func(budget)


    def getImpact(self, budget, pop_size=None, impact_label=None, parser=None, years=None, budget_is_coverage=False):
        # self.func_specs['type'] == 'cost_only':
        # if it is a cost only function, nan will be returned independent of the input values
        if np.isnan(np.sum(self.getCoverage(budget, pop_size))):
            return 0.0

        # if pop_size is provided, update current pop_size
        if pop_size is not None:
            self.pop_size = pop_size

        budget = dcp(budget)  # Just in case.

        # Baseline impact is just coverage.
        if budget_is_coverage:
            imp = budget
        else:
            imp = self.getCoverage(budget, self.pop_size)

        # If impact parameter has an impact function, this is what coverage is scaled by.
        if not impact_label is None:
            if 'f_stack' in self.target_pars[impact_label].keys():
                if parser is None:
                    raise OptimaException('ERROR: Cannot calculate "%s" impact for "%s" without a parser, due to the '
                                          'existence of an impact function.' % (self.label, impact_label))
                f_stack = dcp(self.target_pars[impact_label]['f_stack'])
                attribs = dcp(self.target_pars[impact_label]['attribs'])
                for attrib_label in attribs.keys():
                    if attrib_label in self.attributes.keys():
                        if years is None: years = [max(self.t)]
                        output = self.interpolate(tvec=years, attributes=[attrib_label])
                    else:
                        raise OptimaException(
                            'ERROR: Cannot calculate "%s" impact for "%s" due to a missing reference "%s".'
                            % (self.label, impact_label, attrib_label))
                    attribs[attrib_label] = output[attrib_label]  # [-1]
                #                    if len(output[attrib_label]) > 1:
                #                        print attribs
                new_val = parser.evaluateStack(stack=f_stack, deps=attribs)
                imp *= new_val  # Scale coverage.

        return imp


def constCostCov(*args, **kwargs):
    """Cost only function. Simply returns NaN."""
    return np.nan


def logCostCov(budget, unit_cost, saturation, pop_size, is_fraction, dt_year):
    """Logistic cost coverage function. Returns coverage, expects a budget"""
    val = (2. * saturation / (1. + np.exp(-2. * budget / (pop_size * saturation * unit_cost))) - saturation)
    if is_fraction:
        return val
    else:
        return val * pop_size / (365. * dt_year)


def logCovCost(coverage, unit_cost, saturation, pop_size, is_fraction, dt_year):
    """Logistic coverage cost function. Returns budget, expects a coverage"""
    if is_fraction:
        return - pop_size * saturation * unit_cost / 2. * np.log((saturation - coverage) / (coverage + saturation))
    else:
        return - pop_size * saturation * unit_cost / \
               2. * np.log((saturation * pop_size - coverage) / (coverage + saturation * pop_size)) * 365. * dt_year


def linCostCov(budget, unit_cost, is_fraction):
    """Linear cost coverage function. Returns coverage, expects a budget"""
    if is_fraction:
        return budget * 0.01 / unit_cost
    else:
        return budget / unit_cost


def linCovCost(coverage, unit_cost, is_fraction):
    """Linear coverage cost function. Returns budget, expects a coverage"""
    if is_fraction:
        return coverage * unit_cost / 0.01
    else:
        return coverage * unit_cost
