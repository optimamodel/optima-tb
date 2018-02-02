from optima_tb.utils import OptimaException
from optima_tb.interpolation import interpolateFunc
from optima_tb.costcovfunction import LogisticCCF, LinearCCF, ConstCCF, CCFAttr
from copy import deepcopy as dcp
import numpy as np
from uuid import uuid4 as uuid
import logging


logger = logging.getLogger(__name__)


class ProgramSet:

    def __init__(self, name='default'):
        self.name = name
        self.uid = uuid()

        self.progs = list()
        self.prog_ids = dict()

        self.impacts = dict()

        # stores program dependencies as sorted list. sorted in this context means from independent to dependent
        self.deps = []

        logging.info("Created ProgramSet: %s" % self.name)

    def makeProgs(self, data, settings):
        for l, prog_label in enumerate(data['progs']):
            special = None
            prog_type = data['progs'][prog_label]['prog_type']
            if 'special' in settings.progtype_specs[prog_type]:
                special = settings.progtype_specs[prog_type]['special']
            data['progs'][prog_label]['special'] = special

            # Assume all program coverage/impact that targets transition parameters is 'effective', not 'probabilistic'.
            # For number format coverages/impacts, this will be uniformly distributed across timesteps.
            # For fraction format coverages/impacts, a dt-ly rate is the same as the intended annual rate.
            # This means the two formats result in differing t-dependent distributions, but the total should be the
            # same, assuming no coverage capping.
            if data['progs'][prog_label]['cov_format'] == 'number':
                data['progs'][prog_label]['cov'] *= settings.tvec_dt

            target_pars = settings.progtype_specs[prog_type]['impact_pars']
            for target_par in target_pars.keys():
                if target_par not in self.impacts: self.impacts[target_par] = []
                self.impacts[target_par].append(prog_label)
            data['progs'][prog_label]['target_pars'] = target_pars

            alive_tag = settings.charac_pop_count
            target_pop_size = np.sum([data['characs'][alive_tag][pop]['y'][0]
                                      for pop in data['progs'][prog_label]['target_pops']])
            data['progs'][prog_label]['target_pop_size'] = target_pop_size

            data['progs'][prog_label]['dt'] = settings.tvec_dt

            self.progs.append(Program(prog_label, data['progs'][prog_label]))
            self.prog_ids[prog_label] = l

        # figure out in which order the programs can be processed without violating dependencies:
        # extract all independent programs -> prog.deps is None
        dep_set = set([prog.label for prog in self.progs if prog.deps is None])
        if not dep_set:
            # if there are no dependencies, there is no use in trying to resolve dependencies
            return

        # complementary set which include all dependent programs
        rem_set = set([prog.label for prog in self.progs if prog.label not in dep_set])
        prev_len = -1

        # until all dependencies are resolved..
        while len(dep_set) != len(self.progs):
            # .. check if the number of resolved dependencies has increased:
            if prev_len == len(dep_set):
                # if not, there is a cyclic dependency and further iteration will not help
                raise OptimaException('Cyclic dependency can not be resolved between the two sets'
                                      + '"{}" and "{}"'.format(list(dep_set), list(rem_set)))
            else:
                # if so, update the number of resolved dependencies
                prev_len = len(dep_set)

            # resolve dependencies of all the programs whose dependencies have not be resolved yet
            for prog in filter(lambda x: x.label in rem_set, self.progs):
                # resolve dependency: a dependency is considered resolved if all programs on which it is dependent have
                # their dependencies resolved
                if prog.deps.issubset(dep_set):
                    # add program if dependency was resolved
                    dep_set.update(prog_label)

            # update remaining set by removing all newly resolved dependencies
            rem_set.difference_update(dep_set)

        self.deps = list(dep_set)

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

    def __init__(self, label, ppd):
        """

        :param label: str which denotes the label of the program
        :param ppd: dict which is generate by parsing the databook limited to the section of one program. This
            abbreviation is for 'program property dictionary'.
        """
        self.name = dcp(ppd['name'])
        self.label = dcp(label)
        self.prog_type = dcp(ppd['prog_type'])
        self.uid = uuid()

        self.t = dcp(ppd['t'])  # Time data.

        self.cost = dcp(ppd['cost'])  # Spending data.
        self.cost_format = dcp(ppd['cost_format'])

        self.cov = dcp(ppd['cov'])  # Coverage data.
        self.cov_format = dcp(ppd['cov_format'])

        self.attributes = dcp(ppd['attributes'])

        self.target_pops = dcp(ppd['target_pops'])
        self.target_pars = dcp(ppd['target_pars'])

        self.deps = dcp(ppd['deps']) if 'deps' in ppd else None
        self.dep_types = dcp(ppd['dep_type']) if 'dep_type' in ppd else None

        sat = ppd['sat'] if 'sat' in ppd else None
        unit_cost = ppd['unit_cost'] if 'unit_cost' in ppd else None

        if ppd['special'] == 'cost_only':
            self.ccf = ConstCCF()
        elif any([sat is None, ppd['target_pop_size'] is None]):
            # backward compatibility: linear cost coverage curves
            self.ccf = LinearCCF(unit_cost=unit_cost)
        else:
            cif = True if self.cov_format == 'fraction' else False
            self.ccf = LogisticCCF(unit_cost, sat, ppd['target_pop_size'], ppd['dt'], cif=cif)
        
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
#                logging.info('Inserted/replaced value "%f" at year "%f" for attribute "%s" of program "%s"' % (new_val, t_val, attribute, self.label))
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
        except: raise OptimaException('ERROR: Unable to insert value "%f" at year "%f" for attribute "%s" of program "%s"' % (y, t, attribute, self.label))

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

    def getDefaultBudget(self, year=None):
        '''
        Returns program cost, interpolated either for a given year or the last year that data exists for, as a budget.
        Note that users may enter a unit cost into assumptions, which will provide a flawed budget value in the absence of other data.
        '''

        if year is None: year = max(self.t)
        output = self.interpolate(tvec=[year], attributes=['cost'])
        budget = output['cost'][-1]
        return budget

    # Note that inverse functions may be more challenging than the functions employed in Program.getCoverage().
    def getBudget(self, coverage, pop_size=None):
        '''
        Returns a budget that corresponds to a program coverage. In simplest form, this is coverage times unit cost.
        Note that this method is not typically used during model processing.
        '''

        # If coverage is a scalar, make it a float. If it is a list or array, make it an array of floats.
        try: coverage = float(coverage)
        except: coverage = dcp(np.array(coverage, 'float'))

        # if pop_size is provided, update current pop_size
        if pop_size is not None:
            params = {CCFAttr.pop: pop_size}
        else:
            params = {}

        return self.ccf.costs(coverage, params)

    def getCoverage(self, budget, pop_size=None):
        '''
        Returns prospective coverage for a program. In simplest form, this is budget divided by unit cost.
        Note that this coverage will be dynamically divided amongst populations/compartments in the model prior to scaling into an impact.
        Excess coverage will also be ignored in the model.
        '''

        # If budget is a scalar, make it a float. If it is a list or array, make it an array of floats.
        try: budget = float(budget)
        except: budget = dcp(np.array(budget, 'float'))

        # if pop_size is provided, update current pop_size
        if pop_size is not None:
            params = {CCFAttr.pop: pop_size}
        else:
            params = {}

        return self.ccf.coverage(budget, params)


    def getImpact(self, budget, pop_size=None, impact_label=None, parser=None, years=None, budget_is_coverage=False):
        # if it is a cost only function, nan will be returned independent of the input values
        if np.isnan(np.sum(self.getCoverage(budget, pop_size))):
            return 0.0

        # if pop_size is provided, update current pop_size
        if pop_size is not None:
            self.pop_size = pop_size

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
                f_stack = self.target_pars[impact_label]['f_stack']
                attribs = self.target_pars[impact_label]['attribs']
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
