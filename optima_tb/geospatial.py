from optima_tb.defaults import defaultOptimOptions
from optima_tb.utils import *
from scipy.interpolate import PchipInterpolator
import pickle
from copy import deepcopy as dcp
import copy_reg

class BudgetOutcomeCurve:
    """
    Represents a piecewise cubic hermite interpolating polynomial (pchip) which is used in the grid search algorithm. It
    contains the function itself but also its first derivative and functionality to shift the function
    up/down/left/right.
    """

    def __init__(self, pchip):
        """
        :param: PchipInterpolator(scipy.interpolate.pchip)
        """
        # store the interpolated polynomial
        self._pfun = dcp(pchip)
        # store the first derivative of the interpolated polynomial
        self._pder = self._pfun.derivative(1)
        # initial shifts: none
        self._xshift = 0.
        self._yshift = 0.

    def __call__(self, x):
        """
        Shorthand for function evaluation at point x including previous shifts.
        :param x: float
        :return: Function value of the interpolated polynomial
        """
        return self.eval(x)

    def eval(self, x):
        """
        Evaluates the interpolated polynomial at point x including previous shifts.

        :param x: float
        :return: Function value of the interpolated polynomial
        """
        return self._pfun(x - self._xshift) + self._yshift

    def deriv(self, x):
        """
        Evaluates the first derivative of the interpolated polynomial at point x including previous shifts.

        :param x: float
        :return: Function value of the first derivative of the interpolated polynomial
        """
        return self._pder(x - self._xshift) + self._yshift

    def shift(self, left, down):
        """
        Shifts the interpolated polynomial and its derivative left and down. All shifts are cumulative!

        :param left: float which shifts the curve to the left (i.e. negative value shifts it to the right)
        :param down: float which shifts the curve down (i.e. negative value shifts it up)
        """
        self._xshift += left
        self._yshift -= down


@trace_exception
def _call_project_optimize(task):
    """
    Wrapper around the budget optimisation algorithm Project.optimize(). This wrapper is used to provide an
    interface to parallel computation.


    :param task: tuple of (Project, kwargs), where Project is the object on which the optimisation is performed and
        kwargs are arguments required by the Project.optimize() function
    :return: tuple of (allocation, outcome), where allocation is the budget allocation across interventions and
        outcome is the optimal outcome
    """
    proj = task[0]
    kwargs = task[1]
    params, obj_vals, _ = proj.optimize(**kwargs)
    return params, obj_vals[-1]


class GeospatialOptimization:
    """
    GeospatialOptimization initialises one project with one cascade and many parameter and program sets, which represent
    the geospatially different regions. Everything is set up automatically so that only one of the run*() functions
    needs to be called for the actual optimisation.
    """

    def __init__(self, proj, obj_labels):
        """
        :param proj: Project for the geospatial optimisation
        :param obj_labels: list of parameter labels which is used to measure the outcome
        """
        # sanity checks: obj_labels must be lists
        if not isinstance(obj_labels, list):
            obj_labels = [obj_labels]

        # setup project and corresponding programs
        self._proj = dcp(proj)

        # initialise options for each region and set objective values
        self._opt = {}
        for rname in self._proj.progsets:
            self._opt[rname] = \
                defaultOptimOptions(self._proj.settings, self._proj.progsets[rname], self._proj.settings.tvec_start)
            for obj in obj_labels:
                self._opt[rname]['objectives'].clear()
                # currently only minimisation implemented; for maximisation, a way to pass a weight for the dict below
                # has to be implemented
                self._opt[rname]['objectives'][obj] = {'weight': 1., 'year': self._proj.settings.tvec_end}

        # define class members:
        # - defines how the default budget is scaled for the budget outcome samples.
        # Please NOTE this list is SORTED
        self._budgetScalings = []
        # - budget-outcome-samples
        self._BO = {}
        # - interpolated budget-outcome-curves which are computed from the budget-outcome-samples
        self._BOC = {}
        # - optimal outcomes per region
        self._outcome = {}

    def _calculateBudgetOutcomes(self, budget_scaling, num_threads, filename=None, num_iter=5):
        """
        For the given budget scaling factors, determine the spending for the optimal outcome. If filename is not None,
        the results are pickled in to a file for later reference.

        :param budget_scaling: List of floats. Each number determines how the budget is scaled before the optimal
        resource allocation is computed
        :param filename: Name of the pickle-file in which the results are stored
        :param num_iter: int which specifies how many runs with asd per region are to be computed
        :param num_threads: int which specify how many threads should be used for computation
        """
        # sanity checks:
        # - make sure budget_scaling is a list
        if not isinstance(budget_scaling, list):
            budget_scaling = [budget_scaling]
        # - make sure budget_scaling has no double entries
        if len(set(budget_scaling)) != len(budget_scaling):
            raise OptimaException('There are duplicate entries in the list to scale the budgets. Please remove them!')

        # sort entries: it is required by the interpolation. So do it now to prevent ordering issues after the
        #  interpolation
        self._budgetScalings = sorted(budget_scaling)

        # create list which contains the tasks to process
        tasks = self._createBOCalculationTasks(num_iter)

        best_results = self._processTasks(tasks, num_iter, num_threads)

        # create a dict which contains each region and an entry for each scaled budget; must be in the same order as in
        # _createBOCalculationTasks()
        for region in self._opt:
            self._BO[region] = []
            for scaling in self._budgetScalings:
                self._BO[region].append(best_results.pop(0)[1])

        if filename is not None:
            logger.info('Saving computed budget outcomes in %s' % filename)
            self._writeBudgetOutcomesToFile(filename)

    def _calculateBudgetOutcomeCurves(self):
        """
        Calculate the budget-outcome-curves from the budget-outcome samples
        """
        for region in sorted(self._BO.keys()):
            self._BOC[region] = BudgetOutcomeCurve(PchipInterpolator(
                [x * self._opt[region]['constraints']['total'] for x in self._budgetScalings], self._BO[region]))

    def _createBOCalculationTasks(self, num_iter, budget_scaling=None):
        """
        Create a list of dicts which can be passed to the optimisation function for parallel(or serial)
        processing. In this case it is the Project.optimize-function.

        CAUTION: This function creates a deep copy of the project for each task. This is to ensure no inconsistencies
        occur during parallel processing by modifying the project. If it can be guaranteed that the Project object is
        not modified during Project.optimize, no deep copies are necessary.

        :param num_iter: int which specifies how often the budget outcome should be computed
        :param budget_scaling: None or list of floats which specifies if a specific budget scaling is used (list) or all
            budgetScalings should be used (None). First case is used when an optimal spending already exists, in fact
            the optimal spending scaling is the parameter 'scaling'; the second case is used to determine budget outcome
            samples from which BOCs are computed
        """
        params = []
        if budget_scaling is None:
            # in this case apply each scaling to each region num_iter times and optimise
            for region in self._opt:
                for scaling in self._budgetScalings:
                    opt = self._scaleBudget(scaling, self._opt[region])
                    for it in range(num_iter):
                        params.append((self._proj, {'parset_name': region, 'progset_name': region, 'options': opt}))
        else:
            assert(len(budget_scaling) == len(self._opt))
            for i, region in enumerate(sorted(self._opt.keys())):
                # scale original budget according to optimal budget;
                # creating a new entry in _opt would also have been a choice
                org_budget = self._opt[region]['constraints']['total']
                new_budget = budget_scaling[i]
                scaling = new_budget / org_budget
                opt = self._scaleBudget(scaling, self._opt[region])
                for it in range(num_iter):
                    params.append((self._proj, {'parset_name': region, 'progset_name': region, 'options': opt}))
        return params


    def _getBOCderivatives(self, currentRegionalSpending, extraFunds):
        from numpy import linspace
        numPoints = 2000
        spendingVec = []
        outcomeVec = []
        costEffVecs = []
        numRegions = len(self._BO.keys())
        # regions are stored w.r.t. the region label, not the index. create a mapping which associates an index to the
        # corresponding label
        idx2region = dict(zip(range(len(self._BO.keys())), sorted(self._BO.keys())))

        # calculate cost effectiveness vectors for each region
        for regionIdx in range(numRegions):
            if extraFunds:  # only consider additional funds
                regionalMaxBudget = currentRegionalSpending[regionIdx] + extraFunds
                regionalMinBudget = currentRegionalSpending[regionIdx]
            else:
                regionalMaxBudget = sum(currentRegionalSpending)
                regionalMinBudget = 0.
            spending = linspace(regionalMinBudget, regionalMaxBudget, numPoints)
            spendingThisRegion = spending[1:]  # exclude 0 to avoid division error
            adjustedRegionalSpending = spendingThisRegion - regionalMinBudget  # centers spending if extra funds
            spendingVec.append(adjustedRegionalSpending)
            regionalBOC = self._BOC[idx2region[regionIdx]]
            outcomeThisRegion = regionalBOC(spendingThisRegion)
            outcomeVec.append(outcomeThisRegion)
            costEffThisRegion = outcomeThisRegion / spendingThisRegion
            costEffVecs.append(costEffThisRegion)
        return costEffVecs, spendingVec, outcomeVec

    def getGAResults(self):
        """
        Retrieve the geospatial optimisation result.

        :return: dict with the optimal spendings per region and the corresponding outcome per region
        """
        if not self._outcome:
            raise OptimaException(('No results computed. Call runFromScratch() or runWithBO()'
                                   'before calling getGAResults!'))

        alloc = {}
        outcome = {}
        for region in self._opt:
            alloc[region] = self._opt[region]['opt_alloc']
            outcome[region] = self._outcome[region]

        return alloc, outcome

    def _getSpendingsPerRegion(self, budget):
        """
        :param budget: Total budget available across all regions
        :return: If there are N regions, then create a list of N entries with each entry being budget/N
        """
        # uniformally distribute the budget on the regions
        num_regions = len(self._BO)
        budget /= float(num_regions)
        spendings = [budget] * num_regions
        return spendings

    def _loadBudgetOutcomes(self, filename):
        with open(filename, 'rb') as f:
            self._BO = pickle.loads(pickle.load(f))
            self._budgetScalings = pickle.loads(pickle.load(f))

    def _optimize(self, total_national_budget, num_iter, num_threads):
        """
        Determine the optimal spending per region and then optimise the spendings per region

        :param total_national_budget: total budget which is to be distributed among the regions
        :param num_iter: int which specifies how many runs with asd per region are to be computed
        :param num_threads: int which specify how many threads should be used for computation
        """
        logger.info('Computing budget-outcome-curves..')
        self._calculateBudgetOutcomeCurves()

        logger.info('Finding the optimal budget for each region..')
        spendings = self._getSpendingsPerRegion(total_national_budget)
        opt_regional_budgets = self._runGridSearch(spendings)

        logger.info('Computing optimal outcomes based on the optimal budget')
        self._optimizeInterventions(opt_regional_budgets, num_iter, num_threads)

    def _optimizeInterventions(self, opt_budgets, num_iter, num_threads):
        """
        Optimise the spendings on interventions for a given budget.

        :param opt_budgets: Budget for which the spendings are optimised
        :param num_iter: int which specifies how many runs with asd per region are to be computed
        :param num_threads: int which specify how many threads should be used for computation
        """
        tasks = self._createBOCalculationTasks(num_iter, opt_budgets)

        best_results = self._processTasks(tasks, num_iter, num_threads)

        # create a dict which contains each region and an entry for each scaled budget; must be in the same order as in
        # _createBOCalculationTasks()
        for region in self._opt:
            obj_param, obj_val = best_results.pop(0)
            self._opt[region]['opt_alloc'] = obj_param
            self._outcome[region] = obj_val

    def _processTasks(self, tasks, num_iter, num_threads):
        """
        Enqueues the passed tasks and processes them sequentially, either in parallel or serially. After the tasks
        are done, the best task, i.e. the one with the smallest outcome, per region and/or budget is chosen.

        :param tasks: list of tuple generated by _createBOCalculationTasks()
        :param num_iter: int which specifies how many runs with asd per region are to be computed
        :param num_threads: int which specify how many threads should be used for computation
        :return: list of best outcomes per region and/or budget
        """
        # process tasks serially or in parallel
        results = pmap(_call_project_optimize, tasks, num_threads)

        # extract the best result, i.e. the maximum outcome, from all computed scenarios:
        # due to the process which is used to set up the task, it is possible to split the array of results into arrays
        # which contain the results of only a region with a given budget (that array has a length of num_iter). On each
        # of these arrays the min()-operation is performed to obtain the best outcome.
        return map(lambda x: min(x, key=lambda y: y[1]), np.array_split(results, len(tasks) / num_iter))

    def runFromScratch(self, total_national_budget, budget_scalings, BO_file=None, num_iter=5, num_threads=1):
        """
        Run the geospatial analyses without precomputed budget-outcome samples.

        :param total_national_budget: total budget which is to be distributed among the regions
        :param budget_scalings: list of budget scaling factors which is used to computed budget-outcome samples
        :param BO_file: If not None, the budget-outcome samples along with the budget_scalings are written to that file
        :param num_iter: int which specifies how many runs with asd per region are to be computed
        :param num_threads: int which specify how many threads should be used for computation
        """
        logger.info('Computing budget outcomes:')
        self._calculateBudgetOutcomes(budget_scalings, num_threads, BO_file, num_iter)
        self._optimize(total_national_budget, num_iter, num_threads)

    def runWithBOC(self, total_national_budget, BO_file, num_iter=5, num_threads=1):
        """
        Run the geospatial analyses using previously computed budget-outcome samples.

        :param total_national_budget: total budget which is to be distributed among the regions
        :param BO_file: Filename from which the budget-outcome samples are loaded
        :param num_iter: int which specifies how many runs with asd per region are to be computed
        :param num_threads: int or None which specify how many threads should be used for computation
        """
        logger.info('Loading budget outcomes..')
        self._loadBudgetOutcomes(BO_file)
        self._optimize(total_national_budget, num_iter, num_threads)

    def _runGridSearch(self, currentRegionalSpending, extraFunds=None):
        """ If specified, extraFunds is a scalar value which represents funds to be distributed on top of fixed current
        regional spending """
        from numpy import zeros, inf, nonzero, argmax
        from copy import deepcopy as dcp
        costEffVecs, spendingVec, outcomeVec = self._getBOCderivatives(currentRegionalSpending, extraFunds=extraFunds)
        if extraFunds:
            totalBudget = extraFunds
        else:
            totalBudget = sum(currentRegionalSpending)
        numRegions = len(self._BO.keys())
        maxIters = int(1e6)  # loop should break before this is reached
        remainingFunds = dcp(totalBudget)
        regionalAllocations = zeros(numRegions)
        percentBudgetSpent = 0.

        for i in range(maxIters):
            # find the most cost-effective region to fund
            bestEff = -inf
            bestRegion = None
            for regionIdx in range(numRegions):
                # find most effective spending in each region
                costEffThisRegion = costEffVecs[regionIdx]
                if len(costEffThisRegion):
                    maxIdx = argmax(costEffThisRegion)
                    maxEff = costEffThisRegion[maxIdx]
                    if maxEff > bestEff:
                        bestEff = maxEff
                        bestEffIdx = maxIdx
                        bestRegion = regionIdx

            # once the most cost-effective spending is found, adjust all spending and outcome vectors, update available
            #  funds and regional allocation
            if bestRegion is not None:
                fundsSpent = spendingVec[bestRegion][bestEffIdx]
                remainingFunds -= fundsSpent
                spendingVec[bestRegion] -= fundsSpent
                outcomeVec[bestRegion] -= outcomeVec[bestRegion][bestEffIdx]
                regionalAllocations[bestRegion] += fundsSpent
                # remove funds and outcomes at or below zero
                spendingVec[bestRegion] = spendingVec[bestRegion][bestEffIdx + 1:]
                outcomeVec[bestRegion] = outcomeVec[bestRegion][bestEffIdx + 1:]
                # ensure regional spending doesn't exceed remaining funds
                for regionIdx in range(numRegions):
                    withinBudget = nonzero(spendingVec[regionIdx] <= remainingFunds)[0]
                    spendingVec[regionIdx] = spendingVec[regionIdx][withinBudget]
                    outcomeVec[regionIdx] = outcomeVec[regionIdx][withinBudget]
                    costEffVecs[regionIdx] = outcomeVec[regionIdx] / spendingVec[regionIdx]
                newPercentBudgetSpent = (totalBudget - remainingFunds) / totalBudget * 100.
                if not (i % 100) or (newPercentBudgetSpent - percentBudgetSpent) > 1.:
                    percentBudgetSpent = newPercentBudgetSpent
            else:
                break  # nothing more to allocate

        # scale to ensure correct budget
        scaleRatio = totalBudget / sum(regionalAllocations)
        rescaledAllocation = [x * scaleRatio for x in regionalAllocations]
        return rescaledAllocation

    def _scaleBudget(self, scaling, options):
        """
        Scale all values in the options-dict by scaling which relate to the budget.

        :param scaling: float to scale the budget
        :param options: dict which contains budget allocation information
        :return: options-dict with scaled budget
        """
        opt = dcp(options)

        opt['constraints']['total'] *= scaling
        for prog in options['init_alloc']:
            options['init_alloc'][prog] *= scaling

        return opt

    def _writeBudgetOutcomesToFile(self, filename):
        """
        Write budget outcomes and the scalings in a pickle-file.
        """
        with open(filename, 'wb') as f:
            pickle.dump(pickle.dumps(self._BO), f)
            pickle.dump(pickle.dumps(self._budgetScalings), f)
