from optima_tb.defaults import defaultOptimOptions
from optima_tb.utils import *
from scipy.interpolate import PchipInterpolator
import pickle
from copy import deepcopy as dcp


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

    def _calculateBudgetOutcomes(self, budget_scaling, filename=None, num_iter=5):
        """
        For the given budget scaling factors, determine the spending for the optimal outcome. If filename is not None,
        the results are pickled in to a file for later reference.

        :param budget_scaling: List of floats. Each number determines how the budget is scaled before the optimal
        resource allocation is computed
        :param filename: Name of the pickle-file in which the results are stored
        :param num_iter: int which specifies how many runs with asd per region are to be computed
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

        # compute budget outcome for each region for each scaled budget
        self._BO = {}
        for region in self._opt:
            self._BO[region] = []
            for scaling in budget_scaling:
                logger.info('Optimising region \'%s\' with budget scaled by %f' % (region, scaling))
                opt = self._scaleBudget(scaling, self._opt[region])
                best = np.inf
                for it in range(num_iter):
                    _, obj_vals, _ = self._proj.optimize(parset_name=region, progset_name=region, options=opt)
                    if best > obj_vals[-1]:
                        best = obj_vals[-1]

                self._BO[region].append(best)

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

    def _optimize(self, total_national_budget, num_iter):
        """
        Determine the optimal spending per region and then optimise the spendings per region

        :param total_national_budget: total budget which is to be distributed among the regions
        :param num_iter: int which specifies how many runs with asd per region are to be computed

        """
        logger.info('Computing budget-outcome-curves..')
        self._calculateBudgetOutcomeCurves()

        logger.info('Finding the optimal budget for each region..')
        spendings = self._getSpendingsPerRegion(total_national_budget)
        opt_regional_budgets = self._runGridSearch(spendings)

        logger.info('Computing optimal outcomes based on the optimal budget')
        self._optimizeInterventions(opt_regional_budgets, num_iter)

    def _optimizeInterventions(self, opt_budgets, num_iter):
        """
        Optimise the spendings on interventions for a given budget.

        :param opt_budgets: Budget for which the spendings are optimised
        :param num_iter: int which specifies how many runs with asd per region are to be computed
        """
        # perform optimisation for each region
        for i, region in enumerate(sorted(self._opt.keys())):
            # scale original budget according to optimal budget;
            # creating a new entry in _opt would also have been a choice
            org_budget = self._opt[region]['constraints']['total']
            new_budget = opt_budgets[i]
            scaling = new_budget / org_budget
            opt = self._scaleBudget(scaling, self._opt[region])
            # optimize interventions of current region
            best_val = np.inf
            best_alloc = {}
            for it in range(num_iter):
                alloc, obj_vals, _ = \
                    self._proj.optimize(parset_name=region, progset_name=region, options=opt)
                if best_val > obj_vals[-1]:
                    best_val = obj_vals[-1]
                    best_alloc = alloc

            self._opt[region]['opt_alloc'] = best_alloc
            self._outcome[region] = best_val

    def runFromScratch(self, total_national_budget, budget_scalings, BO_file=None, num_iter=5):
        """
        Run the geospatial analyses without precomputed budget-outcome samples.

        :param total_national_budget: total budget which is to be distributed among the regions
        :param budget_scalings: list of budget scaling factors which is used to computed budget-outcome samples
        :param BO_file: If not None, the budget-outcome samples along with the budget_scalings are written to that file
        :param num_iter: int which specifies how many runs with asd per region are to be computed
        """
        logger.info('Computing budget outcomes:')
        self._calculateBudgetOutcomes(budget_scalings, BO_file, num_iter)
        self._optimize(total_national_budget, num_iter)

    def runWithBOC(self, total_national_budget, BO_file, num_iter=5):
        """
        Run the geospatial analyses using previously computed budget-outcome samples.

        :param total_national_budget: total budget which is to be distributed among the regions
        :param BO_file: Filename from which the budget-outcome samples are loaded
        :param num_iter: int which specifies how many runs with asd per region are to be computed
        """
        logger.info('Loading budget outcomes..')
        self._loadBudgetOutcomes(BO_file)
        self._optimize(total_national_budget, num_iter)

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
        with open(filename, 'wb') as f:
            pickle.dump(pickle.dumps(self._BO), f)
            pickle.dump(pickle.dumps(self._budgetScalings), f)
