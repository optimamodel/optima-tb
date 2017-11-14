from optima_tb.defaults import defaultOptimOptions
from optima_tb.project import Project
from optima_tb.utils import *
from scipy.interpolate import PchipInterpolator
import os, pickle
from copy import deepcopy as dcp


class GeospatialOptimization:
    """
    GeospatialOptimization initialises one project with one cascade and many parameter and program sets, which represent
    the geospatially different regions. Everything is set up automatically so that only one of the run*() functions
    needs to be called for the actual optimisation.
    """

    def __init__(self, cascade_path, region_spreadsheets, region_names, obj_labels, I, pars_interp_num=None):
        """

        :param cascade_path: filename of the cascade spreadsheet. All regions use the same cascade structure
        :param region_spreadsheets: list of filenames which contain the data about each region
        :param region_names: list of names/labels for the regions.
        :param obj_labels: list of parameter labels which is used to measure the outcome
        :param I: Simulation interval (SimInt object) of the geospatial optimisation
        :param pars_interp_num: List of labels (or None) of parameters which have to be reinterpreted as numbers
        """
        # sanity checks:
        # region_spreadsheets, region_names and obj_labels must be lists
        if not isinstance(region_spreadsheets, list):
            region_spreadsheets = [region_spreadsheets]
        if not isinstance(region_names, list):
            region_names = [region_names]
        if not isinstance(obj_labels, list):
            obj_labels = [obj_labels]

        # number of spreadsheets must coincide with the number of names provided
        assert(len(region_spreadsheets) == len(region_names))

        # setup project and corresponding programs
        self._proj = Project('geo_analysis', cascade_path, I, pars_interp_num)

        # initialise options for each region and set objective values
        self._opt = {}
        for rname, rpath in zip(region_names, region_spreadsheets):
            self._proj.addRegion(rpath, rname) # create parameters for the model and setup programs
            self._opt[rname] = defaultOptimOptions(self._proj.settings, self._proj.progsets[rname], I.start)
            for obj in obj_labels:
                self._opt[rname]['objectives'].clear()
                # currently only minimisation implemented; for maximisation, a way to pass a weight for the dict below
                # has to be implemented
                self._opt[rname]['objectives'][obj] = {'weight': 1., 'year': I.stop}

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

    def _calculateBudgetOutcomes(self, budget_scaling, filename=None):
        """
        For the given budget scaling factors, determine the spending for the optimal outcome. If filename is not None,
        the results are pickled in to a file for later reference.

        :param budget_scaling: List of floats. Each number determines how the budget is scaled before the optimal
        resource allocation is computed
        :param filename: Name of the pickle-file in which the results are stored
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
                _, obj_vals, _ = self._proj.optimize(parset_name=region, progset_name=region, options=opt)
                self._BO[region].append(obj_vals[-1])

        if filename is not None:
            logger.info('Saving computed budget outcomes in %s' % (filename))
            self._writeBudgetOutcomesToFile(filename)

    def _calculateBudgetOutcomeCurves(self):
        """
        Calculate the budget-outcome-curves from the budget-outcome samples
        """
        for region in sorted(self._BO.keys()):
            self.BOC[region] = PchipInterpolator(self._budgetScalings, self._BO[region])

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

    def getGAResults(self):
        """
        Retrieve the geospatial optimisation result.
        :return: dict with the optimal spendings per region and the corresponding outcome per region
        """
        if 'opt_alloc' not in self._opt:
            raise OptimaException(('No results computed. Call runFromScratch() or runWithBO()'
                                   'before calling getGAResults!'))

        results = {}
        for region in self._opt:
            results[region] = {'alloc': self._opt[region]['opt_alloc'], 'outcome': self._outcome[region]}

        return results

    def _loadBudgetOutcomes(self, filename):
        with open(filename, 'rb') as f:
            self._BO = pickle.loads(pickle.load(f))
            self._budgetScalings = pickle.loads(pickle.load(f))

    def _optimizeInterventions(self, opt_budgets):
        """
        Optimise the spendings on interventions for a given budget.
        :param opt_budgets: Budget for which the spendings are optimised
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
            self._opt[region]['opt_alloc'], obj_vals, _ = \
                self._proj.optimize(parset_name=region, progset_name=region, options=opt)
            self._outcome[region] = obj_vals[-1]

    def runFromScratch(self, total_national_budget, budget_scalings, BO_file=None):
        """
        Run the geospatial analyses without precomputed budget-outcome samples.
        :param total_national_budget: total budget which is to be distributed among the regions
        :param budget_scalings: list of budget scaling factors which is used to computed budget-outcome samples
        :param BO_file: If not None, the budget-outcome samples along with the budget_scalings are written to that file
        """
        logger.info('Computing budget outcomes:')
        self._calculateBudgetOutcomes(budget_scalings, BO_file)
        self._optimize(total_national_budget)

    def runWithBOC(self, total_national_budget, BO_file):
        """
        Run the geospatial analyses using previously computed budget-outcome samples.
        :param total_national_budget: total budget which is to be distributed among the regions
        :param BO_file: Filename from which the budget-outcome samples are loaded
        """
        logger.info('Loading budget outcomes..')
        self._loadBudgetOutcomes(BO_file)
        self._optimize(total_national_budget)

    def _optimize(self, total_national_budget):
        """
        Determine the optimal spending per region and then optimise the spendings per region
        :param total_national_budget: total budget which is to be distributed among the regions
        """
        logger.info('Computing budget-outcome-curves..')
        self._calculateBudgetOutcomeCurves()

        logger.info('Finding the optimal budget for each region..')
        spendings = self._getSpendingsPerRegion(total_national_budget)
        opt_regional_budgets = self._runGridSearch(spendings)

        logger.info('Computing optimal outcomes based on the optimal budget')
        self._optimizeInterventions(opt_regional_budgets)

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
