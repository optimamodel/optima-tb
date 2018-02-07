import abc
import numpy as np
from enum import Enum


class CCFAttr(Enum):
    """
    Define attributes of cost-coverage-functions.
    """
    cif = 1     # coverage is interpreted as fraction (boolean)
    ine = 2     # ignore numerical errors (boolean). States if CCFValidation.valid_number is checked or not
    uc = 3      # unit cost (non-negative)
    sat = 4     # saturation (between 0 and 1)
    pop = 5     # population size (non-negative)
    dt = 6      # time step size in years (non-negative)


class CostCovFunction(object):
    """
    This class defines the common interface of cost-coverage-functions. Standard functionality includes
     * calculating the coverage w.r.t. costs, which is done in function *coverage()*
     * calculating the costs w.r.t. coverage, which is done in function *costs()*
     * update cost-coverage-function specific attributes, which is done in function *update()*
     * a flag which determines if the coverage is interpreted as fraction or as absolute number

    To implement a concrete cost-coverage-function, derive from this class and implement
     * __costsFraction()
     * __costsNotFraction()
     * __coverageFraction()
     * __coverageNotFraction()
     * store all required attributes for the calculation. To do that use _setAttribute() and provide
        CCFAttr and the value.

    How to implement new CCFs?
     * Derive a class from CostCovFunction.
     * Implement the four abstract functions
     * Are other parameters required? All parameters are defined in CCFAttr. If another variable is required, add it
        to CCFAttr. If this happens, make sure to add an entry to _sanityChecks

    NOTES:
     * __costsFraction() and __coverageFraction(), and __costsNotFraction() and __coverageNotFraction() must be inverse
        functions, respectively
     * Functions to implement in concrete classes are 'private' and cannot be called from outside the class. The way to
        access them is via costs() and coverage() functions (which are 'public'). Sanity checks are performed in these
        public functions. There is thus no need to implement it in the derived classes.
    """
    __metaclass__ = abc.ABCMeta

    class CCFValidation(Enum):
        """
        Define how to sanity-check a particular value.
        """
        nonnegative = 1  # Value must be non-negative
        fraction = 2  # Value must be between 0 and 1
        boolean = 3  # Value must be boolean
        valid_number = 4  # Value must not be Inf or NaN

    # define attributes and valid range (as specified in the comments of CCFAttr). This is valid for all derived classes
    # so there is no need to set these mapping explicitly in each derived class
    _sanityChecks = {CCFAttr.cif:   CCFValidation.boolean,
                     CCFAttr.ine:   CCFValidation.boolean,
                     CCFAttr.uc:    CCFValidation.nonnegative,
                     CCFAttr.sat:   CCFValidation.fraction,
                     CCFAttr.pop:   CCFValidation.nonnegative,
                     CCFAttr.dt:    CCFValidation.nonnegative}

    def __init__(self, cif=True):
        """
        cif determines how the coverage is determined: as fraction (True) or as absolute value (False).
        Numerical errors are not ignored by default.
        """
        # maps attribute labels (CCFAttr) to their respective values
        self.attr = {}

        # ine should be assigned first to prevent errors happening due to other parameters relying on ine
        self._setAttribute(CCFAttr.ine, False)
        self._setAttribute(CCFAttr.cif, cif)

    @abc.abstractmethod
    def _costsFraction(self, coverage):
        """
        Function which returns the costs for a specified coverage, which is given as fraction.

        :param coverage: float representing the coverage.
        :return: float which is the cost to achieve the specified coverage.
        """
        raise NotImplementedError('__costsFraction() of abstract base class "CostCovFunction" must not be called!')

    @abc.abstractmethod
    def _costsNotFraction(self, coverage):
        """
        Function which returns the costs for a specified coverage, which is given as absolute number.

        :param coverage: float representing the coverage.
        :return: float which is the cost to achieve a specified coverage.
        """
        raise NotImplementedError('__costsNotFraction() of abstract base class "CostCovFunction" must not be called!')

    @abc.abstractmethod
    def _coverageFraction(self, costs):
        """
        Function which returns the fractional coverage for specified costs.

        :param costs: float representing the costs.
        :return: float which is the fractional coverage (value between 0 and 1) obtained by the specified costs.
        """
        raise NotImplementedError('__coverageFraction() of abstract base class "CostCovFunction" must not be called!')

    @abc.abstractmethod
    def _coverageNotFraction(self, costs):
        """
        Function which returns the absolute coverage for specified costs.

        :param costs: float representing the costs.
        :return: float which is the absolute coverage (value greater than 0) obtained by the specified costs.
        """
        raise NotImplementedError('__coverageNotFraction() of abstract base class "CostCovFunction" must not be called!')

    def _checkParameter(self, value, validation_type):
        """
        Check if the passed value satisfies the passed criterion. This function returns nothing and only raises an
        exception if anything goes wrong.

        :param value: float or boolean
        :param validation_type: CostCovFunction.CCFValidation which defined valid range and type
        """
        all_good = True
        if validation_type == CostCovFunction.CCFValidation.nonnegative and value < 0.:
            all_good = False
        elif validation_type == CostCovFunction.CCFValidation.fraction and not (0. <= value <= 1.):
            all_good = False
        elif validation_type == CostCovFunction.CCFValidation.boolean and not isinstance(value, bool):
            all_good = False

        if not self.attr[CCFAttr.ine] and (np.isnan(value) or np.isinf(value)):
            # this is done to pass the information of the error to _createErrorMessage
            validation_type = CostCovFunction.CCFValidation.valid_number
            all_good = False
        elif not all_good and self.attr[CCFAttr.ine] and (np.isnan(value) or np.isinf(value)):
            # in this case one of the previous validations fails due to NaN or Inf. Since ine is set in this case,
            # this error is ignored -> all good :)
            all_good = True

        if not all_good:
            raise ValueError(self._createErrorMessage(value, validation_type))

    def _createErrorMessage(self, value, violation_type):
        """
        Generates a string containing error information. This function is not meant to be called separately but from
        _checkParamter().

        NOTE: No error handling is performed, merely the error message is generated.

        :param value: float or boolean
        :param violation_type: CostCovFunction.CCFValidation which is violated
        :return: string with error message
        """
        error_msg = 'In {}: {} was passed, but '.format(self.__class__.__name__, value)

        if violation_type == CostCovFunction.CCFValidation.nonnegative:
            error_msg += 'a non-negative value is expected.'
        elif violation_type == CostCovFunction.CCFValidation.fraction:
            error_msg += 'a value between 0 and 1 is expected.'
        elif violation_type == CostCovFunction.CCFValidation.boolean:
            error_msg += 'a boolean is expected.'
        elif violation_type == CostCovFunction.CCFValidation.valid_number:
            error_msg += 'an actual numerical value is expected.'

        return error_msg

    def _setAttribute(self, attr_label, attr_value):
        """
        Set the variable ('label') to the passed value. Sanity check is performed.

        Technically, a dict-entry with the key *attr_label* and value *attr_value* is generated.

        :param attr_label: CCFAttr which defines the attribute
        :param attr_value: float/boolean which is assigned to the specified attribute
        """
        self.attr[attr_label] = attr_value
        self._checkParameter(attr_value, self._sanityChecks[attr_label])

    def costs(self, coverage, params={}):
        """
        Determine the costs w.r.t. the passed coverage. If params is defined, the attributes in params are updated
        *before* the coverage is determined. Any change in an attribute will be saved until it is changed again.

        :param coverage: float representing the coverage.
        :param params: dict of {attribute: value} with new attributes for the CostCovFunction.
        :return: float which represents the costs of the specified coverage.
        """
        self.update(params)

        if self.attr[CCFAttr.cif]:
            self._checkParameter(coverage, CostCovFunction.CCFValidation.fraction)
            cost = self._costsFraction(coverage)
        else:
            self._checkParameter(coverage, CostCovFunction.CCFValidation.nonnegative)
            cost = self._costsNotFraction(coverage)

        self._checkParameter(cost, CostCovFunction.CCFValidation.nonnegative)

        return cost

    def coverage(self, costs, params={}):
        """
        Determine the coverage w.r.t. the passed costs. If params is defined, the attributes in params are updated
        *before* the coverage is determined. Any change in an attribute will be saved until it is changed again.

        :param costs: float representing the costs.
        :param params: dict of {attribute: value} with new attributes for the CostCovFunction.
        :return: float which represents the coverage for the specified costs.
        """
        self.update(params)

        self._checkParameter(costs, CostCovFunction.CCFValidation.nonnegative)

        if self.attr[CCFAttr.cif]:
            coverage = self._coverageFraction(costs)
            self._checkParameter(coverage, CostCovFunction.CCFValidation.fraction)
        else:
            coverage = self._coverageNotFraction(costs)
            self._checkParameter(coverage, CostCovFunction.CCFValidation.nonnegative)

        return coverage

    def update(self, params):
        for attr in filter(lambda x: x in params, self.attr):
            self._checkParameter(params[attr], self._sanityChecks[attr])
            self.attr[attr] = params[attr]


class ConstCCF(CostCovFunction):
    """
    Implementation of a constant cost-coverage-function. It only returns NaN.
    """
    def __init__(self, cif=True):
        super(ConstCCF, self).__init__(cif)
        self._setAttribute(CCFAttr.ine, True)

    def _costsFraction(self, coverage):
        return np.nan

    def _costsNotFraction(self, coverage):
        return np.nan

    def _coverageFraction(self, costs):
        return np.nan

    def _coverageNotFraction(self, costs):
        return np.nan


class LinearCCF(CostCovFunction):
    """
    Implementation of a linear cost-coverage-function.
    """
    def __init__(self, unit_cost, cif=True):
        """
        :param unit_cost: float representing the unit costs.
        :param cif: boolean which specifies if the coverage is treated as fraction (True) or absolute value (False).
        """
        super(LinearCCF, self).__init__(cif)
        self._setAttribute(CCFAttr.uc, unit_cost)

    def _costsFraction(self, coverage):
        return self._costsNotFraction(coverage) / 0.01

    def _costsNotFraction(self, coverage):
        return coverage * self.attr[CCFAttr.uc]

    def _coverageFraction(self, costs):
        return self._coverageNotFraction(costs) * 0.01

    def _coverageNotFraction(self, costs):
        return costs / self.attr[CCFAttr.uc]


class LogisticCCF(CostCovFunction):
    """
    Implementation of a logistic cost-coverage-function.
    """
    def __init__(self, unit_cost, saturation, pop_size, dt, cif=True):
        """
        :param unit_cost: float representing the unit costs.
        :param saturation: float representing the saturation of the function as a fraction.
        :param pop_size: float representing the population size this function is applied to.
        :param dt: float which represents the simulated time step size in years.
        :param cif: boolean which specifies if the coverage is treated as fraction (True) or absolute value (False).
        """
        super(LogisticCCF, self).__init__(cif)

        self._setAttribute(CCFAttr.uc, unit_cost)
        self._setAttribute(CCFAttr.sat, saturation)
        self._setAttribute(CCFAttr.pop, pop_size)
        self._setAttribute(CCFAttr.dt, dt)

    def _costsFraction(self, coverage):
        return - self.attr[CCFAttr.pop] * self.attr[CCFAttr.sat] * self.attr[CCFAttr.uc] / 2. \
               * np.log((self.attr[CCFAttr.sat] - coverage) / (coverage + self.attr[CCFAttr.sat]))

    def _costsNotFraction(self, coverage):
        return - self.attr[CCFAttr.pop] * self.attr[CCFAttr.sat] * self.attr[CCFAttr.uc] \
               / 2. * np.log((self.attr[CCFAttr.sat] * self.attr[CCFAttr.pop] - coverage)
                             / (coverage + self.attr[CCFAttr.sat] * self.attr[CCFAttr.pop])) \
               * 365. * self.attr[CCFAttr.dt]

    def _coverageFraction(self, costs):
        return - self.attr[CCFAttr.sat] + 2. * self.attr[CCFAttr.sat] / \
                                          (1. + np.exp(-2. * costs / (self.attr[CCFAttr.pop] * self.attr[CCFAttr.sat] * self.attr[CCFAttr.uc])))

    def _coverageNotFraction(self, costs):
        return self._coverageFraction(costs) * self.attr[CCFAttr.pop] / (365. * self.attr[CCFAttr.dt])
