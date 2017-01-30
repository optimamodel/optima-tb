

from optima_tb.programs import CostCoverageCurve

def checkCostCovCurve(curve, budget):
    print curve
    cov = curve.getValue(budget)
    print "Budget  :$%g \t\t==> Coverage : %g"%(budget,cov)
    bud =  curve.getInverseValue(cov)
    print "Coverage: %g \t==> Budget   : $%g"%(cov,bud)
    # Note that the recovered budget should be identical to the original budget
    print " Original budget : $%g"%budget
    print "Recovered budget : $%g"%bud
        
    
    
costcurve1= CostCoverageCurve()
budget = 40
checkCostCovCurve(costcurve1, budget)

params={'exponent':0.045,'saturation':3e5}
costcurve2 = CostCoverageCurve(params=params)
checkCostCovCurve(costcurve2, budget)