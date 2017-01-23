

from optima_tb.programs import CostCoverageCurve

ccc= CostCoverageCurve()
budget = 40
cov = ccc.getValue(budget)
print "Budget  :$%g \t\t==> Coverage : %g"%(budget,cov)
bud =  ccc.getInverseValue(cov)
print "Coverage: %g \t==> Budget   : $%g"%(cov,bud)
# Note that the recovered budget should be identical to the original budget
print " Original budget : $%g"%budget
print "Recovered budget : $%g"%bud

