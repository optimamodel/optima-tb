"""
Demonstrates rationale for combining potentially incompatible transition rates.
"""

from pylab import zeros, rand, log
npeople = 100000 # Number of agents
infrate = 0.75 # Infection rate
vacrate = 0.25 # Vaccination probability
bothallowed = False # Whether both infection and vaccination can happen

z = [0,0]
people = zeros((npeople,2)) # List of agends, tracking infection and vaccination status
for p in range(npeople): # Loop over people
    infmove = rand()
    vacmove = rand()
    doinf = True if infmove<infrate else False # See if infection happens
    dovac = True if vacmove<vacrate else False # See if vaccination happens
    if doinf and dovac and not bothallowed: # Whether or not to disallow both transitions
        timetoinf = -1./infrate*log(1-rand())
        timetovac = -1./vacrate*log(1-rand())
        if timetoinf<timetovac: 
            dovac = False
            z[0] += 1
        else: 
            doinf = False
            z[1] += 1
    if doinf: people[p,0] = 1 # Update state
    if dovac: people[p,1] = 1 # Update state

fractions = zeros((2,2)) # Count how many people are in each state
for p in range(npeople):
    i = people[p,0]
    j = people[p,1]
    fractions[i,j] += 1./npeople
print(fractions)

analytic = zeros((2,2)) # Calculate equivalent analytic result
analytic[0,0] = (1-infrate)*(1-vacrate)
analytic[1,0] = (  infrate)*(1-vacrate)
analytic[0,1] = (1-infrate)*(  vacrate)
if bothallowed:
    analytic[1,1] = (infrate)*(vacrate)
else:
    proportion = (infrate)*(vacrate)
    denominator = (infrate+vacrate)
    analytic[1,0] += proportion * infrate/denominator
    analytic[0,1] += proportion * vacrate/denominator
print(analytic)

print(fractions/analytic)
    
    
    