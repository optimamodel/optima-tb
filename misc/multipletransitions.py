"""
Demonstrates rationale for combining potentially incompatible transition rates.
"""

from pylab import zeros, rand, log, exp
npeople = 100000 # Number of agents
infrate = 0.5 # Infection rate
vacrate = 0.5 # Vaccination probability
timestep = exp(1)/2 # The timestep being considered
bothallowed = True # Whether both infection and vaccination can happen

z = [0,0]
people = zeros((npeople,2)) # List of agends, tracking infection and vaccination status
inftimes = [] # Store times to infection for e.g. mean(array(inftimes<timestep))
vactimes = []
for p in range(npeople): # Loop over people
    inftime = -1./infrate*log(1-rand())
    vactime = -1./vacrate*log(1-rand())
    inftimes.append(inftime)
    vactimes.append(vactime)
    doinf = True if inftime<timestep else False # See if infection happens
    dovac = True if vactime<timestep else False # See if vaccination happens
    if not bothallowed: # Whether or not to disallow both transitions
        if doinf and dovac: # Check that both transitions are on the cards
            if inftime<vactime: dovac = False # Infection takes priority
            if vactime<inftime: doinf = False # Vaccination takes priority
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
    analytic[1,0] += (infrate)*(vacrate) * (infrate/(infrate+vacrate))
    analytic[0,1] += (infrate)*(vacrate) * (vacrate/(infrate+vacrate))
print(analytic)

print(fractions/analytic)
    
    
    