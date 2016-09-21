"""
Demonstrates rationale for combining potentially incompatible transition rates.
"""

# Define things
from pylab import zeros, rand, log, exp
npeople = 1000000 # Number of agents
infrate = 2.0 # Infection rate (average number of events per year, can be over 1)
vacrate = 0.5 # Vaccination rate (average number of events per year, can be over 1)
bothallowed = False # Whether both infection and vaccination can happen
timestep = 0.5 # Timestep for which probability of transition is calculated 
infprob = 1-exp(-infrate*timestep)
vacprob = 1-exp(-vacrate*timestep)

# Run simulation
people = zeros((npeople,2)) # List of agends, tracking infection and vaccination status
for p in range(npeople): # Loop over people
    timetoinf = -1./infrate*log(1-rand()) # Calculate time to infection
    timetovac = -1./vacrate*log(1-rand()) # Calculate time to vaccination

    doinf = False
    dovac = False
    if timetoinf < timestep: doinf = True
    if timetovac < timestep: dovac = True
    if not bothallowed:
        if timetoinf > timetovac: doinf = False
        else: dovac = False
    
#    doinf = True if rand()<infrate else False # See if infection happens
#    dovac = True if rand()<vacrate else False # See if vaccination happens
#    if doinf and dovac and not bothallowed: # Whether or not to disallow both transitions
#        timetoinf = -1./infrate*log(1-rand()) # Calculate time to infection
#        timetovac = -1./vacrate*log(1-rand()) # Calculate time to vaccination
#        if timetoinf<timetovac: dovac = False # Infection came first
#        else:                   doinf = False # Vaccination came first
    if doinf: people[p,0] = 1 # Update state
    if dovac: people[p,1] = 1 # Update state

# Count how many people are in each state
simulation = zeros((2,2)) 
for p in range(npeople):
    i = people[p,0]
    j = people[p,1]
    simulation[i,j] += 1./npeople
print('Simulation:\n%s\n' % simulation)

# Calculate equivalent analytic result
analytic = zeros((2,2)) 
analytic[0,0] = (1-infprob)*(1-vacprob) # People who don't move
analytic[1,0] = (  infprob)*(1-vacprob) # People who are infected
analytic[0,1] = (1-infprob)*(  vacprob) # People who are vaccinated
if bothallowed:
    analytic[1,1] = (infprob)*(vacprob) # People who are both infected and vaccinated
else:
    proportion = infprob + vacprob - infprob*vacprob # People who are either infected or vaccinated... or both
    infsplit = infrate
    vacsplit = vacrate
    denominator = (infsplit+vacsplit) # Normalization for rates
    analytic[1,0] = proportion * infsplit/denominator # People who are infected
    analytic[0,1] = proportion * vacsplit/denominator # People who are vaccinated
print('Analytic:\n%s\n' % analytic)

# Check how well they match up
print('Ratios:\n%s\n' % (simulation/analytic))  
    
    