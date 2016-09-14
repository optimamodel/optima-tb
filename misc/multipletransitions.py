"""
Demonstrates rationale for combining potentially incompatible transition rates.
"""

# Define things
from pylab import zeros, rand, log
npeople = 100000 # Number of agents
infrate = 0.75 # Infection rate
vacrate = 0.25 # Vaccination probability
bothallowed = False # Whether both infection and vaccination can happen

# Run simulation
people = zeros((npeople,2)) # List of agends, tracking infection and vaccination status
for p in range(npeople): # Loop over people
    doinf = True if rand()<infrate else False # See if infection happens
    dovac = True if rand()<vacrate else False # See if vaccination happens
    if doinf and dovac and not bothallowed: # Whether or not to disallow both transitions
        timetoinf = -1./infrate*log(1-rand()) # Calculate time to infection
        timetovac = -1./vacrate*log(1-rand()) # Calculate time to vaccination
        if timetoinf<timetovac: dovac = False # Infection came first
        else:                   doinf = False # Vaccination came first
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
analytic[0,0] = (1-infrate)*(1-vacrate) # People who don't move
analytic[1,0] = (  infrate)*(1-vacrate) # People who are infected
analytic[0,1] = (1-infrate)*(  vacrate) # People who are vaccinated
if bothallowed:
    analytic[1,1] = (infrate)*(vacrate) # People who are both infected and vaccinated
else:
    proportion = (infrate)*(vacrate) # People who would've been both infected and vaccinated
    denominator = (infrate+vacrate) # Normalization for rates
    analytic[1,0] += proportion * infrate/denominator # Additional people who are infected
    analytic[0,1] += proportion * vacrate/denominator # Additional people who are vaccinated
print('Analytic:\n%s\n' % analytic)

# Check how well they match up
print('Ratios:\n%s\n' % (simulation/analytic))  
    
    