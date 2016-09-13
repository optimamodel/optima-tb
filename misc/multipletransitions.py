"""
Demonstrates rationale for combining potentially incompatible transition rates.
"""

from pylab import zeros, rand
npeople = 1000 # Number of agents
people = zeros((npeople,2)) # List of agends, tracking infection and vaccination status
infrate = 0.5 # Infection rate
vacrate = 0.5 # Vaccination probability
bothallowed = True # Whether both infection and vaccination can happen
for p in range(npeople): # Loop over people
    infmove = rand()
    vacmove = rand()
    doinf = True if infmove<infrate else False # See if infection happens
    dovac = True if vacmove<vacrate else False # See if vaccination happens
    if not bothallowed: # Whether or not to disallow both transitions
        if doinf and dovac: # Check that both transitions are on the cards
            doinf = True if infmove>vacmove else False # Infection takes priority
            dovac = True if vacmove>infmove else False # Vaccination takes priority
    if doinf: people[p,0] = 1 # Update state
    if dovac: people[p,1] = 1 # Update state
counts = zeros((2,2)) # Count how many people are in each state
for p in range(npeople):
    i = people[p,0]
    j = people[p,1]
    counts[i,j] += 1
print(counts)
        
    
    
    