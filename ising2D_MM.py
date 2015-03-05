#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
#
# Code to compute average magnetization of 2D Ising magnet with no applied field
#
# The Onsager analytical solution gives Tc=2.3J/k
# The Mean-field model predicts Tc=4J/k
#

import scipy
import numpy as np
import math

N     = 10            # size of 1d lattice.  total lattice size will be NxN
nIter = 10000         # number of monte carlo iterations to perform
T     = 10            # temperature in units of J/kB
deltaItt = 100           # how often to write output
invT = 1.0/float(T)

#################################################################################################
###############################            SUBROUTINES          #################################
#################################################################################################

# compute the total energy of a lattice configuration for the Ising model with zero applied field
def TotalEnergy(lattice):
	totalEnergy=0
	for i in range(len(lattice)):
		for j in range(len(lattice)):
			S = lattice[i,j]
			# sum spins of nearest neighbors (mod N comes from applying periodic boundary conditions)
			WF = lattice[(i+1)%N, j] + lattice[i,(j+1)%N] + lattice[(i-1)%N,j] + lattice[i,(j-1)%N]
			totalEnergy += -WF*S # Each neighbor gives energy 1.0
	return totalEnergy/2. # Each par counted twice

# generate a random configuation of square 2D lattice with binary occupation
def RandomL(N):
	lattice = np.zeros((N,N), dtype=int)
	for i in range(N):
		for j in range(N):
			lattice[i,j] = np.sign(2*np.random.rand()-1)
	return lattice


#################################################################################################
###############################            MAIN PROGRAM         #################################
#################################################################################################


# We need to define an initial configuration.  In this case we will do so by randomly assigning spins in the lattice
lattice = RandomL(N)
#print lattice
# initial lattice energy
energy = TotalEnergy(lattice)
#print energy
N2 = N*N

#HUH ?
#wT = np.linspace(4,0.5,100)
#print wT

# Measurements
Mn = sum(sum(lattice))
print "Initial magnetism:", Mn
avgE = 0
avgMn = 0
count = 0
# Run Monte Carlo iterations
for itt in range (0, nIter):
	# pick a random integer between 1 and NxN
	t = int(N2*np.random.rand())
	# determine x and y values of random integer
	i = t%N
	j = int(t/N)
#	print i, j
	# determine spin of randomly selected particle
	S = lattice[i,j]
        # spins of  nearest neighbors 
	WF = lattice[(i+1)%N, j] + lattice[i,(j+1)%N] + lattice[(i-1+N)%N,j] + lattice[i,(j-1+N)%N]
	deltaE = 2*S*WF
	rand = np.random.rand()
#	print "deltaE (units of J):", deltaE, "metropolis term:", math.exp(-float(deltaE)*invT), "rand:", rand
	# check to see if we should accept the move
	if deltaE<0 or math.exp(-float(deltaE)*invT)>rand: # accept move based on metropolis criteria
		# change spin
		lattice[i,j] = -S;
		# modify total energy
 	        energy += deltaE
#		print "New Energy:", energy, "Check=", TotalEnergy(lattice)
		# update magnetism
		Mn -= 2*S;
#		print "New magnetism:", Mn
	if itt%deltaItt==0:
#		print "Step:", itt
#		print lattice
		avgE += energy
		avgMn += Mn
		count += 1
		print itt, float(avgE)/float(count), float(avgMn)/float(count)



# Finish averages
	avgE = float(avgE)/float(count)
	avgMn = float(avgMn)/float(count)
