"""
Monte Carlo simulation of the Ising model
"""

import numpy as np
import copy

class Ising_lattice:

	  def __init__(self, xsize, ysize, zsize, I, T):
    	self.xsize     = xsize  	# length of x edge
    	self.ysize     = ysize  	# length of y edge
    	self.zsize     = zsize  	# length of z edge

    	self.const	  = I 			# interaction constant
    	self.T        = T     		# temperature
    	self.M        = 0     		# total magnetization
    	self.H        = 0     		# total energy
    
    	self.initLattice() # a new lattice of atoms
  
     	self.calcH()          		# total energy 		
    	self.calcM()          		# total magnetization

      def initLattice(self):
    	lattice = np.zeros((self.xsize, self.ysize, self.zsize))
	    for x in range(self.xsize):
	      	for y in range(self.ysize):
	      		for z in range(self.zsize):
		        	if (random() < p):
		          		lattice[(x,y,z)] = 1
		        	else:
		          		lattice[(x,y,z)] = -1	    
    	self.lattice = lattice

	    def calcH(self):
	    	numAdjacentAlignedSpins = 0
		    for x in range(self.xsize):
		      	for y in range(self.ysize):
		      		for z in range(self.zsize):
		      			numAdjacentAlignedSpins += (xAlignedSpins(x,y,z) + 
		      										yAlignedSpins(x,y,z) +
		      										zAlignedSpins(x,y,z))
		    self.H = -self.const*numAdjacentAlignedSpins

		def xAlignedSpins(x,y,z):
			return self.lattice(x,y,z)*self.lattice[(x-1,y,z)]

		def yAlignedSpins(x,y,z):
			return self.lattice(x,y,z)*self.lattice[(x,y-1,z)]

		def zAlignedSpins(x,y,z):
			return self.lattice(x,y,z)*self.lattice[(x,y,z-1)]

	    def calcM(self):
	    	magnetization = 0
		    for x in range(self.xsize):
		      	for y in range(self.ysize):
		      		for z in range(self.zsize):
		      			if self.lattice[(x,y,z)] == 1:
		      				magnetization += 1
		      			else:
		      				magnetization += -1
		    self.M = magnetization

		  def iterate(self, its):
		  	for i in range(its):
		      	# choose random atom
		    	x = randrange(self.xsize)
		    	y = randrange(self.ysize)
		    	z = randrange(self.zsize)

		      	testlattice = copy.deepcopy(self)
		      	testlattice.lattice[(x,y,z)] *= -1
		      	newEnergy = testlattice.calcH()

		      	deltaH  = self.H - newEnergy
		      
		      	if (deltaH <= 0):
		        	self.lattice[(x,y,z)] *= -1

		        	self.M += 2 * self.lattice[(x,y)]
		        	self.H += deltaH

		    else:
		    	# progress depending on algorithm


