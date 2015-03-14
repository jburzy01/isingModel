"""
Monte Carlo simulation of the 2D Ising model
"""

import numpy as np


class Ising_lattice:

	  def __init__(self, xsize, ysize, zsize, T):
    	self.size     = size  		# length of edge
    	self.dim 	  = dimension	# lattice dimension

    	self.const	  = I 			# interaction constant
    	self.T        = T     		# temperature
    	self.M        = 0     		# total magnetization
    	self.H        = 0     		# total energy
    
    	self.initLattice(dimension) # a new lattice of atoms
  
     	self.calcH()          		# total energy 		
    	self.calcM()          		# total magnetization

      def initLattice(self, dimension):
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



