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

    # Initializes the lattice with random spins 
	def initLattice(self):
		lattice = np.zeros((self.xsize, self.ysize, self.zsize))
		for x in range(self.xsize):
			for y in range(self.ysize):
				for z in range(self.zsize):
					if (random() % 2 == 1):
						lattice[(x,y,z)] = 1
					else:
						lattice[(x,y,z)] = -1
		self.lattice = lattice

		# Calculates the Hamiltonian of the lattice
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

	# Calculates the net magnetization of the lattice
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

	# Calculates the boltzmann probability of acceptance in the case where the flipped lattice has a higher energy
	def calculateBoltzmannProbability(initialHamiltonian,flippedHamiltonian):
		randNum = random()%probabilityPrecision
		proportionOfBoltzmannProbabilities = math.exp(-(1/(BoltzmannConstant*self.T))*(flippedHamiltonian - initialHamiltonian))*probabilityPrecision
		if (randNum < proportionOfBoltzmannProbabilities):
			return true
		else:
			return false

		# Goes through the iterations
	def iterate(self, its, algorithm):
		for i in range(its):
			# choose random atom
			x = randrange(self.xsize)
			y = randrange(self.ysize)
			z = randrange(self.zsize)

			if algorithm == "Metropolis":
				self.Metropolis(x,y,z)
			elif algorithm == "Wolff":
				self.Wolff()
			elif algorithm == "Swedsen-Wang":
				self.Swedsen-Wang()


	def Metropolis(self):
		testlattice = copy.deepcopy(self)
		testlattice.lattice[(x,y,z)] *= -1
		flippedEnergy = testlattice.calcH()

		deltaH  = self.H - newEnergy
			      
		if (deltaH <= 0 or calculateBoltzmannProbability(self.H, newEnergy)):
			self.lattice[(x,y,z)] *= -1

			self.M += 2 * self.lattice[(x,y)]
			self.H += deltaH


