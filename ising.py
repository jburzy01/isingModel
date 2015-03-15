"""
Monte Carlo simulation of the Ising model
"""

import numpy as np
import copy
import random as rd
import math

probabilityPrecision = 10000; # what is this?
BoltzmannConstant = 8.617e-5;

class Ising_lattice:

	def __init__(self, xsize, ysize, zsize, I, T):
		self.xsize    = xsize  			# length of x edge
		self.ysize    = ysize  			# length of y edge
		self.zsize    = zsize  			# length of z edge

		self.initLattice() 				# a new lattice of atoms

		self.const	  = I 				# interaction constant
		self.T        = T     			# temperature
		self.M        = self.calcM()  	# total magnetization
		self.H        = self.calcH()    # total energy


    # Initializes the lattice with random spins 
	def initLattice(self):
		lattice = np.zeros((self.xsize, self.ysize, self.zsize))
		for x in range(self.xsize):
			for y in range(self.ysize):
				for z in range(self.zsize):
					if (rd.randint(1,1000) % 2 == 1):
						lattice[(x,y,z)] = 1
					else:
						lattice[(x,y,z)] = -1
		self.lattice = lattice

	def xAlignedSpins(self,x,y,z):
		return self.lattice[(x,y,z)]*self.lattice[(x-1,y,z)]

	def yAlignedSpins(self,x,y,z):
		return self.lattice[(x,y,z)]*self.lattice[(x,y-1,z)]

	def zAlignedSpins(self,x,y,z):
		return self.lattice[(x,y,z)]*self.lattice[(x,y,z-1)]

	# Calculates the Hamiltonian of the lattice
	def calcH(self):
		numAdjacentAlignedSpins = 0
		for x in range(self.xsize):
			for y in range(self.ysize):
				for z in range(self.zsize):
					numAdjacentAlignedSpins += (self.xAlignedSpins(x,y,z) + 
		      									self.yAlignedSpins(x,y,z) +
		      									self.zAlignedSpins(x,y,z))
		return -self.const*numAdjacentAlignedSpins

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
		return magnetization

	# Calculates the boltzmann probability of acceptance in the case where the flipped lattice has a higher energy
	def calculateBoltzmannProbability(self,initialHamiltonian,flippedHamiltonian):
		randNum = rd.random()%probabilityPrecision 

		# this calculation results in overflow for small T
		proportionOfBoltzmannProbabilities = math.exp(-(1/(BoltzmannConstant*self.T))*(flippedHamiltonian - initialHamiltonian))*probabilityPrecision
		if (randNum < proportionOfBoltzmannProbabilities):
			return True
		else:
			return False

		# Goes through the iterations
	def iterate(self, its, algorithm):

		if algorithm == "Metropolis":
			self.Metropolis(its)
		elif algorithm == "Wolff":
			self.Wolff(its)
		elif algorithm == "Swedsen-Wang":
			self.Swedsen-Wang(its)

	def Metropolis(self,its):

		for i in range(its):
			# choose random atom
			x = rd.randint(1,1000) % self.xsize
			y = rd.randint(1,1000) % self.ysize # what if ysize is zero?
			z = rd.randint(1,1000) % self.zsize

			testlattice = copy.deepcopy(self)
			testlattice.lattice[(x,y,z)] *= -1
			flippedEnergy = testlattice.calcH() 
			initEnergy = self.calcH()

			deltaH  = initEnergy - flippedEnergy
			      
			if (deltaH <= 0 or self.calculateBoltzmannProbability(initEnergy, flippedEnergy)):
				self.M += 2 * self.lattice[(x,y,z)]
				self.H += deltaH

				self.lattice[(x,y,z)] *= -1


x=Ising_lattice(5,5,5,1,100)

print "energy before "

print x.calcH()

x.iterate(500,"Metropolis")

print "energy after "

print x.calcH()
