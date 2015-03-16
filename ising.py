"""
Monte Carlo simulation of the Ising model
"""

import numpy as np
import copy
import random as rd
import math

probabilityPrecision = 100000000; # what is this?
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
		randNum = rd.randint(1,1000)%probabilityPrecision # what range should randint be in ?

		# this calculation results in overflow for small T
		proportionOfBoltzmannProbabilities = probabilityPrecision*math.exp(-(1/(BoltzmannConstant*self.T))*(flippedHamiltonian - initialHamiltonian))
		if (randNum < proportionOfBoltzmannProbabilities):
			return True
		else:
			return False

		# Goes through the iterations
	def iterate(self, its, algorithm):

		for i in range(its):
			# choose random 
			x = rd.randint(1,1000) % self.xsize
			if self.ysize == 0:
				y = 0
			else:
				y = rd.randint(1,1000) % self.ysize

			if self.zsize == 0:
				z = 0
			else:
				z = rd.randint(1,1000) % self.zsize

			p = (x,y,z)

			if algorithm == "Metropolis":
				self.Metropolis(p)
			elif algorithm == "Wolff":
				self.Wolff(p)
			elif algorithm == "Swedsen-Wang":
				self.Swedsen_Wang(p)

	def Metropolis(self,randPoint):

		testlattice = copy.deepcopy(self)
		testlattice.lattice[randPoint] *= -1
		flippedEnergy = testlattice.calcH() 
		initEnergy = self.calcH()
     
		if (flippedEnergy <= initEnergy or self.calculateBoltzmannProbability(initEnergy, flippedEnergy)):
			self.lattice[randPoint] *= -1

	def Wolff(self,randPoint):

		# define the cluster and perimeter list
		cluster = []
		perimeterList = []

		# choose random spin to be seed of cluster
		cluster.append(randPoint)

		# examine all neighbors
		# add parallel neighbors to perimeter list
		neighbors = self.get_Neighbors(randPoint)
		self.addToPerimeter(perimeterList, randPoint, neighbors)

		while not perimeterList:
			# pop spin from perimeter list
			neighbor = perimeterList.pop()

			if not self.addToCluster(cluster,randPoint,neighbors):
				perimeterList.append(neighbor)
			else:
				self.inspectNeighbors()

		self.flipCluster(cluster)

		# for each of its neighbors that is already in the cluster,
		# 	add spin to cluster w/ prob = 1 - exp(-2*beta*J)

		# returns a list of neighbors to a given point
	def get_Neighbors(self,point):
		x = point[0]
		y = point[1]
		z = point[2]

		neighbors = []

		neighbors.append((x-1,y,z))
		neighbors.append((x,y-1,z))
		neighbors.append((x,y,z-1))

		if (x == self.xsize):
			neighbors.append((0,y,z))
		else:
			neighbors.append((x+1,y,z))
		if (y == self.ysize):
			neighbors.append((x,0,z))
		else:
			neighbors.append((x,y+1,z))
		if (z == self.zsize):
			neighbors.append((x,y,0))
		else:
			neighbors.append((x,y,z+1))

		return neighbors

	def addToPerimeter(self, perimeterList, randPoint, neighbors):

		spin = self.lattice[point]

		for i in range(len(neighbors)):
			neighbor = neighbors[i]
			if spin == self.lattice[neighbor]:
				perimeterList.append(neighbor)

	# unsure about this function
	# how do we "add to the cluster with prob = 1 - exp(-2*beta*J)"?
	def addToCluster(self, cluster, point, neighbors):
		# for each of its neighbors that is already in the cluster,
		# 	add spin to cluster w/ prob = 1 - exp(-2*beta*J)

		probability = probabilityPrecision*(1-math.exp(-(2/(BoltzmannConstant*self.T))*self.const))

		for i in range(len(cluster)):
			randNum = rd.randint(1,1000)%probabilityPrecision
			if randNum < probability:
				cluster.append(point)
				return True

		return False

	def inspectNeighbors(self,cluster, perimeterList, point, neighbors):
		for i in range(len(neighbors)):
			neighbor = neighbors[i]
			if spin == self.lattice[neighbor]:
				if not (neighbor in cluster) or not (neighbor in perimeterList):
					perimeterList.append(neighbor)

	# flips all spins in the cluster
	def flipCluster(self,cluster):
		for i in range(len(cluster)):
			point = cluster[i]
			self.lattice[point] *= -1



#	def Swedsen_Wang(self,randPoint):


x=Ising_lattice(10,10,10,1,1000)

print "energy before "

print x.calcH()

x.iterate(50,"Wolff")

print "energy after "

print x.calcH()
