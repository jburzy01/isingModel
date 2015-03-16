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
		randNum = rd.uniform(0,probabilityPrecision)

		proportionOfBoltzmannProbabilities = probabilityPrecision*math.exp(-(1/(BoltzmannConstant*self.T))*(flippedHamiltonian - initialHamiltonian))
		if (randNum < proportionOfBoltzmannProbabilities):
			return True
		else:
			return False

		# Goes through the iterations
	def iterate(self, its, algorithm):

		for i in range(its):
			# choose random 
			x = rd.randint(0,self.xsize-1)
			if self.ysize == 0:
				y = 0
			else:
				y = rd.randint(0,self.ysize-1)

			if self.zsize == 0:
				z = 0
			else:
				z = rd.randint(0,self.zsize-1) 

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

		# iterate through the perimeter list
		while perimeterList:

			# pop spin from perimeter list
			point = perimeterList.pop()
			pointsNeighbors = self.get_Neighbors(point)
			# attempt to add the spin to the cluster
			# if spin is not added to cluster, move to next perimeter spin
			if not self.addToCluster(cluster,point,pointsNeighbors):
				continue
			# if spin is added to cluster, inspect neighbors with parrallel spins
			# and add them to the perimeter list
			else:
				self.inspectNeighbors(cluster,perimeterList,point,pointsNeighbors)

		# flip all spins in the cluster
		self.flipCluster(cluster)

	# returns a list of neighbors to a given point
	# doesn't work for 1D and 2D cases
	def get_Neighbors(self,point):
		x = point[0]
		y = point[1]
		z = point[2]

		neighbors = []

		if (x == 0):
			neighbors.append((self.xsize-1,y,z))
		else:
			neighbors.append((x-1,y,z))
		if (x == self.xsize-1):
			neighbors.append((0,y,z))
		else:
			neighbors.append((x+1,y,z))

		if self.ysize > 0:
			if (y == 0):
				neighbors.append((x,self.ysize-1,z))
			else:
				neighbors.append((x,y-1,z))
			if (y == self.ysize-1):
				neighbors.append((x,0,z))
			else:
				neighbors.append((x,y+1,z))

		if self.zsize > 0:
			if (z == 0):
				neighbors.append((x,y,self.zsize-1))
			else:
				neighbors.append((x,y,z-1))
			if (z == self.zsize-1):
				neighbors.append((x,y,0))
			else:
				neighbors.append((x,y,z+1))

		return neighbors

	# adds neighbors with parallel spins to the perimeter list
	def addToPerimeter(self, perimeterList, point, neighbors):

		spin = self.lattice[point]

		for i in range(len(neighbors)):
			neighbor = neighbors[i]
			if spin == self.lattice[neighbor]:
				perimeterList.append(neighbor)

	# this function doesn't work properly
	def addToCluster(self, cluster, point, neighbors):
		# for each of its neighbors that is already in the cluster,
		# 	add spin to cluster w/ prob = 1 - exp(-2*beta*J)

		probability = probabilityPrecision*(1-math.exp(-(2/(BoltzmannConstant*self.T))*self.const))

		for i in range(len(neighbors)):
			if neighbors[i] in cluster:

				randNum = rd.uniform(0,probabilityPrecision)

				if randNum < probability:
					cluster.append(point)
					return True

		return False

	# if a parallel neighboring spin is neither in the cluster nor in the 
	# perimeter list, we add it to the perimeter list
	def inspectNeighbors(self,cluster, perimeterList, point, neighbors):

		spin = self.lattice[point]

		for i in range(len(neighbors)):
			neighbor = neighbors[i]

			if spin == self.lattice[neighbor]:
				if not (neighbor in cluster) and not (neighbor in perimeterList):
					perimeterList.append(neighbor)

	# flips all spins in the cluster
	def flipCluster(self,cluster):
		for i in range(len(cluster)):
			point = cluster[i]
			self.lattice[point] *= -1


#	def Swedsen_Wang(self,randPoint):


x=Ising_lattice(10,10,10,1,1000000)

print "energy before "

print x.calcH()

x.iterate(2,"Wolff")

print "energy after "

print x.calcH()
