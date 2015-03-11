//
//  2D_Lattice.h
//  2D Ising Model
//
//  Created by Andrew Mascioli on 6/3/14.
//  Copyright (c) 2014 Andrew Mascioli. All rights reserved.
//

#ifndef ___D_Ising_Model___D_Lattice__
#define ___D_Ising_Model___D_Lattice__

#include <iostream>
#include <ctime>
#include <cmath>
#include <unistd.h>

const int probabilityPrecision = 100000000;
const double BoltzmannConstant = 8.617e-5;

const double defaultInteraction = 1;
const int UP_SPIN = 1;
const int DOWN_SPIN = -1;

class Lattice {
    
public:
    Lattice();
    Lattice(int numXSites, int numYSites);
    Lattice(double interactionConst, int numXSites, int numYSites);

    void createInitialSpinConfiguration();
    void printSpinConfiguration();
    
    void setNumLatticeSites(int inputX, int inputY);
    void setInteractionConstant(double input);
    
    Lattice flipOneSpin();
    int pickRandomXSite();
    int pickRandomYSite();
    
    bool isAccepted(double initialHamiltonian, double flippedHamiltonian, double temp);
    void acceptFlippedLattice(Lattice inputLattice);
    
    double calculateHamiltonian();
    double calculateNetMagnetization();
    
    int rightAlignedSpins(int x, int y);
    int bottomAlignedSpins(int x, int y);
    
    void runSimulation(int numIterations, double temp);
    void runSimulationWithVisualization(int numIterations, double temp);
    
    
    
private:
    constexpr static const int possibleSpinStates[2] = {DOWN_SPIN, UP_SPIN};
    int numXLatticeSites;
    int numYLatticeSites;
    double interactionConstant;
    
    int** spinConfiguration;
    
    int** copySpinConfiguration();
    bool calculateBoltzmannProbability(double initialHamiltonian, double flippedHamiltonian, double temp);
};
#endif /* defined(___D_Ising_Model___D_Lattice__) */
