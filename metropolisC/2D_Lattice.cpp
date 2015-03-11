//
//  2D_Lattice.cpp
//  2D Ising Model
//
//  Created by Andrew Mascioli on 6/3/14.
//  Copyright (c) 2014 Andrew Mascioli. All rights reserved.
//

#include "2D_Lattice.h"

//Default constructor creates a 10x10 lattice, with an interaction constant of 0
//Not used for the most part
Lattice::Lattice(){
    interactionConstant= 0;
    setNumLatticeSites(10, 10);
    spinConfiguration = new int*[numXLatticeSites];
    for(int i = 0; i < numXLatticeSites;i++){
        spinConfiguration[i] = new int[numYLatticeSites];
    }
}

//Constructor that uses the default interaction of 1 to create the lattice
Lattice::Lattice(int numXSites, int numYSites){
    setInteractionConstant(defaultInteraction);
    setNumLatticeSites(numXSites, numYSites);
    spinConfiguration = new int*[numXLatticeSites];
    for(int i = 0; i < numXLatticeSites;i++){
        spinConfiguration[i] = new int[numYLatticeSites];
    }
}

//Constructor that takes in the size of the lattice, as well as the interaction constant
Lattice::Lattice(double interactionConst, int numXSites, int numYSites){
    setInteractionConstant(interactionConst);
    setNumLatticeSites(numXSites, numYSites);
    spinConfiguration = new int*[numXLatticeSites];
    for(int i = 0; i < numXLatticeSites;i++){
        spinConfiguration[i] = new int[numYLatticeSites];
    }
}

//Simply sets the iteraction constant, kind of useless
void Lattice::setInteractionConstant(double input){
    interactionConstant = input;
}

//Sets the number of lattice sites for a lattice, not really used unless changing the lattice size
void Lattice::setNumLatticeSites(int inputX, int inputY){
    numXLatticeSites = inputX;
    numYLatticeSites = inputY;
}

//Goes throughout the entire 2D lattice and assigns a random spin to it, either +1 or -1
void Lattice::createInitialSpinConfiguration(){
    for(int i = 0; i < numXLatticeSites; i++){
        for(int j = 0; j < numYLatticeSites; j++){
            if(rand() % 2 == 1)
                spinConfiguration[i][j] = UP_SPIN;
            else
                spinConfiguration[i][j] = DOWN_SPIN;
        }
    }
}

//Prints out the spin config in a two dimensional way, ANSI escape codes are used to clear the terminal window
//Note that the clearing of the terminal only works on Unix based terminals as far as I know
//It uses ones for the positive spin and blanks for the negative spins, in order to more easily view clustering
void Lattice::printSpinConfiguration(){
    printf("\033[2J\033[1;1H");
    for(int i = 0; i < numXLatticeSites; i++){
        std::cout << "[ ";
        for(int j = 0; j < numYLatticeSites; j++){
            if(j >= numYLatticeSites){
                if(spinConfiguration[i][j] == UP_SPIN)
                    std::cout << " 1";
                else
         //           std::cout << spinConfiguration[i][j];
                    std::cout << "  ";
                std::cout << " ";
            }
            else{
                if(spinConfiguration[i][j] == UP_SPIN)
                    std::cout << " 1";
                else
          //          std::cout << spinConfiguration[i][j];
                    std::cout << "  ";
                std::cout << " ";
            }
        }
        std::cout << "]" << std::endl;
    }
    std::cout << calculateNetMagnetization() << std::endl;
}

//Simply copies the spin configuration from one lattice and returns a pointer to where the new one is
//Necessary due to the use of dynamic arrays, if not used any changes to the array would modify the original
int** Lattice::copySpinConfiguration(){
    int** copiedConfig = new int*[numXLatticeSites];
    for(int i = 0; i < numXLatticeSites;i++){
        copiedConfig[i] = new int[numYLatticeSites];
    }
    for(int i = 0; i < numXLatticeSites; i++){
        for(int j = 0; j < numYLatticeSites; j++){
            copiedConfig[i][j] = spinConfiguration[i][j];
        }
    }
    return copiedConfig;
}

//Generates two random numbers to be used as indices and flips the spin of that index, will return the
//lattice with the flipped spin as a new lattice
Lattice Lattice::flipOneSpin(){
    Lattice newConfig(interactionConstant, numXLatticeSites, numYLatticeSites);
    newConfig.spinConfiguration = copySpinConfiguration();
    
    int indexX = pickRandomXSite();
    int indexY = pickRandomYSite();
    
    if(newConfig.spinConfiguration[indexX][indexY] == UP_SPIN)
        newConfig.spinConfiguration[indexX][indexY] = DOWN_SPIN;
    else
        newConfig.spinConfiguration[indexX][indexY] = UP_SPIN;
    
    return newConfig;
}

//Calculates the energy level of the current lattice and returns it, for each lattice site,
//only takes into account the right and bottom neighbors, as the periodic boundary conditions
//will cuase the left and top to be taken care of as the calculation starts in the top left and finishes in the bottom right
double Lattice::calculateHamiltonian(){
    int numAdjacentAlignedSpins = 0;
    for(int i = 0; i < numXLatticeSites; i++){
        for(int j = 0; j < numYLatticeSites; j++){
            numAdjacentAlignedSpins += (rightAlignedSpins(i, j) + bottomAlignedSpins(i, j));
        }
    }
    return -interactionConstant*numAdjacentAlignedSpins;
}


//Calculates the net magnitization of the lattice by iterating through each point and adding in the spin configuartion
//at each point
double Lattice::calculateNetMagnetization(){
    int netMagnetization = 0;
    for(int i = 0; i < numXLatticeSites; i++){
        for(int j = 0; j < numYLatticeSites; j++){
            if(spinConfiguration[i][j] == UP_SPIN)
                netMagnetization += UP_SPIN;
            else
                netMagnetization += DOWN_SPIN;
        }
    }
    return netMagnetization;
}

//returns positive if the site and its right neight have aligned spins, negative otherwise, deals with the boundary conditions
int Lattice::rightAlignedSpins(int x, int y){
    if((x+1) >= numXLatticeSites)
        return spinConfiguration[x][y]*spinConfiguration[0][y];
    else
        return spinConfiguration[x][y]*spinConfiguration[x+1][y];
}

//same as above except with the bottom neighbor
int Lattice::bottomAlignedSpins(int x, int y){
    if((y+1) >= numYLatticeSites)
        return spinConfiguration[x][y]*spinConfiguration[x][0];
    else
        return spinConfiguration[x][y]*spinConfiguration[x][y+1];
}

//Generates a random lattice site for the x direction
int Lattice::pickRandomXSite(){
    return (int) rand() % numXLatticeSites;
}

//Generates a random lattice site for the y direction
int Lattice::pickRandomYSite(){
    return (int) rand() % numYLatticeSites;
}

//assuming that the conditions for acceptance are met, this will delete the current spin config and accept the flipped one
void Lattice::acceptFlippedLattice(Lattice inputLattice){
    delete [] spinConfiguration;
    spinConfiguration = inputLattice.spinConfiguration;
}

//Calculates the boltzmann probability of acceptance in the case where the flipped lattice has a higher energy
bool Lattice::calculateBoltzmannProbability(double initialHamiltonian, double flippedHamiltonian, double temp){
    double randNum = (double)(rand()%probabilityPrecision);
    double proportionOfBoltzmannProbabilities = exp(-(1/(BoltzmannConstant*temp))*(flippedHamiltonian - initialHamiltonian))*probabilityPrecision;
    if(randNum < proportionOfBoltzmannProbabilities)
        return true;
    else
        return false;
}

//Will determine whether the flipped lattice should be accepted, first checks if the flipped energy is less, if so accept it, otherwise calculate the boltzmann probability and try to accept it based on that, otherwise reject it
bool Lattice::isAccepted(double initialHamiltonian, double flippedHamiltonian, double temp){
    if(flippedHamiltonian <= initialHamiltonian){
        return true;
    }
    else if(calculateBoltzmannProbability(initialHamiltonian, flippedHamiltonian, temp)){
        return true;
    }
    else
        return false;
}

//Runs the simulation without any visualization, useful for obtaining results from the command line and scripting, will return the net magnetization with each iteration
void Lattice::runSimulation(int numIterations, double temp){
    for(int i = 0; i < numIterations; i++){
        double currentHamiltonian = calculateHamiltonian();
        
        Lattice flippedLattice = flipOneSpin();
        double flippedHamiltonian = flippedLattice.calculateHamiltonian();
        
        if(isAccepted(currentHamiltonian, flippedHamiltonian, temp))
            acceptFlippedLattice(flippedLattice);
        std::cout << calculateNetMagnetization() << std::endl;
    }
}

//This function actually goes through running the simulation, at each iteration will flip one site and try to accept it and then print out the configuration after the attempted acceptance, it will continue this for numIterations iterations
void Lattice::runSimulationWithVisualization(int numIterations, double temp){    
    for(int i = 0; i < numIterations; i++){
        double currentHamiltonian = calculateHamiltonian();
        
        Lattice flippedLattice = flipOneSpin();
        double flippedHamiltonian = flippedLattice.calculateHamiltonian();
        
        if(isAccepted(currentHamiltonian, flippedHamiltonian, temp))
            acceptFlippedLattice(flippedLattice);
        
        std::cout << std::endl;
        printSpinConfiguration();
        usleep(27500);
        
    }        

}








