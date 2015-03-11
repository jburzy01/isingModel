//
//  main.cpp
//  2D Ising Model
//
//  Created by Andrew Mascioli on 6/4/14.
//  Copyright (c) 2014 Andrew Mascioli. All rights reserved.
//

#include <iostream>
#include "2D_Lattice.h"
using namespace std;

void runWithDefaultInteraction(int argc, const char* argv[]);
void runWithCustomInteraction(int argc, const char* argv[]);
void runInteractively();

int main(int argc, const char * argv[]){
    srand(arc4random());
    if (argc == 5)     runWithDefaultInteraction(argc, argv);
    else if(argc == 6) runWithCustomInteraction(argc, argv);
    else               runInteractively();
    
    return 0;
}

//Runs the program with a default interaction constant of 1
//runs when provided four arguments in the command line
//order is numX, numY, numIterations, Temp in terms of 1/k
void runWithDefaultInteraction(int argc, const char* argv[]){
    double inputs[4];
    for(int i = 1; i < argc; i++){
        inputs[i-1] = atof(argv[i]);
    }
    Lattice l(inputs[0], inputs[1]);
    l.createInitialSpinConfiguration();
    l.runSimulation(inputs[2], inputs[3]/BoltzmannConstant);
}

//Runs the program with a custom interaction constant
//Runs when provided five arguments in the command line
//Order is interactionConst, numX, numY, numIterations, Temp in terms of 1/k
void runWithCustomInteraction(int argc, const char* argv[]){
    double inputs[5];
    for(int i = 1; i < argc; i++){
        inputs[i-1] = atof(argv[i]);
    }
    Lattice l(inputs[0], inputs[1], inputs[2]);
    l.createInitialSpinConfiguration();
    l.runSimulation(inputs[3], inputs[4]/BoltzmannConstant);
}

//Allows the program to run interactively with a terminal graphics visualization, gets called by default when there are no command line arguments or an invalid number of arguments
void runInteractively(){
    int numSitesX;
    int numSitesY;
    double constantValue;
    double temp;
    int numIterations;
    cout << "Enter the number of X Lattice Sites: ";
    cin >> numSitesX;
    cout << "Enter the number of Y Lattice Sites: ";
    cin >> numSitesY;
    cout << "Enter the Interaction Constant: ";
    cin >> constantValue;
    cout << "Enter the Temperature (In terms of 1/k): ";
    cin >> temp;
    cout << "Enter the Number of Iterations: ";
    cin >> numIterations;
    
    Lattice l(constantValue, numSitesX, numSitesY);
    l.createInitialSpinConfiguration();
    l.printSpinConfiguration();
    l.runSimulationWithVisualization(numIterations, temp/BoltzmannConstant);
}
