#include <iostream>
#include <fstream>
#include "system.h"
#include "particle.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/simplegaussianinteraction.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/interactionharmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "InitialStates/gaussiandistribution.h"
#include "Math/random.h"
#include <string>
#include <omp.h>
#include <cmath>

using namespace std;

void alphaListRun(string filename, int numberOfP);
string setMethodName(bool analyticOrNot);

int main() {

// omp_set_num_threads(4);
// int thread_num = omp_get_max_threads ();

/* The standard set-up */

    bool analyticOrNot       = true;
    bool importanceOrNot     = true;
    bool allEnergiesOrNot    = true;
    int equilibration        = 0;          // Number of the total steps used for equilibration
    int numberOfDimensions   = 2;
    int numberOfParticles    = 2;
    int numberOfSteps        = (int) pow(2.0,19.0);
    double omega             = 2.0;          // Oscillator frequency.
    double stepLength        = 0.5;          // Metropolis step length.
    int firstCriteria        = 0;            // print header in file
    double alpha             = 1.0;
    double beta              = 1.0;   
    double inititalizingStep = stepLength;

    double spinFactor  = 1;
    // int numberOfBins = 800;
    // double densityLength = 10.0;

    // interaction or spherical trap (2.82843 or 1.0)
 

    /* Set-up to run and save local energies for every step to file*/

    // allEnergiesOrNot    = false;

    // clock_t start, end;
    // // Recording the starting clock tick.
    // start = clock();

    
    // System* system = new System();
    // system->setHamiltonian                (new HarmonicOscillator(system, omega));
    // system->setWaveFunction               (new SimpleGaussian(system, alpha));

    // system->setInitialState               (new RandomUniform(system, numberOfDimensions, 
    //                                             numberOfParticles, inititalizingStep));

    // system->setEquilibration              (equilibration);

    // system->setAnalytical                 (analyticOrNot);

    // // system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
    // system->setFileName                   ("Output/test");

    // system->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
    //                                                         allEnergiesOrNot, stepLength);

    // cout << "number of steps: " << numberOfSteps << endl;

    // cout << "energy: " << system->getSampler()->getEnergy()/((double) numberOfParticles*numberOfDimensions) << endl;
   
    // end = clock();
    // double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    // cout << "CPU time: " << time_taken << " seconds" << endl;



    /* Same but for the interacting case */   


    importanceOrNot = false;
    stepLength = 0.5;
    inititalizingStep = stepLength;
    allEnergiesOrNot = false;

    spinFactor  = 1.0;
    // numberOfBins = 800;
    // densityLength = 10.0;

    // interaction or spherical trap (2.82843 or 1.0)
    beta = 1.0;    // omega_normal^2/omega_ho^2

    clock_t start, end;
    // Recording the starting clock tick.
    start = clock();


    System* system = new System();
    system->setHamiltonian                (new InteractionHarmonicOscillator(system, omega));
    system->setWaveFunction               (new SimpleGaussianInteraction(system, alpha, beta, spinFactor));

    system->setInitialState               (new RandomUniform(system, numberOfDimensions, 
                                                numberOfParticles, inititalizingStep));
    system->setEquilibration              (equilibration);
    system->setAnalytical                 (analyticOrNot);
    // system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
    system->setFileName                   ("Output/test_");

    system->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
                                                            allEnergiesOrNot, stepLength);

    cout << "energy: " << system->getSampler()->getEnergy()/((double) numberOfParticles*numberOfDimensions) << endl;

    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    cout << "CPU time: " << time_taken << " seconds" << endl;



}








string setMethodName(bool analyticOrNot){
    string methodName;

    if (analyticOrNot == true){
        methodName = "analytical_";
    }else{
        methodName = "numerical_";
    }
return methodName;
}