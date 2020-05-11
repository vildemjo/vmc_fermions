#include <iostream>
#include <fstream>
#include "system.h"
#include "particle.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/slaterdeterminant.h"
#include "WaveFunctions/slaterdeterminantinteraction.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/interactionharmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "InitialStates/gaussiandistribution.h"
#include "Math/random.h"
#include <string>
// #include <mpi.h>
#include <cmath>

using namespace std;

void alphaListRun(string filename, int numberOfP);
string setMethodName(bool analyticOrNot);

int main() {
/* The standard set-up */

    bool analyticOrNot       = true;
    bool importanceOrNot     = true;
    bool allEnergiesOrNot    = true;
    int equilibration        = 1e5;          // Number of the total steps used for equilibration
    int numberOfDimensions   = 2;
    int numberOfParticles    = 6;
    int numberOfSteps        = (int) pow(2.0,20.0);
    double omega             = 1.0;          // Oscillator frequency.
    double stepLength        = 0.5;          // Metropolis step length.
    int firstCriteria        = 0;            // print header in file
    double alpha             = 1.0;
    // double beta              = 0;   
    double inititalizingStep = stepLength;
    // int my_rank = 0;

    // double spinFactor  = 1;
    // int numberOfBins = 400;
    // double densityLength = 5.0;

    // interaction or spherical trap (2.82843 or 1.0)
 

   /* Set-up to run and save local energies for every step to file*/

    allEnergiesOrNot    = true;
    importanceOrNot     = false;

    clock_t start, end;
    // Recording the starting clock tick.
    start = clock();

    
    System* system = new System();
    system->setHamiltonian                (new HarmonicOscillator(system, omega));
    system->setWaveFunction               (new SlaterDeterminant(system, alpha));

    system->setInitialState               (new RandomUniform(system, numberOfDimensions, 
                                                numberOfParticles, inititalizingStep));

    system->setEquilibration              (equilibration);

    system->setAnalytical                 (analyticOrNot);

    // system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
    system->setFileName                   ("Output/slater_test");

    system->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
                                                            allEnergiesOrNot, stepLength);

    cout << "number of steps: " << numberOfSteps << endl;

    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    cout << "CPU time: " << time_taken << " seconds" << endl;



    /* Same but for the interacting case */   


    // importanceOrNot = true;
    // stepLength = 0.05;
    // inititalizingStep = stepLength;
    // allEnergiesOrNot = true;

    // spinFactor  = 1.0;
    // string file_name = "Output/test_interaction_" + to_string(my_rank);

    // // interaction or spherical trap (2.82843 or 1.0)
    // beta = 0.380869;    // omega_normal^2/omega_ho^2
    // alpha = 0.994229;

    // clock_t start2, end2;

    // // Recording the starting clock tick.
    // start2 = clock();


    // System* system2 = new System();
    // system2->setHamiltonian                (new InteractionHarmonicOscillator(system2, omega));
    // system2->setWaveFunction               (new SlaterDeterminantInteraction(system2, alpha, beta, spinFactor));

    // system2->setInitialState               (new GaussianDistribution(system2, numberOfDimensions, 
    //                                             numberOfParticles, inititalizingStep));
    // system2->setEquilibration              (equilibration);
    // system2->setAnalytical                 (analyticOrNot);
    // system2->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
    // system2->setFileName                   (file_name);

    // system2->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
    //                                                         allEnergiesOrNot, stepLength);

    // // cout << "energy: " << system->getSampler()->getEnergy() << endl;

    // end2 = clock();
    // double time_taken2 = double(end2 - start2) / double(CLOCKS_PER_SEC); 
    // cout << "CPU time: " << time_taken2 << " seconds" << endl;

    // MPI_Finilize();

// ---------------------------------------------------------------------------------------------------------------------

/* Gradient descent */

    // ofstream  file2, file;

    // double energyChange = 1.0;
    // double stopCriteria = 1e-6;
    // double energyNew = 0.0;
    // std::vector<double> energyDerivative(2,0);
    // double alphaNew = 0;
    // double betaNew = 0;
    // double minimizationRate_alpha = 0.01;
    // double minimizationRate_beta = 0.01;
    // allEnergiesOrNot = false;
    // importanceOrNot = true;
    // alpha = 1.0;
    // beta = 0.3;
    // stepLength = 0.05;
    // inititalizingStep = stepLength;

    // double energy       = 0;

    // numberOfDimensions  = 2;
    // numberOfParticles   = 2;
    // numberOfSteps       = (int) std::pow(2,19.0);

    // string file_name = "Output/test_gradient_descent_imp.txt";
    // string energy_file = "Output/test_energy_descent_imp.txt";

    // file.open (file_name, ios::out | ios::trunc);
    // file << "Alpha: \t Beta: \t Energy: \t Derivative (alpha): \t Derivative (beta): \n";
    // file.close();
    // cout << "Alpha: \t Beta: \t Energy: \t Derivative (alpha): \t Derivative (beta): \n";

    // file2.open (energy_file, ios::out | ios::trunc);
    // file2 << "Alpha: \t Beta: \t Energy: \t Kinetic Energy: \t Potential Energy:  \t Interaction Energy:\n";
    // file2.close();

    // for (int k=0;  energyChange > stopCriteria; k++){
    

    //     System* system = new System();
    //     system->setHamiltonian              (new InteractionHarmonicOscillator(system, omega));
    //     system->setWaveFunction             (new SlaterDeterminantInteraction(system, alpha, beta, spinFactor));
    //     system->setInitialState             (new GaussianDistribution(system, numberOfDimensions, 
    //                                                 numberOfParticles, inititalizingStep));
    //     system->setEquilibration            (equilibration);
    //     system->setAnalytical               (analyticOrNot);
    //     system->runMetropolisSteps          (numberOfSteps, firstCriteria, 
    //                                         importanceOrNot, allEnergiesOrNot, stepLength);

    //     firstCriteria = 1;
        
    //     energyNew = system->getSampler()->getEnergy();
    //     energyDerivative = system->getSampler()->getDerivative();
    //     alphaNew = alpha - minimizationRate_alpha*energyDerivative[0];
    //     betaNew = beta - minimizationRate_beta*energyDerivative[1];

        
    //     file.open (file_name, ios::out | ios::app);
    //     file << alpha << "\t" << beta << "\t" << energy << "\t" << energyDerivative[0] << "\t" << energyDerivative[1] << "\n";
    //     file.close();

    //     file2.open (energy_file, ios::out | ios::app);
    //     file2 << alpha << "\t" << beta << "\t" << energy << "\t" << system->getSampler()->getKineticEnergy() 
    //     << "\t" << system->getSampler()->getPotentialEnergy()  << "\t" << system->getSampler()->getInteractionEnergy()<< "\n";
    //     file2.close();

    //     cout << alpha << "\t" << beta << "\t" << energy << "\t" << energyDerivative[0] << "\t" << energyDerivative[1] << "\n";

    //     energyChange = std::abs(energyNew - energy);
    //     alpha = alphaNew;
    //     beta = betaNew;
    //     energy = energyNew;
    // } 





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