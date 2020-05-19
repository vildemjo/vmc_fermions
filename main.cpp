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
#include <cmath>

using namespace std;

void alphaListRun(string filename, int numberOfP);
string setMethodName(bool analyticOrNot);
string setSamplingName(bool importanceOrNot);

int main() {
/* The standard set-up */

    bool analyticOrNot       = true;
    bool importanceOrNot     = false;
    bool allEnergiesOrNot    = true;
    int equilibration        = 1e5;          // Number of the total steps used for equilibration
    int numberOfDimensions   = 2;
    int numberOfParticles    = 2;
    int numberOfSteps        = (int) pow(2.0,22.0);
    double omega             = 1.0;          // Oscillator frequency.
    double stepLength        = 1.0;//0.001;          // Metropolis step length.
    int firstCriteria        = 0;            // print header in file
    double alpha             = 0.9;
    double beta              = 0;   
    double inititalizingStep = 0.5;

    double spinFactor  = 1;
    int numberOfBins = 500;
    double densityLength = 4.0;
 

//    /* Set-up to run and save local energies for every step to file*/

//     allEnergiesOrNot    = true;
//     importanceOrNot     = true;

//     clock_t start, end;
//     // Recording the starting clock tick.
//     start = clock();

    
//     System* system = new System();
//     system->setHamiltonian                (new HarmonicOscillator(system, omega));
//     system->setWaveFunction               (new SlaterDeterminant(system, alpha));

//     system->setInitialState               (new GaussianDistribution(system, numberOfDimensions, 
//                                                 numberOfParticles, inititalizingStep));

//     system->setEquilibration              (equilibration);
//     system->setAnalytical                 (analyticOrNot);
//     system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
//     system->setFileName                   ("Output/exercise_f/slater_2p");

//     system->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
//                                                             allEnergiesOrNot, stepLength);

//     cout << "number of steps: " << numberOfSteps << endl;
//     cout << "energy: " << system->getSampler()->getEnergy() << endl;

//     end = clock();
//     double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
//     cout << "CPU time: " << time_taken << " seconds" << endl;



    /* Same but for the interacting case */   


    // importanceOrNot = false;
    // stepLength = 0.5;
    // inititalizingStep = 0.5;
    // allEnergiesOrNot = true;
    // omega = 0.1;

    // spinFactor  = 1.0;
    // string file_name = "Output/exercise_e/interaction_ground_state_2p_omega_10";

    // // interaction or spherical trap (2.82843 or 1.0)
    // beta = 0.17917; // 0.34998; //0.39975;    // omega_normal^2/omega_ho^2
    // alpha = 0.07937; // 0.48521; //0.98834;

    // clock_t start2, end2;

    // // Recording the starting clock tick.
    // start2 = clock();


    // System* system2 = new System();
    // system2->setHamiltonian                (new InteractionHarmonicOscillator(system2, omega));
    // system2->setWaveFunction               (new SlaterDeterminantInteraction(system2, alpha, beta));

    // system2->setInitialState               (new RandomUniform(system2, numberOfDimensions, 
    //                                             numberOfParticles, inititalizingStep));
    // system2->setEquilibration              (equilibration);
    // system2->setAnalytical                 (analyticOrNot);
    // system2->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
    // system2->setFileName                   (file_name);

    // system2->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
    //                                                         allEnergiesOrNot, stepLength);

    // // cout << "energy: " << system2->getSampler()->getEnergy() << endl;

    // end2 = clock();
    // double time_taken2 = double(end2 - start2) / double(CLOCKS_PER_SEC); 
    // cout << "CPU time: " << time_taken2 << " seconds" << endl;

// ---------------------------------------------------------------------------------------------------------------------
/* Printing all local energy measurements for some alphas */
// With and without analytical
// ---------------------------------------------------------------------------------------------------------------------

    // double alphaStart = 0.5;
    // double alphaStop = 1.5;
    // double alphaStep = 0.1;
    
    // allEnergiesOrNot = true;
    // numberOfSteps = pow(2,21);
    // analyticOrNot = true;

    // ofstream myfile;

    // string methodName = setMethodName(analyticOrNot);
    // string cpufilename = "Output/exercise_c/"+ setMethodName(analyticOrNot)+ ".txt";

    // std::vector <int> Ns = {2};//, };
    // std::vector <int> ds = {2};


    // myfile.open( cpufilename, ios::out | ios::trunc);
    // myfile.close();     

    // alpha = alphaStart;
    // for(int a=0; alpha <= alphaStop; a++){

    //     double alphaPrint = alpha*100.0;
    //     int alphaPrintable = ceil(alphaPrint);    
    //     string filename = "Output/exercise_c/allEnergies/" + 
    //         methodName + "2d_" +"2p"+ "_alpha_" + to_string(alphaPrintable) + "_MC_21";

    //     numberOfDimensions  = 2;
    //     numberOfParticles   = 2;
        
    //     allEnergiesOrNot    = true;
    //     importanceOrNot     = false;

    //     clock_t start, end;
    //     // Recording the starting clock tick.
    //     start = clock();

        
    //     System* system = new System();
    //     system->setHamiltonian                (new HarmonicOscillator(system, omega));
    //     system->setWaveFunction               (new SlaterDeterminant(system, alpha));

    //     system->setInitialState               (new RandomUniform(system, numberOfDimensions, 
    //                                                 numberOfParticles, inititalizingStep));

    //     system->setEquilibration              (equilibration);
    //     system->setAnalytical                 (analyticOrNot);
    //     system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
    //     system->setFileName                   (filename);

    //     system->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
    //                                                             allEnergiesOrNot, stepLength);

    //     end = clock();
    //     double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 

    //     myfile.open (cpufilename, ios::out | ios::app);
    //     myfile << time_taken << "\n";
    //     myfile.close();
        
    //     alpha += alphaStep;
    //     firstCriteria = 1;

    // }
// ---------------------------------------------------------------------------------------------------------------------
/* Printing all local energy measurements for some alphas */
// With and without importance
// ---------------------------------------------------------------------------------------------------------------------
    // analyticOrNot    = true;
    // allEnergiesOrNot = true;
    // importanceOrNot  = true;

    // numberOfDimensions  = 2;
    // numberOfParticles   = 6;
    // alpha               = 0.9;
    // numberOfSteps       = (int) pow(2.0, 21.0);

    // string methodName = setMethodName(analyticOrNot);
    // string samplingType;

    // std::vector <double> dls = { 1.0, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001};

    // ofstream thisfile;

    // samplingType = setSamplingName(importanceOrNot);

    // string file_name = "Output/exercise_d/" + methodName + 
    //                     to_string(numberOfParticles) + "p_" +
    //                     to_string(numberOfDimensions) + 
    //                     "d_" + samplingType + ".txt";

    // thisfile.open(file_name, ios::out | ios::trunc);
    // thisfile << " Step lenght: \t acceptance [%]: \t energy: \t CPU time: \n";
    // thisfile.close();

    // for (int dl = 0; dl < dls.size(); dl++ ){


    //     double dlPrint = dls[dl]*1000.0;
    //     int dlPrintable = ceil(dlPrint);
    //     string filename = "Output/exercise_d/allEnergies/" + 
    //     methodName + "2d_" +"6p"+ "_stepsize_" + to_string(dlPrintable) + "_MC_21_" + samplingType;

    //     clock_t start, end;
    //     // Recording the starting clock tick.
    //     start = clock();

    //     double timeStep     = dls[dl];          // Metropolis step length.
    //     inititalizingStep = 0.5;

    //     System* system = new System();
    //     system->setHamiltonian                (new HarmonicOscillator(system, omega));
    //     system->setWaveFunction               (new SlaterDeterminant(system, alpha));
    //     if (importanceOrNot == true){
    //         system->setInitialState               (new RandomUniform(system, numberOfDimensions, 
    //                                                 numberOfParticles, inititalizingStep));
    //     }else{
    //         system->setInitialState               (new GaussianDistribution(system, numberOfDimensions, 
    //                                                     numberOfParticles, inititalizingStep));
    //     }
    //     system->setEquilibration              (equilibration);
    //     system->setAnalytical                 (analyticOrNot);
    //     system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
    //     system->setFileName                   (filename);

    //     system->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
    //                                                             allEnergiesOrNot, timeStep);

    //     end = clock();
    //     double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 

    //     thisfile.open(file_name, ios::out | ios::app);

    //     thisfile << dls[dl] << "\t" << 
    //                 system->getSampler()->getAcceptance() << "\t" <<
    //                 system->getSampler()->getEnergy() << "\t"<<
    //                 time_taken << "\n"; 
    //     thisfile.close();

    // }





/* Gradient descent */

/* One parameter */

    // ofstream file;
    
    // std::vector <double> A = {0.9, 0.95, 1.0, 1.01, 1.1, 1.2};// minrate: { 1.0, 0.5, 0.1, 0.05, 0.01, 0.005};

    // for (int aa = 0; aa < A.size(); aa++){
    //     double energyChange = 1.0;
    //     double stopCriteria = 1e-7;
    //     double energyNew = 0.0;
    //     double energyDerivative = 1.0;
    //     double alphaNew = 0;
    //     double minimizationRate = 1.0;
    //     allEnergiesOrNot = false;
    //     importanceOrNot = false;
    //     alpha = 0.9;
    //     stepLength = 0.5;
    //     inititalizingStep = stepLength;

    //     double energy       = 0;

    //     numberOfDimensions  = 2;
    //     numberOfParticles   = 2;
    //     numberOfSteps       = (int) std::pow(2,19.0);
            
    //     alpha = A[aa];

    //     double alphaPrint = alpha*100.0;
    //     int alphaPrintable = ceil(alphaPrint);
    //     double gammaPrint = minimizationRate*1000.0;
    //     int gammaPrintable = ceil(gammaPrint);


    //     string file_name = "Output/exercise_e/gradient_descent_p2_alphastart_"+ to_string(alphaPrintable) + "_gamma_" + to_string(gammaPrintable) + ".txt";

    //     file.open (file_name, ios::out | ios::trunc);
    //     file << "Alpha: \t Energy: \t Derivative: \n";
    //     file.close();
    //     cout << "Alpha: \t Energy: \t Derivative: \n";

    //     for (int k=0;  energyChange > stopCriteria && k <= 80; k++){

    //         clock_t start, end;
    //     // Recording the starting clock tick.
    //         start = clock();
        

    //         System* system = new System();
    //         system->setHamiltonian              (new HarmonicOscillator(system, omega));
    //         system->setWaveFunction             (new SlaterDeterminant(system, alpha));
    //         system->setInitialState             (new RandomUniform(system, numberOfDimensions, 
    //                                                     numberOfParticles, inititalizingStep));
    //         system->setEquilibration            (equilibration);
    //         system->setAnalytical               (analyticOrNot);
    //         system->runMetropolisSteps          (numberOfSteps, firstCriteria, 
    //                                             importanceOrNot, allEnergiesOrNot, stepLength);

    //         firstCriteria = 1;
            
    //         energyNew = system->getSampler()->getEnergy();
    //         energyDerivative = system->getSampler()->getDerivative()[0];
    //         alphaNew = alpha - minimizationRate*energyDerivative/numberOfParticles;

            
    //         file.open (file_name, ios::out | ios::app);
    //         file << alpha << "\t" << energy << "\t" << energyDerivative << "\n";
    //         file.close();

    //         cout << alpha << "\t" << energy << "\t" << energyDerivative << "\n";

    //         energyChange = std::abs(energyNew - energy);
    //         alpha = alphaNew;
    //         energy = energyNew;

    //         end = clock();
    //         double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    //         // cout << "CPU time: " << time_taken << " seconds" << endl;
        
    //     }
        
    // } 









/* Two parameters */

    ofstream  file; //file2, 

    double energyChange = 1.0;
    double stopCriteria = 1e-8;
    double energyNew = 0.0;
    std::vector<double> energyDerivative(2,0);
    double alphaNew = 0;
    double betaNew = 0;
    double minimizationRate_alpha = 0.02; //0.01, 0.05, 0.5, 0.1
    double minimizationRate_beta = 0.01;
    allEnergiesOrNot = false;
    importanceOrNot = false;
    alpha = 0.5;
    beta = 0.2;
    stepLength = 0.5;
    inititalizingStep = 0.5;
    omega = 0.05;

    double energy       = 0;

    numberOfDimensions  = 2;
    numberOfParticles   = 2;
    numberOfSteps       = (int) std::pow(2,19.0);

    double alphaPrint = alpha*100.0;
    int alphaPrintable = ceil(alphaPrint);
    double betaPrint = beta*100.0;
    int betaPrintable = ceil(betaPrint);
    double gammaPrint = minimizationRate_alpha*1000.0;
    int gammaPrintable = ceil(gammaPrint);
    double omegaPrint = omega*100.0;
    int omegaPrintable = ceil(omegaPrint);


    string file_name = "Output/exercise_e/gradient_descent_interaction_p2_omega_"+ to_string(omegaPrintable) +"_alphastart_"+ to_string(alphaPrintable) + "_betastart_"+ to_string(betaPrintable) + "_gamma_" + to_string(gammaPrintable) + ".txt";
    string energy_file = file_name + "_energy";

    file.open (file_name, ios::out | ios::trunc);
    file << "Alpha: \t Beta: \t Energy: \t Derivative (alpha): \t Derivative (beta): \n";
    file.close();
    cout << "Alpha: \t Beta: \t Energy: \t Derivative (alpha): \t Derivative (beta): \n";

    // file2.open (energy_file, ios::out | ios::trunc);
    // file2 << "Alpha: \t Beta: \t Energy: \t Kinetic Energy: \t Potential Energy:  \t Interaction Energy:\n";
    // file2.close();

    for (int k=0;  energyChange > stopCriteria && k <= 200; k++){
    

        System* system = new System();
        system->setHamiltonian              (new InteractionHarmonicOscillator(system, omega));
        system->setWaveFunction             (new SlaterDeterminantInteraction(system, alpha, beta));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, 
                                                    numberOfParticles, inititalizingStep));
        system->setEquilibration            (equilibration);
        system->setAnalytical               (analyticOrNot);
        system->runMetropolisSteps          (numberOfSteps, firstCriteria, 
                                            importanceOrNot, allEnergiesOrNot, stepLength);

        firstCriteria = 1;
        
        energyNew = system->getSampler()->getEnergy();
        energyDerivative = system->getSampler()->getDerivative();
        alphaNew = alpha - minimizationRate_alpha*energyDerivative[0];
        betaNew = beta - minimizationRate_beta*energyDerivative[1];

        
        file.open (file_name, ios::out | ios::app);
        file << alpha << "\t" << beta << "\t" << energy << "\t" << energyDerivative[0] << "\t" << energyDerivative[1] << "\n";
        file.close();

        // file2.open (energy_file, ios::out | ios::app);
        // file2 << alpha << "\t" << beta << "\t" << energy << "\t" << system->getSampler()->getKineticEnergy() 
        // << "\t" << system->getSampler()->getPotentialEnergy()  << "\t" << system->getSampler()->getInteractionEnergy()<< "\n";
        // file2.close();

        cout << alpha << "\t" << beta << "\t" << energy << "\t" << energyDerivative[0] << "\t" << energyDerivative[1] << "\n";

        energyChange = std::abs(energyNew - energy);
        alpha = alphaNew;
        beta = betaNew;
        energy = energyNew;
    } 





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

string setSamplingName(bool importanceOrNot){
    string samplingName;

    if (importanceOrNot == true){
        samplingName = "importance";
    }else{
        samplingName = "brute_force";
    }
return samplingName;
}