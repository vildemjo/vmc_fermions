#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include <string>
#include <omp.h>

using namespace std;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}


void Sampler::sample(bool acceptedStep) {
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
        m_cumulativeEnergyDerivativeAlpha = 0;
        m_cumulativeEnergyDerivativeBeta = 0;
        m_cumulativeAlphaDerivative = 0;
        m_cumulativeBetaDerivative = 0;
        m_cumulativeDistance = 0;
    }

    double localEnergy = 0;

    // Sampling if the equilibrium stage is passed
    if (m_stepNumber > m_system->getEquilibration()){
        
        // Counting number of accepted steps
        if (acceptedStep == 1) { m_numberOfAcceptedSteps += 1; }

        localEnergy = m_system->getHamiltonian()->computeLocalEnergy();
        // m_cumulativeDistance += m_system->getHamiltonian;

        m_cumulativeEnergy  += localEnergy;
        m_cumulativeEnergySquared += localEnergy*localEnergy;
        // m_cumulativeEnergyDerivativeAlpha += localEnergy*m_system->getWaveFunction()
        //                                                 ->computeAlphaDerivative();
                                                        
        // m_cumulativeEnergyDerivativeBeta += localEnergy*m_system->getWaveFunction()
        //                                                 ->computeBetaDerivative();
        // m_cumulativeAlphaDerivative += m_system->getWaveFunction()->computeAlphaDerivative();
        // m_cumulativeBetaDerivative += m_system->getWaveFunction()->computeBetaDerivative();
    }

    m_stepNumber++;
}

void Sampler::sampleAllEnergies(bool acceptedStep) {
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
        m_cumulativeEnergyDerivativeAlpha = 0;
        m_cumulativeEnergyDerivativeBeta = 0;
        m_cumulativeAlphaDerivative = 0;
        m_cumulativeBetaDerivative = 0;
        m_localEnergyVector = std::vector <double>();
        m_localEnergyVector.reserve(m_numberOfMetropolisSteps);
    }
    
    // Sampling if the equilibrium stage is passed
    if (m_stepNumber >= m_system->getEquilibration()){

        // Counting number of accepted steps
        if (acceptedStep == 1) { m_numberOfAcceptedSteps += 1; }

        // When all energies are saved, the one-body density data is also acquired 
        m_system->getWaveFunction()->updateOneBodyDensity();

        double localEnergy = m_system->getHamiltonian()->computeLocalEnergy();
        
        m_cumulativeEnergy  += localEnergy;
        m_localEnergyVector.push_back(localEnergy);
        m_cumulativeEnergySquared += localEnergy*localEnergy;
        // m_cumulativeEnergyDerivativeAlpha += localEnergy*m_system->getWaveFunction()
        //                                                 ->computeAlphaDerivative();
                                                        
        // m_cumulativeEnergyDerivativeBeta += localEnergy*m_system->getWaveFunction()
        //                                                 ->computeBetaDerivative();
        // m_cumulativeAlphaDerivative += m_system->getWaveFunction()->computeAlphaDerivative();
        // m_cumulativeBetaDerivative += m_system->getWaveFunction()->computeBetaDerivative();
        
    }

    m_stepNumber++;
    
}

void Sampler::printOutputToTerminal() {
    /* This function prints the relevant information about the set-up and the results
    to the terminal. */

    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibration();
    double  ac = 100*m_numberOfAcceptedSteps/ms;
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(ef) << endl;
    cout << " Accepted steps: " << ac << " %" << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int o=0; o < p; o++) {
        cout << " Parameter " << o+1 << " : " << pa.at(o) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Variance : " <<  m_energySquared - m_energy*m_energy << endl;
    cout << endl;
}

void Sampler::printOutputToEnergyAlphaFile(){
    /* This function makes a file of the expectation energy 
    and the alphas used in the calculation and it also includes 
    the percent of steps that are accepted. */
    
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    double  ef = m_system->getEquilibration();

    ofstream myfile;
    double alpha = m_system->getWaveFunction()->getParameters()[0];

    string filename = "Output/" + m_system->getFileName() + "_energy_alpha.txt";

    if (m_firstCriteria == 0) { 
        
        myfile.open (filename, ios::out | ios::trunc);
        myfile << "  -- System info -- " << endl;
        myfile << " Number of particles  : " << np << endl;
        myfile << " Number of dimensions : " << nd << endl;
        myfile << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
        myfile << " Number of equilibration steps  : 10^" << std::log10(ef) << endl;
        myfile << " Step length/ time step: " << m_system->getStepLength() << endl;
        myfile << "===================================================" << endl;
        myfile << "Energy: \t Alpha: \t Acceptance [%]: \n"; 
        myfile << m_energy << "\t" << alpha << "\t" << 100.0*(double)m_numberOfAcceptedSteps/(double)m_numberOfMetropolisSteps <<"\n";
        myfile.close(); 
    }else{
        // std::cout << "printing" << std::endl;
    myfile.open (filename, ios::out | ios::app);
        
    myfile << m_energy << "\t" << alpha << "\t" << 100.0*(double)m_numberOfAcceptedSteps/(double)m_numberOfMetropolisSteps <<"\n";
    myfile.close();
    std::cout << m_energy << "\t" << alpha <<"\n";
    }
}

void Sampler::printOutputToEnergyFile(){
    /* This function prints the local energy in every step to file 
    if this is the type of run used. */

    ofstream myfile;

    
    string filename = m_system->getFileName() + "_energy.txt";

    myfile.open (filename, ios::out | ios::trunc);

    // adding some extra information to the energy file
    myfile << " # mean: " << m_energy << "\n";
    myfile << " # std: " << "sqrt(" << m_energySquared << " - " << m_energy << "^2) = " << sqrt(m_energySquared - m_energy*m_energy) << "\n";


    for (int n3 = 0; n3<(int) m_localEnergyVector.size(); n3++){
        myfile << m_localEnergyVector[n3] << "\n";
    }
    myfile.close();
}

void Sampler::printOneBodyDensityToFile(){
    /* This function prints the one-body density data to file for 
    the runs where all local energies are saved to file */

    std::vector <std::vector <double>> oneBodyDensity = m_system->getWaveFunction()->getOneBodyDensity();
    
    // The number of cycles included are the number of cycles where the energy 
    // and the one-body density data is acquired. The function is called here to
    // update m_numberOfCyclesIncluded which is used to normalize the data.
    evaluateNumberOfCyclesIncluded();

    ofstream myfile;

    string filename = m_system->getFileName() + "_density.txt";

    myfile.open (filename, ios::out | ios::trunc);

    int numberOfParticles = (double) m_system->getNumberOfParticles();

    for (int n5 = 0; n5 < (int)oneBodyDensity[0].size(); n5++){
        myfile << oneBodyDensity[0][n5] << "\t";
        for (int m5 = 1; m5 < m_system->getNumberOfDimensions()+1; m5++){
            myfile << oneBodyDensity[m5][n5]/((double) m_numberOfCyclesIncluded*(double) numberOfParticles) << "\t";
        }
        myfile << "\n";
    }

    myfile.close();
}

void Sampler::computeAverages() {
    /* This function computes the averages of the sampled quantities. */

    // The number of cycles included are the number of cycles where the energy 
    // and the one-body density data is acquired. The function is called here to
    // update m_numberOfCyclesIncluded which is used to normalize the data.
    evaluateNumberOfCyclesIncluded();

    m_energy = (m_cumulativeEnergy / (double) m_numberOfCyclesIncluded);
    m_energySquared = m_cumulativeEnergySquared / (double) m_numberOfCyclesIncluded;

    // The derivative is used by the gradient descent methods
    // m_derivativeAlpha = 2* m_cumulativeEnergyDerivativeAlpha / (double) m_numberOfCyclesIncluded
                    // - 2*(m_cumulativeAlphaDerivative / (double) m_numberOfCyclesIncluded)*m_energy;
    // m_derivativeBeta =  2* m_cumulativeEnergyDerivativeBeta / (double) m_numberOfCyclesIncluded
                    // - 2*(m_cumulativeBetaDerivative / (double) m_numberOfCyclesIncluded)*m_energy;
    // m_derivative[0] = m_derivativeAlpha;
    // m_derivative[1] = m_derivativeBeta;
}

void Sampler::setFileOutput(int firstCriteria){
    /* This function is used to print the header in the printOutputToEnergyAlphaFile()
    only once. */
    m_firstCriteria = firstCriteria;
}

void Sampler::evaluateNumberOfCyclesIncluded(){
    /* This function finds the number of cycles where the energy is sampled. */
    m_numberOfCyclesIncluded = m_numberOfMetropolisSteps - m_system->getEquilibration();
}
