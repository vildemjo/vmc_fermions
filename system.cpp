#include "system.h"
#include <cassert>
#include <cmath>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <string>
#include <iostream>
#include <omp.h>

using namespace std;

bool System::metropolisStep() {
    /* This function contains the actual metropolis step. Here a random particle is
     chosen and moved randomly in the dimensions included in the simulation. Afterwards
     the acceptance of the move is checked using the standard Metropolis algorithm check. 
     If the move is accepted, the function returns true and lets function runMetropolisstep
     know that the step was accepted. If it is not accepted the particle is moved back to 
     its original position and the function return false. The implies that the sampling is 
     done with the wavefunction made up of particles at the original position.*/


    // A vector to save the distances the particle is moved in
    // case it has to be moved back
    std::vector<double> randomAmount(m_numberOfDimensions);
    double   oldWaveFunction      = m_waveFunction->evaluate();
    int      randomParticleIndex  = Random::nextInt(m_numberOfParticles);



    for(int m1=0;m1<m_numberOfDimensions; m1++){
        randomAmount[m1] = m_stepLength*(Random::nextDouble()-0.5);
        m_particles[randomParticleIndex]->adjustPosition(randomAmount[m1], m1);
    }


    m_waveFunction->updateSlaterRelatedThings(randomParticleIndex);

    // cout << "ok after update \n";

    double newWaveFunction = m_waveFunction->evaluate();
    
    // std::cout << "New wavefunction: " << newWaveFunction << endl;

    // Determening if step is accepted (return true) or not (move particle back and return false)
    if (Random::nextDouble() <= m_waveFunction->computeRatio(oldWaveFunction, newWaveFunction)){
                                    // cout << "ok" << endl;
        return true;
        }

    for(int m2=0;m2<m_numberOfDimensions; m2++){
        m_particles[randomParticleIndex]->adjustPosition(-randomAmount[m2], m2);
    }

    m_waveFunction->updateSlaterRelatedThings(randomParticleIndex);
    // cout << "no move" << endl;
    return false;
}

bool System::metropolisStepImportance() {
    /* This function contains the actual metropolis step. Here a random particle is
     chosen and moved randomly in the dimensions included in the simulation. Afterwards
     the acceptance of the move is checked using importance sampling/Metropolis-Hastings
     algorithm check which includes Green's function and a quantum force/drift force 
     If the move is accepted, the function returns true and lets function 
     runMetropolisstepImportance know that the step was accepted. If it is not accepted the
     particle is moved back to its original position and the function return false. The 
     implies that the sampling is done with the wavefunction made up of particles at the 
     original position.*/


    // A vector to save the distances the particle is moved in
    // case it has to be moved back
    std::vector<double> importanceAmount(m_numberOfDimensions);

    int particleIndex = Random::nextInt(m_numberOfParticles);

    // std::cout << "Particle index: " << particleIndex << endl;

    double oldWaveFunction    = m_waveFunction->evaluate();
    auto   oldQuantumForce    = m_hamiltonian->computeQuantumForce(particleIndex);
    auto   oldPosition        = m_particles[particleIndex]->getPosition();

    for(int m1=0;m1<m_numberOfDimensions; m1++){
        importanceAmount[m1] = Random::nextGaussian(0,1)*sqrt(m_timeStep)+oldQuantumForce[m1]*m_timeStep*m_diffConstant;
        m_particles[particleIndex]->adjustPosition(importanceAmount[m1], m1);
    }
 
    m_waveFunction->updateSlaterRelatedThings(particleIndex);
    double newWaveFunction  = m_waveFunction->evaluate();
    auto   newQuantumForce  = m_hamiltonian->computeQuantumForce(particleIndex);
    auto   newPosition      = m_particles[particleIndex]->getPosition();
    

    double greensFunctionFrac = greensFunctionFraction(oldPosition, oldQuantumForce,
                                                         newPosition, newQuantumForce);
    // Determening if step is accepted (return true) or not (move particle back and return false)
    if (Random::nextDouble() <= greensFunctionFrac*m_waveFunction->computeRatio(oldWaveFunction, newWaveFunction)){
        return true;
        }

    
    for(int m4=0;m4<m_numberOfDimensions; m4++){
        m_particles[particleIndex]->adjustPosition(-importanceAmount[m4], m4);
    }

    m_waveFunction->updateSlaterRelatedThings(particleIndex);

    return false;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps, int firstCriteria, bool importanceOrNot,
                                                         bool allEnergiesOrNot, double stepLength) {
    /* This function runs through the Monte Carlo cycles and performs the metropolis steps
    through the function metropolisStep or metropolisStepImportance. Here the energy and the 
    information needed to evaluate the one-body density is sampled in the Sampler class and 
    the result is printed to file. */

    m_particles                             = m_initialState->getParticles();
    m_sampler                               = new Sampler(this);
    // The eqilibration is set to be a number of steps instead of a fraction of the steps
    // and here it is added onto the number of steps.
    m_numberOfMetropolisSteps               = numberOfMetropolisSteps + m_equilibration; 
    m_stepLength                            = stepLength;
    m_sampler->setNumberOfMetropolisSteps   (m_numberOfMetropolisSteps);
    m_sampler->setFileOutput                (firstCriteria);
    setImportance                           (importanceOrNot);
    setAllEnergies                          (allEnergiesOrNot);
    m_waveFunction->setupSlaterRelatedThings();

 



    for (int i = 0; i < m_numberOfMetropolisSteps; i++) {
        
        m_steps += 1;

        bool acceptedStep;

        if (getImportance() == true){
            m_timeStep = getStepLength();
            acceptedStep = metropolisStepImportance();
        }else{
            acceptedStep = metropolisStep();
        }

        // There are two different runs, one where the local enery in every step 
        // is saved to file (true) and another where only the the expectation value
        // is saved (false).
        if (getAllEnergies() == true){
            m_sampler->sampleAllEnergies(acceptedStep);
        }else{
            m_sampler->sample(acceptedStep);
        }
    }
    // std::cout << "finished MC loop for alpha "<< getWaveFunction()->getParameters()[0] << std::endl;
    
    m_sampler->computeAverages();

    //  std::cout << "computes averages" << std::endl;

    if (getAllEnergies() == true){
        m_sampler->printOutputToEnergyFile();
        // m_sampler->printOneBodyDensityToFile();
    }else{
        m_sampler->printOutputToTerminal();
    }

}



void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibration(double equilibration) {
    assert(equilibration >= 0);
    m_equilibration = equilibration;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setFrequency(double omega) {
    m_omega = omega;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

void System::setAnalytical(bool statement){
    m_analytical = statement;
}

void System::setImportance(bool statement){
    m_importance = statement;
}

void System::setAllEnergies(bool statement){
    m_allEnergies = statement;
}

double System::greensFunctionFraction(std::vector<double> posNew, std::vector<double> posOld, std::vector<double> forceNew, std::vector<double> forceOld){
    /* This function calculates the fraction between the Green's function for the transition
    from the old state to the new and the transition from the new to the old. This expression
    is used to determine wether a move is accepted or not when importance sampling is used. */
    
    double exponent = 0;

    for (int n10=0; n10<m_numberOfDimensions; n10++){
        exponent += 0.5*(forceOld[n10]+forceNew[n10])*
	            (m_diffConstant*m_timeStep*0.5*(forceOld[n10]-forceNew[n10])-posNew[n10]+posOld[n10]);
    }
    return exp(exponent);
}

void System::setSpinFactor(double spinFactor){
    m_spinFactor = spinFactor;
}