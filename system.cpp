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

System::System(){
    m_random = new Random();
}

System::System(int seed){
    m_random = new Random(seed);
}

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
    int      randomParticleIndex  = m_random->nextInt(0, m_numberOfParticles-1);

    // std::cout << "\n ------------- \n Particle index: " << randomParticleIndex << "\n ------------ " << endl;
    // std::cout << "Old wavefunction: " << oldWaveFunction << endl;


    for(int m1=0;m1<m_numberOfDimensions; m1++){
        randomAmount[m1] = m_stepLength*(m_random->nextDouble()-0.5);
        m_particles[randomParticleIndex]->adjustPosition(randomAmount[m1], m1);
    }

    // m_waveFunction->updateSlaterRelatedThings(randomParticleIndex);
    m_waveFunction->updateSlaterMatrix(randomParticleIndex);
    // cout << "ok after update \n";

    double newWaveFunction = m_waveFunction->evaluate();
    
    // std::cout << "New wavefunction: " << newWaveFunction << endl;
    

    // cout << "Old ratio: " << newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction) << "\n";
   
    // Determening if step is accepted (return true) or not (move particle back and return false)
    if (m_random->nextDouble() <= m_waveFunction->computeRatio(newWaveFunction, oldWaveFunction)){
                                    // cout << "ok" << endl;
        m_waveFunction->updateInverseSlaterMatrix(randomParticleIndex);
        return true;
        }

    for(int m2=0;m2<m_numberOfDimensions; m2++){
        m_particles[randomParticleIndex]->adjustPosition(-randomAmount[m2], m2);
    }

    m_waveFunction->updateSlaterMatrix(randomParticleIndex);
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

    int particleIndex = m_random->nextInt(m_numberOfParticles-1);



    // cout << "----- \n new test \n -----\n";
    // std::cout << "Particle index: " << particleIndex << endl;

    double oldWaveFunction    = m_waveFunction->evaluate();
    auto   oldQuantumForce    = m_waveFunction->computeQuantumForce(particleIndex, true); // true = old quantum force
    auto   oldPosition        = m_particles[particleIndex]->getPosition();

    for(int m1=0;m1<m_numberOfDimensions; m1++){
        importanceAmount[m1] = m_random->nextGaussian(0,1)*sqrt(m_timeStep)+oldQuantumForce[m1]*m_timeStep*m_diffConstant;
        m_particles[particleIndex]->adjustPosition(importanceAmount[m1], m1);
    }
 
    m_waveFunction->updateSlaterMatrix(particleIndex);
    m_waveFunction->updateInverseSlaterMatrix(particleIndex);
    double newWaveFunction  = m_waveFunction->evaluate();
    auto   newQuantumForce  = m_waveFunction->computeQuantumForce(particleIndex, false); // false = new quantum force
    auto   newPosition      = m_particles[particleIndex]->getPosition();
    

    double greensFunctionFrac = greensFunctionFraction(newPosition, oldPosition, newQuantumForce, oldQuantumForce);

    // cout << greensFunctionFrac << "\t" << m_waveFunction->computeRatio(oldWaveFunction, newWaveFunction) << "\t" << oldWaveFunction << "\t" << newWaveFunction << "\n";
    // Determening if step is accepted (return true) or not (move particle back and return false)
    if (m_random->nextDouble() <= greensFunctionFrac*m_waveFunction->computeRatio(newWaveFunction, oldWaveFunction)){
                // cout << "move accepted \n";
        // m_waveFunction->updateInverseSlaterMatrix(particleIndex);
        return true;
        }
    // cout << "move not accepted \n";
                // std::cout << "Particle index: " << particleIndex << endl;
    
    for(int m4=0;m4<m_numberOfDimensions; m4++){
        m_particles[particleIndex]->adjustPosition(-importanceAmount[m4], m4);
    }

    m_waveFunction->updateSlaterMatrix(particleIndex);
    m_waveFunction->updateInverseSlaterMatrix(particleIndex);

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
    setSpinFactor();
    //  std::cout << "ok before setup";
    m_waveFunction->setupSlaterRelatedThings();
    //  std::cout << "ok after setup"<< " \n";
 



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
            // std::cout << "sample ok: " << m_hamiltonian->computeLocalEnergy() << " \n";
        }
    }
    // std::cout << "finished MC loop for alpha "<< getWaveFunction()->getParameters()[0] << std::endl;
    
    m_sampler->computeAverages();

    //  std::cout << "computes averages" << std::endl;

    if (getAllEnergies() == true){
        m_sampler->printOutputToEnergyFile();
        m_sampler->printOneBodyDensityToFile();
        m_sampler->printOutputToTerminal();
    }else{
        // m_sampler->printOutputToTerminal();
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
            // std::cout << "eq:" << m_equilibration <<" \n";
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
    
    // double exponentTop = 0;
    // double exponentLow = 0;
    double exponent = 0;
    double sigma = m_diffConstant*m_timeStep;

    for (int n10=0; n10<m_numberOfDimensions; n10++){
        // auto posChange = posOld[n10]-posNew[n10];
        // exponentTop += -(posOld[n10]-posNew[n10]-sigma*forceNew[n10])*(posOld[n10]-posNew[n10]-sigma*forceNew[n10])/(4.0*sigma);
        // exponentLow += -(posNew[n10]-posOld[n10]-sigma*forceOld[n10])*(posNew[n10]-posOld[n10]-sigma*forceOld[n10])/(4.0*sigma);

        exponent += 0.5*(forceOld[n10]+forceNew[n10])*
	            (m_diffConstant*m_timeStep*0.5*(forceOld[n10]-forceNew[n10])-posNew[n10]+posOld[n10]); // fra lecture notes
        // exponent += 0.5*(forceOld[n10]-forceNew[n10])*(posNew[n10]-posOld[n10]); // Fra Evens master (+ legg til number of dim)
        // exponent += 0.25*m_diffConstant*m_timeStep*(forceNew[n10]*forceNew[n10]+forceOld[n10]*forceOld[n10]) 
        //             + 0.5*(posOld[n10]-posNew[n10])*(forceOld[n10]+forceNew[n10]); // Fra meg selv
    }
    return exp(exponent);///exp(exponentLow);// + m_numberOfDimensions;
}

void System::setSpinFactor(){
    
    
    std::vector <std::vector <double> > spinFactor(m_numberOfParticles, std::vector<double>(m_numberOfParticles, (double) 0));
    // m_spinFactor = spinFactor;

    for (int k5 = 0; k5 < m_numberOfParticles; k5++){
        for (int k6 = 0; k6 < m_numberOfParticles; k6++){
            if (k6 < m_numberOfParticles/2 && k5 >= m_numberOfParticles/2){
                spinFactor[k5][k6] = 1;
            }if (k6 < m_numberOfParticles/2 && k5 < m_numberOfParticles/2){
                spinFactor[k5][k6] = (double)1/3;
            }if (k6 >= m_numberOfParticles/2 && k5 >= m_numberOfParticles/2){
                spinFactor[k5][k6] = (double)1/3;
            }if (k6 >= m_numberOfParticles/2 && k5 < m_numberOfParticles/2){
                spinFactor[k5][k6] = 1;
            }
            // cout << "a_" <<k5 << k6 << ": " << spinFactor[k5][k6] << "\n";
        }
    }
    
    m_spinFactor = spinFactor;
   
}