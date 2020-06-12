#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_system->setFrequency(omega);
    m_omega = m_system->getFrequency();
}

double HarmonicOscillator::computePotentialEnergy() {
    /* This function calculates the potential energy from the 
    external harmonic oscillator potential. */
    
    double rSum2 = 0.0;
    double potentialEnergy = 0;

    auto m_particles = m_system->getParticles();

    int numberOfParticles = (int) m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector <double> particlePosition(numberOfDimensions);

    for (int p1 = 0; p1 < numberOfParticles; p1++){
        particlePosition = m_particles[p1]->getPosition();

        for (int p2 = 0; p2 < numberOfDimensions; p2++){
            rSum2 += particlePosition[p2]*particlePosition[p2];
        }
    }

    potentialEnergy = 0.5*m_omega*m_omega*rSum2;  // Mass is in atomic units i.e. = 1. 


    return potentialEnergy;
}

double HarmonicOscillator::computeKineticEnergy() {
    /* This function calculates the kinetic energy. It can be calculated both analytically and numerically.
    This is determined by a bool statement parameter in the System class. */

    double doubleDerivative = 0.0;
    double kineticEnergy = 0;

    if (m_system->getAnalytical() == true){
        doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative();
    }else{
        doubleDerivative = computeDoubleDerivativeNumerically();
    }

    kineticEnergy = -0.5*doubleDerivative; // Plack's constant and mass are in atomic units i.e. = 1. 

    return kineticEnergy;
}

double HarmonicOscillator::computeLocalEnergy(){
    /* This function calculates the local energy. */

    double kineticEnergy = computeKineticEnergy();
    double potentialEnergy = computePotentialEnergy();

    return kineticEnergy + potentialEnergy;
}


