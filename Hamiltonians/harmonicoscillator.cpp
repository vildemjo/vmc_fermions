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
    /* This function calculated the local energy for the wavefunction with particles 
    in the positions given by the input. The potential energy is always calculated the 
    same way, but the kinetic energy can be calcualted both analytically and numerically.
    This is determined by a bool statement parameter in the System class. */

    double m = 1.0;           // The boson mass, but in natural units
    double rSum2 = 0.0;
    double potentialEnergy = 0;

    auto m_particles = m_system->getParticles();

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector <double> particlePosition(numberOfDimensions);

    for (int p1 = 0; p1 < numberOfParticles; p1++){
        particlePosition = m_particles[p1]->getPosition();

        for (int p2 = 0; p2 < numberOfDimensions; p2++){
            rSum2 += particlePosition[p2]*particlePosition[p2];
        }
    }

    potentialEnergy = 0.5*m*m_omega*m_omega*rSum2;


    return potentialEnergy;
}

double HarmonicOscillator::computeKineticEnergy() {
    /* This function calculated the local energy for the wavefunction with particles 
    in the positions given by the input. The potential energy is always calculated the 
    same way, but the kinetic energy can be calcualted both analytically and numerically.
    This is determined by a bool statement parameter in the System class. */

    double hbar = 1.0;        // Planck's constant, but in natural units
    double m = 1.0;           // The boson mass, but in natural units

    double doubleDerivative = 0.0;
    double kineticEnergy = 0;

    auto m_particles = m_system->getParticles();

    if (m_system->getAnalytical() == true){
        doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative();
    }else{
        doubleDerivative = computeDoubleDerivativeNumerically();
    }

    kineticEnergy = -0.5*(hbar*hbar/m)*doubleDerivative;

    return kineticEnergy;
}

double HarmonicOscillator::computeLocalEnergy(){
    double kineticEnergy = computeKineticEnergy();
    double potentialEnergy = computePotentialEnergy();

    return kineticEnergy + potentialEnergy;
}

std::vector<double> HarmonicOscillator::computeQuantumForce(int particleIndex){
    /* This function calculates the quantum force/drift force with is used for importance
        sampling. The quantum force is given by the derivative of the wavefunction. */
    
     auto derivative = m_system->getWaveFunction()->computeDerivative(particleIndex);

    for (int m=0;m<m_system->getNumberOfDimensions();m++){
        derivative[m] *= 2;
    }

    return derivative;
}
