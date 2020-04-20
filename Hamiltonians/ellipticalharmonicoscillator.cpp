#include "ellipticalharmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

EllipticalHarmonicOscillator::EllipticalHarmonicOscillator(System* system, double omega, double gamma) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
    m_gamma   = gamma;
}

double EllipticalHarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* This function calculated the local energy for the wavefunction with particles 
    in the positions given by the input. The potential energy is always calculated the 
    same way, but the kinetic energy can be calcualted both analytically and numerically.
    This is determined by a bool statement parameter in the System class. */

    double hbar = 1.0;        // Planck's constant, but in natural units
    double m = 1.0;           // The boson mass, but in natural units

    double rSum2 = 0.0;
    double doubleDerivative = 0.0;

    double kineticEnergy = 0;
    double potentialEnergy = 0;

    particles = m_system->getParticles();

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector <double> particlePosition(numberOfDimensions);

    for (int p1 = 0; p1 < numberOfParticles; p1++){
        particlePosition = particles[p1]->getPosition();

        for (int p2 = 0; p2 < numberOfDimensions; p2++){
            if (p2 == 2){ // This is for the elliptical trap case, else beta = 1
                rSum2 += m_gamma*particlePosition[p2]*particlePosition[p2];
            }else{
                rSum2 += particlePosition[p2]*particlePosition[p2];
            }
        }
    }

    potentialEnergy = 0.5*m*m_omega*m_omega*rSum2;

    // if (m_system->getSteps() == 4533+(int)1e5 || m_system->getSteps() == 4532+(int)1e5){
    //     std::cout << "potentialEnergy: " << potentialEnergy << std::endl;
    // }

    if (m_system->getAnalytical() == true){
        doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative();
    }else{
        doubleDerivative = computeDoubleDerivativeNumerically(particles);
    }

    // if (m_system->getSteps() >= 4532+(int)1e5-10 && m_system->getSteps() <= 4532+(int)1e5+10){
    //     std::cout << "doubleDerivative: " << doubleDerivative << "  (step: " << m_system->getSteps()-1e5 << ")" << std::endl;
    // }

    kineticEnergy = -0.5*(hbar*hbar/m)*doubleDerivative;

    return potentialEnergy + kineticEnergy; // Do not need to add the interaction energy
                                            // because the energy only sampled when the
                                            // interaction energy is zero. This is fixed by
                                            // the wavefunction becoming zero at the same 
                                            // time.
}

std::vector<double> EllipticalHarmonicOscillator::computeQuantumForce(int particleIndex,std::vector<class Particle*> particles){
    /* This function calculates the quantum force/drift force with is used for importance
        sampling. The quantum force is given by the derivative of the wavefunction. */
    
     auto derivative = m_system->getWaveFunction()->computeDerivative(particleIndex);

    for (int m=0;m<m_system->getNumberOfDimensions();m++){
        derivative[m] *= 2;
    }

    return derivative;
}
