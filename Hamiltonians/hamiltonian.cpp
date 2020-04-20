#include "hamiltonian.h"
#include "../particle.h"
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include <iostream>

using std::cout;
using std::endl;

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeDoubleDerivativeNumerically(std::vector<Particle*> particles) {
    /* This function calculates the double derivative of the wavefunction numerically.
    This is used to calculate the local energy for the system. The function hence returns 
    the double derivative divided by the wavefunction. */

    double step = 1e-3;
    double dpsidr2 = 0.0;


    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double waveNext = 0;
    double waveLast = 0;
    double waveCurrent = 0;

    for(int i4=0; i4<numberOfParticles; i4++){
         for(int n4=0; n4<numberOfDimensions;n4++){
            // Wave function at forward step
            particles[i4]->adjustPosition(step, n4);
            waveNext = m_system->getWaveFunction()->evaluate();

            // Wave function at backward step
            particles[i4]->adjustPosition(-2*step, n4);
            waveLast = m_system->getWaveFunction()->evaluate();

            // Calculating the part of the double derivative which involves psi(x+dx) and psi(x-dx)
            dpsidr2 += (waveNext+waveLast)/(step*step);

            // Resetting the particle positions so that a new particle and spesific dimension can be calculated
            particles[i4]->adjustPosition(step, n4); 
        }
    }

    waveCurrent = m_system->getWaveFunction()->evaluate();

    // Calculating the part of the double derivative which involves psi(x)
    dpsidr2 += -2*numberOfParticles*numberOfDimensions*waveCurrent/(step*step);

    return (1/waveCurrent)*dpsidr2;
 
}

