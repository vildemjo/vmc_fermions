#include "interactionharmonicoscillator.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

InteractionHarmonicOscillator::InteractionHarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double InteractionHarmonicOscillator::computePotentialEnergy() {
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
// cout << "calcu pot" << endl;

    potentialEnergy = 0.5*m*m_omega*m_omega*rSum2;

    return potentialEnergy;
}


double InteractionHarmonicOscillator::computeKineticEnergy() {
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

// cout << "calcu kin" << endl;

    return kineticEnergy;
}


double InteractionHarmonicOscillator::computeInteractionEnergy() {
    /* This function calculated the local energy for the wavefunction with particles 
    in the positions given by the input. The potential energy is always calculated the 
    same way, but the kinetic energy can be calcualted both analytically and numerically.
    This is determined by a bool statement parameter in the System class. */
    int     numberOfParticles       = m_system->getNumberOfParticles();
    auto distances = calculateInterparticleDistances();
    double  distanceInverted        = 0;
    

    for (int j5 = 0; j5 < numberOfParticles-1; j5++){
        for (int j6 = j5+1; j6 <numberOfParticles; j6++){
            distanceInverted += 1/distances[j5][j6];

        }
    }

    // cout << "calcu int" << endl;

    return distanceInverted;
}


double InteractionHarmonicOscillator::computeLocalEnergy(){
    double kineticEnergy = computeKineticEnergy();
    double potentialEnergy = computePotentialEnergy();
    double interactionEnergy = computeInteractionEnergy();

    return kineticEnergy + potentialEnergy + interactionEnergy;
}


std::vector <std::vector <double> >  InteractionHarmonicOscillator::calculateInterparticleDistances(){
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector <std::vector <double> > distances(numberOfParticles, std::vector<double>(numberOfParticles, (double) 0)); // a matrix of the distance between all particles 
    
    std::vector <class Particle*> m_particles = m_system->getParticles();

    std::vector <double> r1(numberOfDimensions), r2(numberOfDimensions);

    for (int k1 = 0; k1 < numberOfParticles; k1++){

        r1 = m_particles[k1]->getPosition();
        
        for (int k2 = 0; k2 <numberOfParticles; k2++){
            
            if (k1 != k2){
                r2 = m_particles[k2]->getPosition();
                for (int k3 = 0; k3<numberOfDimensions; k3++){
                    distances[k1][k2] +=  (r1[k3]-r2[k3])*(r1[k3]-r2[k3]);    
                }
                distances[k1][k2] = sqrt(distances[k1][k2]);      
            }else{
                distances[k1][k2] = 0;
            }
        }
    }

    return distances;
}