#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include "InitialStates/initialstate.h"
#include <iostream>

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate() {
    /* This function calculates the trial wavefunction. */

    auto rSum2 = calculatePositionSumSquared();
    double norm = 1;
    double omega = m_system->getFrequency();

    return norm*exp(-0.5*m_parameters[0]*omega*rSum2);

}

double SimpleGaussian::computeDoubleDerivative() {
    /* This function calculates double derivative of the trial wavefunction
    analytically. */

    int         numberOfParticles   = m_system->getNumberOfParticles();
    int         numberOfDimensions  = m_system->getNumberOfDimensions();
    double      omega               = m_system->getFrequency();

    auto rSum2 = calculatePositionSumSquared(); 
    
    return (-omega*m_parameters[0]*numberOfParticles*numberOfDimensions 
                        + omega*omega*m_parameters[0]*m_parameters[0]*rSum2);
}

std::vector<double> SimpleGaussian::computeDerivative(int particleIndex){
    /* This function calculates the derivative of the wavefunction with 
    regards to one spesific particle. This is used to calcualte the drift
    force used in importance sampling. */
    
    int numberOfDimensions = m_system->getNumberOfDimensions();
    auto m_particles = m_system->getParticles();

    std::vector<double> vectorSum(numberOfDimensions);
    double omega = m_system->getFrequency();

    auto r = m_particles[particleIndex]->getPosition();

    for (int n8=0; n8<numberOfDimensions; n8++){
        vectorSum[n8] = -omega*getParameters()[0]*r[n8];
    }

    return vectorSum;

}

double SimpleGaussian::computeAlphaDerivative(){
    /* This function calculates the normalized derivative of the wavefunction 
    with regards to the parameter alpha. This is used to perform optimization
    by the use of gradient descent methods.*/

    auto vectorSumSquared = calculatePositionSumSquared();
    double omega = m_system->getFrequency();

    return (-0.5)*omega*vectorSumSquared;

}

double SimpleGaussian::getDistance(){
    auto m_particles = m_system->getParticles();
    auto r0 = m_particles[0]->getPosition();
    auto r1 = m_particles[1]->getPosition();

    double rLength = sqrt((r0[0]-r1[0])*(r0[0]-r1[0]) + (r0[1]-r1[1])*(r0[1]-r1[1]));

    return rLength;
}