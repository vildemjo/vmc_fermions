#include "gaussiandistribution.h"
#include "wavefunction.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"
#include <cmath>

using std::cout;
using std::endl;

GaussianDistribution::GaussianDistribution(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles, 
                             double     stepLength)  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;
    m_stepLength         = stepLength;

    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState();
 
}

void GaussianDistribution::setupInitialState() {
    /* This function initializes the particles with a random position 
    according to a Gaussian distribution. This initialization is used 
    when importance sampling is used in the metropolis step because the 
    gaussain distribution is also used to find the new position.*/

    bool positionCheck = false;
    while ( positionCheck == false){

        for (int m3=0; m3 < m_numberOfParticles; m3++) {
            std::vector<double> position = std::vector<double>();

            for (int m4=0; m4 < m_numberOfDimensions; m4++) {
                position.push_back(Random::nextGaussian(0,1)*sqrt(m_stepLength));
            }

            m_particles.push_back(new Particle());
            m_particles.at(m3)->setNumberOfDimensions(m_numberOfDimensions);
            m_particles.at(m3)->setPosition(position);
            m_particles.at(m3)->setParticleIndex(m3);
        }

        positionCheck = m_system->getWaveFunction()->getDistanceCheck(m_particles);

    }

}
