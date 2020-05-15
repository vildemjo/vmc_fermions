#include "randomuniform.h"
#include "wavefunction.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"
#include <cmath>

using std::cout;
using std::endl;

RandomUniform::RandomUniform(System*    system,
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

void RandomUniform::setupInitialState() {
    /* This function initializes the particles with a random position 
    according to a uniform distribution.*/



    for (int m3=0; m3 < m_numberOfParticles; m3++) {
        std::vector<double> position = std::vector<double>();

        for (int m4=0; m4 < m_numberOfDimensions; m4++) {
            position.push_back(m_stepLength*(m_system->getRandomEngine()->nextDouble()-0.5));
        }

        m_particles.push_back(new Particle());
        m_particles.at(m3)->setNumberOfDimensions(m_numberOfDimensions);
        m_particles.at(m3)->setPosition(position);
        m_particles.at(m3)->setParticleIndex(m3);
    }

  
}
