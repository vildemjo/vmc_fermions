#include "wavefunction.h"
#include "system.h"
#include "particle.h"
#include <cmath>
#include <iostream>


WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

void WaveFunction::setOneBodyDensityBins(int numberOfBins, double densityLength){
    m_densityLength = densityLength;

    std::vector<double> bins(m_numberOfBins);
    double binStartValue = 0;
    m_oneBodyDensityRadial = std::vector<std::vector<double> >(2, std::vector<double>(numberOfBins/2,0)); 


    for (int n4 = 0; n4 < numberOfBins/2; n4++){
        binStartValue = n4*(m_densityLength/(double) m_numberOfBins);
        m_oneBodyDensityRadial[0][n4] = binStartValue;     // Creating an array to have control over what values the bins represent
        m_oneBodyDensityRadial[1][n4] = 0;
    }


}

void WaveFunction::updateOneBodyDensity(){
    /* This function collects the data needed to get the radial one-body density.*/
    
    auto m_particles = m_system->getParticles();

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double binLength = m_densityLength/(double) m_numberOfBins;

    for (int l = 0; l < numberOfParticles; l++){

        auto r = m_particles[l]->getPosition();

        double rLength = 0;

        for (int j4 = 0; j4 < numberOfDimensions; j4++){
            rLength += r[j4]*r[j4];
        }
        rLength = sqrt(rLength);

        for (int j6 = 0; j6 < m_numberOfBins/2; j6++){
            if (rLength >= (j6)*binLength && rLength < (j6+1)*binLength){    
                m_oneBodyDensityRadial[1][j6] += 1;
            }
        }
    
    }

}


double WaveFunction::calculatePositionSumSquared(){
    /* This function calcualtes the positions squared and summed over all
    particles. This is used in several of the functions in the SlaterDeterminant 
    class. The other have it's own to include beta.*/

    double rSum = 0.0;
    auto m_particles = m_system->getParticles();
    
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    for(int i1=0; i1<numberOfParticles; i1++){
        auto r = m_particles[i1]->getPosition();
        
        for(int n1=0; n1<numberOfDimensions; n1++){
            rSum += r[n1]*r[n1];
        }
    }

    return rSum;
}

