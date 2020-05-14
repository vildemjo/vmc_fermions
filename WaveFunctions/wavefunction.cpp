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

    // Making sure that the sum of positive and negative values are correct
    if (numberOfBins % 2 == 0){
        m_numberOfBins = numberOfBins;
    }else{
        std::cout << "Number of bins must be able to be divided by 2" << std::endl;
    }
    
    // int numberOfDimensions = m_system->getNumberOfDimensions();
    std::vector<double> bins(m_numberOfBins);
    // std::vector<double> radialBins(m_numberOfBins);

    for (int n3 = 0; n3<numberOfBins; n3++){
        bins[n3] = 0;
        // binStartValue[n3] = n3*(m_densityLength/(double) m_numberOfBins);// (-m_numberOfBins/2+n3)*(densityLength/(double)numberOfBins);
    }
    // Creating an array to have control over what values the bins represent
    // m_oneBodyDensity.push_back(radialBins);

    for (int n6 = 0; n6<m_numberOfBins; n6++){
        m_oneBodyDensity.push_back(bins);
    }

}

void WaveFunction::updateOneBodyDensity(){
    /* This function collects the data needed to get the the one-body density.*/
    
    auto m_particles = m_system->getParticles();

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double binLength = m_densityLength/(double) m_numberOfBins;

    //     for (int l = 0; l < numberOfParticles; l++){

    //     auto r = m_particles[l]->getPosition();
    //     double rLength = 0;

    //     for (int j4 = 0; j4 < numberOfDimensions; j4++){
    //         rLength += r[j4]*r[j4];
    //     }
    //     rLength = sqrt(rLength);

    //     for (int j6 = 0; j6 < m_numberOfBins; j6++){
    //         if (rLength > (j6)*binLength && rLength <= (j6+1)*binLength){    
    //             m_oneBodyDensity[0][j6] += 1;
    //         }
            
    //     }
        
    // }


    for (int l = 0; l < numberOfParticles/2; l++){

        auto r = m_particles[l]->getPosition();

              
        for (int j5 = -m_numberOfBins/2; j5 < m_numberOfBins/2; j5++){
            if (j5 < 0){
                if (r[0] > (j5)*binLength && r[0] <= (j5+1)*binLength){

                    for (int j3 = -m_numberOfBins/2; j3 < m_numberOfBins/2; j3++){
                        if (j3 < 0){
                            if (r[1] > (j3)*binLength && r[1] <= (j3+1)*binLength){
                                m_oneBodyDensity[j3+m_numberOfBins/2][j5+m_numberOfBins/2] += 1;
                            }
                        }else{
                            if (r[1] >= (j3)*binLength && r[1] < (j3+1)*binLength){
                                m_oneBodyDensity[j3+m_numberOfBins/2][j5+m_numberOfBins/2] += 1;
                            }
                        }
                    }
                }
            }else{
                if (r[0] > (j5)*binLength && r[0] <= (j5+1)*binLength){
                    for (int j3 = -m_numberOfBins/2; j3 < m_numberOfBins/2; j3++){
                        if (j3 < 0){
                            if (r[1] > (j3)*binLength && r[1] <= (j3+1)*binLength){
                                m_oneBodyDensity[j3+m_numberOfBins/2][j5+m_numberOfBins/2] += 1;
                            }
                        }else{
                            if (r[1] >= (j3)*binLength && r[1] < (j3+1)*binLength){
                                m_oneBodyDensity[j3+m_numberOfBins/2][j5+m_numberOfBins/2] += 1;
                            }
                        }
                    }
                }
            }
            
        }
        //         if (j5 < 0 && j3 < 0){
        //             if (r[0] > (j5)*binLength && r[0] <= (j5+1)*binLength && r[1] > (j3)*binLength && r[1] <= (j3+1)*binLength){
        //                 // std::cout << r[0] << " < " << (j5)*binLength << " && " << r[0] << " <= " << (j5+1)*binLength << "\n";
        //                 m_oneBodyDensity[j3+m_numberOfBins/2][j5+m_numberOfBins/2] += 1;
        //             }
        //         }if (j5 > 0 && j3 > 0){
        //             if (r[0] >= (j5)*binLength && r[0] < (j5+1)*binLength && r[1] >= (j3)*binLength && r[1] < (j3+1)*binLength){
        //                 m_oneBodyDensity[j3+m_numberOfBins/2][j5+m_numberOfBins/2] += 1;
        //             }
        //         }if (j5 > 0 && j3 < 0){
        //             if (r[0] >= (j5)*binLength && r[0] < (j5+1)*binLength && r[1] > (j3)*binLength && r[1] <= (j3+1)*binLength){
        //                 m_oneBodyDensity[j3+m_numberOfBins/2][j5+m_numberOfBins/2] += 1;
        //             }
        //         }if (j5 < 0 && j3 > 0){
        //             if (r[0] > (j5)*binLength && r[0] <= (j5+1)*binLength && r[1] >= (j3)*binLength && r[1] < (j3+1)*binLength){
        //                 m_oneBodyDensity[j3+m_numberOfBins/2][j5+m_numberOfBins/2] += 1;
        //             }
        //         }
        //     }
        // }
    
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

