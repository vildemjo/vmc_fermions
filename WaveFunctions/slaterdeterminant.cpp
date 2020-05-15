#include "slaterdeterminant.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include "InitialStates/initialstate.h"
#include <iostream>

SlaterDeterminant::SlaterDeterminant(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    m_alpha = getParameters()[0];
    m_omega = m_system->getFrequency();
}

double SlaterDeterminant::evaluate() {
    /* This function calculates the trial wavefunction. */

        return m_slaterDeterminantSpinDown*m_slaterDeterminantSpinUp;
}

double SlaterDeterminant::computeDoubleDerivative() {
    /* This function calculates double derivative of the trial wavefunction
    analytically. */

    int         numberOfParticles   = m_system->getNumberOfParticles();
    double      doubleDerivative    = 0;



        // std::cout << "until double deriv ok \n";

    for(int p0 = 0; p0 < numberOfParticles/2; p0++){
        
        int p0_ = p0 + numberOfParticles/2;
        auto phi_00_p0 = phi_00(p0);
        auto phi_00_p0_ = phi_00(p0_);

        for (int p1 = 0; p1 < numberOfParticles/2; p1++){
            if (p1 == 0){
                doubleDerivative += phi_00_double_der(p0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerivative += phi_00_double_der(p0_)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 1){
                doubleDerivative += phi_10_double_der(p0,0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerivative += phi_10_double_der(p0_,0)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 2){
                doubleDerivative += phi_10_double_der(p0,1)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerivative += phi_10_double_der(p0_,1)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 3){
                doubleDerivative += phi_20_double_der(p0,0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerivative += phi_20_double_der(p0_,0)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 4){
                doubleDerivative += phi_11_double_der(p0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerivative += phi_11_double_der(p0_)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 5){
                doubleDerivative += phi_20_double_der(p0, 1)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerivative += phi_20_double_der(p0_,1)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }
            // std::cout << doubleDerivative << " p1 = "<< p1 << " p0 = "<< p0 << "\n";
        }
    }

    // std::cout << "double deriv ok \n";

    return doubleDerivative;


}

std::vector<double> SlaterDeterminant::computeDerivative(int particleIndex){
    /* This function calculates the derivative of the wavefunction with 
    regards to one spesific particle. This is used to calcualte the drift
    force used in importance sampling. */
    
    int numberOfDimensions = m_system->getNumberOfDimensions();
    int numberOfParticles = m_system->getNumberOfParticles();

    int cor = 0;
    int pI = particleIndex;

    std::vector<double> vectorSum(numberOfDimensions, 0);


        if(particleIndex < numberOfParticles/2){
        
            for (int p3 = 0; p3<numberOfDimensions; p3++){
                for (int p1 = 0; p1 < numberOfParticles/2; p1++){
                    auto phi_00_pI = phi_00(pI);
                    if (p1 == 0){
                        vectorSum[p3] += phi_00_der(pI)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 1){
                        vectorSum[p3] += phi_10_der(pI,0)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 2){
                        vectorSum[p3] += phi_10_der(pI,1)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 3){
                        vectorSum[p3] += phi_20_der(pI,0)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 4){
                        vectorSum[p3] += phi_11_der(pI)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 5){
                        vectorSum[p3] += phi_20_der(pI,1)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }
                }
            }

        }if (particleIndex >= numberOfParticles/2){
            pI = particleIndex - numberOfParticles/2;
            cor = numberOfParticles/2;

            for (int p3 = 0; p3<numberOfDimensions; p3++){
                for (int p1 = 0; p1 < numberOfParticles/2; p1++){
                    auto phi_00_pI = phi_00(pI+cor);
            
                    if (p1 == 0){
                        vectorSum[p3] += phi_00_der(pI+cor)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 1){
                        vectorSum[p3] += phi_10_der(pI+cor,0)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 2){
                        vectorSum[p3] += phi_10_der(pI+cor,1)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 3){
                        vectorSum[p3] += phi_20_der(pI+cor,0)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 4){
                        vectorSum[p3] += phi_11_der(pI+cor)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 5){
                        vectorSum[p3] += phi_20_der(pI+cor,1)[p3]*phi_00_pI*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }
                }
            }
        }

    return vectorSum;

}

double SlaterDeterminant::computeAlphaDerivative(){
    /* This function calculates the normalized derivative of the wavefunction 
    with regards to the parameter alpha. This is used to perform optimization
    by the use of gradient descent methods.*/

    auto vectorSumSquared = calculatePositionSumSquared();
    // double wf = 0;
    // int numberOfParticles = m_system->getNumberOfParticles();

    // for (int k3 = 0; k3<numberOfParticles; k3++){
    //     wf *= phi_00(k3);
    // }


    return (-0.5)*m_omega*vectorSumSquared;//*m_slaterDeterminantSpinDown*m_slaterDeterminantSpinUp/wf;

}

double SlaterDeterminant::getDistance(){
    auto m_particles = m_system->getParticles();
    std::vector<double> r0 = {0.0,0.0};//m_particles[0]->getPosition();
    auto r1 = m_particles[0]->getPosition();

    double rLength = sqrt((r0[0]-r1[0])*(r0[0]-r1[0]) + (r0[1]-r1[1])*(r0[1]-r1[1]));

    return rLength;
}

void SlaterDeterminant::setupSlaterMatrix(){
    int numberOfParticles = m_system->getNumberOfParticles();

    std::vector <std::vector <double> > slaterMatrixSpinUp(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > slaterMatrixSpinDown(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    // Test if number of particles are 2, 6 or 12?

    // std::cout << "calc wf: " << phi_00(0) <<"\n";

    for (int m1 = 0; m1 < numberOfParticles/2; m1++){ // going through the different states (single particle wavefunctions)
        for (int m2 = 0; m2 < numberOfParticles/2; m2++){ // going through the different particles
            
            int m3 = m2 + numberOfParticles/2;
            auto phi_00_m2 = phi_00(m2);
            auto phi_00_m3 = phi_00(m3);
            
            if (m1 == 0){
                slaterMatrixSpinUp[m1][m2] = phi_00_m2;
                slaterMatrixSpinDown[m1][m2] = phi_00_m3;
            }if(m1 == 1){
                slaterMatrixSpinUp[m1][m2] = phi_10(m2, 0)*phi_00_m2;
                slaterMatrixSpinDown[m1][m2] = phi_10(m3, 0)*phi_00_m3;
            }if(m1 == 2){
                slaterMatrixSpinUp[m1][m2] = phi_10(m2,1)*phi_00_m2;
                slaterMatrixSpinDown[m1][m2] = phi_10(m3,1)*phi_00_m3;
            }if(m1 == 3){
                slaterMatrixSpinUp[m1][m2] = phi_20(m2,0)*phi_00_m2;
                slaterMatrixSpinDown[m1][m2] = phi_20(m3,0)*phi_00_m3;
            }if(m1 == 4){
                slaterMatrixSpinUp[m1][m2] = phi_11(m2)*phi_00_m2;
                slaterMatrixSpinDown[m1][m2] = phi_11(m3)*phi_00_m3;
            }if(m1 == 5){
                slaterMatrixSpinUp[m1][m2] = phi_20(m2,1)*phi_00_m2;
                slaterMatrixSpinDown[m1][m2] = phi_20(m3,1)*phi_00_m3;
            }
        }
    }

    // std::cout << "ok here \n";

    m_slaterMatrixSpinDown = slaterMatrixSpinDown;
    m_slaterMatrixSpinUp = slaterMatrixSpinUp;
}

void SlaterDeterminant::updateSlaterMatrix(int particleNumber){
    int numberOfParticles = m_system->getNumberOfParticles();
    std::vector <std::vector <double> > m_slaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    int correction = 0;
    int pN = particleNumber;
    double determinant = 0;



    if(particleNumber < numberOfParticles/2){
        m_slaterMatrix = m_slaterMatrixSpinUp;
        m_oldSlaterMatrixSpinUp = m_slaterMatrixSpinUp; // saving the "old" one to calculate the inverse
        correction = 0;
    }else{
        m_slaterMatrix = m_slaterMatrixSpinDown;
        m_oldSlaterMatrixSpinDown = m_slaterMatrixSpinDown; // saving the "old" one to calculate the inverse
        pN = particleNumber - numberOfParticles/2;
        correction = numberOfParticles/2;
    }

    // std::cout << "particle number: " << particleNumber << "\n";

    // std::cout << "pN: " << pN << "\n";
    
    // std::cout << "correction: " << correction << "\n";

    for (int m3 = 0; m3 < numberOfParticles/2; m3++){ // going through the different states (for particle pN)
        auto phi_00_pN = phi_00(pN+correction);
        // std::cout << "particle number: "<< pN+correction << "\n";
        if (m3 == 0){
            m_slaterMatrix[m3][pN] = phi_00_pN;
        }if(m3 == 1){
            m_slaterMatrix[m3][pN] = phi_10(pN+correction,0)*phi_00_pN;
        }if(m3 == 2){
            m_slaterMatrix[m3][pN] = phi_10(pN+correction,1)*phi_00_pN;
        }if(m3 == 3){
            m_slaterMatrix[m3][pN] = phi_20(pN+correction,0)*phi_00_pN;
        }if(m3 == 4){
            m_slaterMatrix[m3][pN] = phi_11(pN+correction)*phi_00_pN;
        }if(m3 == 5){
            m_slaterMatrix[m3][pN] = phi_20(pN+correction,1)*phi_00_pN;
        }
        // std::cout << m_slaterMatrix[m3][pN] << "\t";
    }

            // std::cout << "update ok \n";

        for(int i_ = 0; i_ < numberOfParticles/2; i_++)
            if (numberOfParticles == 6){
                determinant += (m_slaterMatrix[0][i_] * (m_slaterMatrix[1][(i_+1)%3] * m_slaterMatrix[2][(i_+2)%3] - m_slaterMatrix[1][(i_+2)%3] * m_slaterMatrix[2][(i_+1)%3]));
            }if (numberOfParticles == 2){
                determinant = m_slaterMatrix[0][0];
            }
    // std::cout << "In update matrix: " << "\n ---------- \n";
    // std::cout << "determinant before: "<< m_slaterDeterminantSpinDown << " \n";

    if(particleNumber < numberOfParticles/2){
        m_slaterMatrixSpinUp = m_slaterMatrix;
        m_slaterDeterminantSpinUp = determinant;
    }else{
        m_slaterMatrixSpinDown = m_slaterMatrix;
        m_slaterDeterminantSpinDown = determinant;
    }
        // std::cout << "determinant after: "<< m_slaterDeterminantSpinDown << " \n  -------\n";

}

void SlaterDeterminant::updateSlaterRelatedThings(int particleNumber){
    updateSlaterMatrix(particleNumber);
    updateInverseSlaterMatrix(particleNumber);
}

void SlaterDeterminant::setupSlaterRelatedThings(){
    setupSlaterMatrix();
    // std::cout << "Slater matrix ok" << std::endl;
    calculateInverseSlaterMatrix();
    // std::cout << "Slater matrix inverse ok" << std::endl;
}

void SlaterDeterminant::calculateInverseSlaterMatrix(){
    // Need to find out how to calculate this.
    // Should I do this in the initiation of this class?
    // Should I have another class for the slater determinant parts?
    int numberOfParticles = m_system->getNumberOfParticles();
    int i_, j_;

    std::vector <std::vector <double> > m_slaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > m_inverseSlaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
	float determinant = 0;
	
  

    for (int s = 0; s < 2; s++){

        if (s == 0){
            m_slaterMatrix = m_slaterMatrixSpinUp;
        }if(s == 1){
            m_slaterMatrix = m_slaterMatrixSpinDown;
        }


        // finding determinant
        for(i_ = 0; i_ < numberOfParticles/2; i_++)
            if (numberOfParticles == 6){
                determinant += (m_slaterMatrix[0][i_] * (m_slaterMatrix[1][(i_+1)%3] * m_slaterMatrix[2][(i_+2)%3] - m_slaterMatrix[1][(i_+2)%3] * m_slaterMatrix[2][(i_+1)%3]));
            }if (numberOfParticles == 2){
                determinant = m_slaterMatrix[0][0];
            }
        // std::cout << "\n ------- \n determinant: " << determinant << "\n ---------- \n";

        for(i_ = 0; i_ < numberOfParticles/2; i_++){
            // std::cout << "\n";
            for(j_ = 0; j_ < numberOfParticles/2; j_++){
                if (numberOfParticles == 6){
                    m_inverseSlaterMatrix[i_][j_] = ((m_slaterMatrix[(j_+1)%3][(i_+1)%3] * m_slaterMatrix[(j_+2)%3][(i_+2)%3]) - (m_slaterMatrix[(j_+1)%3][(i_+2)%3] * m_slaterMatrix[(j_+2)%3][(i_+1)%3]))/ determinant;
                // std::cout << m_inverseSlaterMatrix[i_][j_] << "\t";
                }if (numberOfParticles == 2){
                    m_inverseSlaterMatrix[i_][j_] = 1/m_slaterMatrix[i_][j_];
                }
            }
        }
    
            // std::cout << "In inverse matrix: " << "\n ---------- \n";
            // std::cout << "determinant before: "<< m_slaterDeterminantSpinDown << " \n";

        if (s == 0){
            m_inverseSlaterMatrixSpinUp = m_inverseSlaterMatrix;
            m_slaterDeterminantSpinUp = determinant;
            m_oldInverseSlaterMatrixSpinUp = m_inverseSlaterMatrixSpinUp;
        }else{
            m_inverseSlaterMatrixSpinDown = m_inverseSlaterMatrix;
            m_slaterDeterminantSpinDown = determinant;
            m_oldInverseSlaterMatrixSpinDown = m_inverseSlaterMatrixSpinDown;
        }
    
            // std::cout << "determinant after: "<< m_slaterDeterminantSpinDown << " \n ------------ \n";


    }


}

void SlaterDeterminant::updateInverseSlaterMatrix(int particleNumber){
    // Need to find out how to calculate this.
    // Should I do this in the initiation of this class?
    // Should I have another class for the slater determinant parts?
    int numberOfParticles = m_system->getNumberOfParticles();
    int i_, j_;

    double correction = 0;
    int pN = particleNumber;
	
    // if(particleNumber < numberOfParticles/2){
    //         m_oldSlaterMatrixSpinUp = m_slaterMatrixSpinUp; // saving the "old" one to calculate the inverse
    //         m_oldInverseSlaterMatrixSpinUp = m_inverseSlaterMatrixSpinUp;
    //         correction = 0;
    //         for(i_ = 0; i_ < numberOfParticles/2; i_++){
    //         // std::cout << "\n";
    //         for(j_ = 0; j_ < numberOfParticles/2; j_++){
    //             if (numberOfParticles == 6){
    //                 m_inverseSlaterMatrixSpinUp[i_][j_] = ((m_slaterMatrixSpinUp[(j_+1)%3][(i_+1)%3] * m_slaterMatrixSpinUp[(j_+2)%3][(i_+2)%3]) 
    //                 - (m_slaterMatrixSpinUp[(j_+1)%3][(i_+2)%3] * m_slaterMatrixSpinUp[(j_+2)%3][(i_+1)%3]))/ m_slaterDeterminantSpinUp;
    //             // std::cout << inverseSlaterMatrix[i_][j_] << "\t";
    //             }if (numberOfParticles == 2){
    //                 m_inverseSlaterMatrixSpinUp[i_][j_] = 1/m_slaterMatrixSpinUp[i_][j_];
    //             }if(numberOfParticles != 2 && numberOfParticles != 6){
    //                 std::cout << "ERROR: invalid number of particles \n";
    //             }
    //         }
    //     }
            
    //     }else{
    //         m_oldSlaterMatrixSpinDown = m_slaterMatrixSpinDown; // saving the "old" one to calculate the inverse
    //         pN = particleNumber - numberOfParticles/2;
    //         correction = numberOfParticles/2;
    //         m_oldInverseSlaterMatrixSpinDown = m_inverseSlaterMatrixSpinDown;

    //         for(i_ = 0; i_ < numberOfParticles/2; i_++){
    //             // std::cout << "\n";
    //             for(j_ = 0; j_ < numberOfParticles/2; j_++){
    //                 if (numberOfParticles == 6){
    //                     m_inverseSlaterMatrixSpinDown[i_][j_] = ((m_slaterMatrixSpinDown[(j_+1)%3][(i_+1)%3] * m_slaterMatrixSpinDown[(j_+2)%3][(i_+2)%3]) 
    //                     - (m_slaterMatrixSpinDown[(j_+1)%3][(i_+2)%3] * m_slaterMatrixSpinDown[(j_+2)%3][(i_+1)%3]))/ m_slaterDeterminantSpinDown;
    //                 // std::cout << inverseSlaterMatrix[i_][j_] << "\t";
    //                 }if (numberOfParticles == 2){
    //                     m_inverseSlaterMatrixSpinDown[i_][j_] = 1/m_slaterMatrixSpinDown[i_][j_];
    //                 }if(numberOfParticles != 2 && numberOfParticles != 6){
    //                     std::cout << "ERROR: invalid number of particles \n";
    //                 }
    //             }
    //         }
    //     }   

    double determinant = 0;
    std::vector <std::vector <double> > slaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > inverseSlaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));


    if(particleNumber < numberOfParticles/2){
            slaterMatrix = m_slaterMatrixSpinUp;
            m_oldSlaterMatrixSpinUp = m_slaterMatrixSpinUp; // saving the "old" one to calculate the inverse
            correction = 0;
            determinant = m_slaterDeterminantSpinUp;
            
        }else{
            slaterMatrix = m_slaterMatrixSpinDown;
            m_oldSlaterMatrixSpinDown = m_slaterMatrixSpinDown; // saving the "old" one to calculate the inverse
            pN = particleNumber - numberOfParticles/2;
            correction = numberOfParticles/2;
            determinant = m_slaterDeterminantSpinDown;
        }


        for(i_ = 0; i_ < numberOfParticles/2; i_++){
            // std::cout << "\n";
            for(j_ = 0; j_ < numberOfParticles/2; j_++){
                if (numberOfParticles == 6){
                    inverseSlaterMatrix[i_][j_] = ((slaterMatrix[(j_+1)%3][(i_+1)%3] * slaterMatrix[(j_+2)%3][(i_+2)%3]) - (slaterMatrix[(j_+1)%3][(i_+2)%3] * slaterMatrix[(j_+2)%3][(i_+1)%3]))/ determinant;
                // std::cout << inverseSlaterMatrix[i_][j_] << "\t";
                }if (numberOfParticles == 2){
                    inverseSlaterMatrix[i_][j_] = 1/slaterMatrix[i_][j_];
                }if(numberOfParticles != 2 && numberOfParticles != 6){
                    std::cout << "ERROR: invalid number of particles \n";
                }
            }
        }
    
            // std::cout << "In inverse matrix: " << "\n ---------- \n";
            // std::cout << "determinant before: "<< m_slaterDeterminantSpinDown << " \n";

        if(particleNumber < numberOfParticles/2){
            m_slaterMatrixSpinUp = slaterMatrix;
            m_inverseSlaterMatrixSpinUp = inverseSlaterMatrix;
            m_slaterDeterminantSpinUp = determinant;
        }else{
            m_slaterMatrixSpinDown = slaterMatrix;
            m_inverseSlaterMatrixSpinDown = inverseSlaterMatrix;
            m_slaterDeterminantSpinDown = determinant;
        }
    
            // std::cout << "determinant after: "<< m_slaterDeterminantSpinDown << " \n ------------ \n";


}


double SlaterDeterminant::computeRatio(double oldWaveFunction, double newWaveFunction){
    m_metropolisRatio = newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction);
    return m_metropolisRatio;
}


std::vector<double> SlaterDeterminant::computeQuantumForce(int particleIndex, bool oldOrNew){
    /* This function calculates the quantum force/drift force with is used for importance
        sampling. The quantum force is given by the derivative of the wavefunction. */
        
    double R_inv = 1;

    if (oldOrNew == true){
        R_inv = 1;
    }if(oldOrNew == false){
        R_inv = 1;// /m_metropolisRatio;
    }

  
    auto derivative = computeDerivative(particleIndex);

    for (int m=0;m<m_system->getNumberOfDimensions();m++){
        derivative[m] *= 2;
        // std::cout << derivative[m] << "\n";
    }


    return derivative;
}

// -----------------------------------
// Basis:
// -----------------------------------

double SlaterDeterminant::phi_00(int particleNumber){
    double A = 1;
    double factor = 0;
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double rSquared = 0;

    auto m_alpha = m_system->getWaveFunction()->getParameters()[0];
   
    auto r = m_system->getParticles()[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        rSquared += r[b0]*r[b0];
    }

    factor = -0.5*m_alpha*m_omega*rSquared;

    return A*exp(factor);
}

double SlaterDeterminant::phi_10(int particleNumber, int dimension){
    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    return 2*sqrt(m_omega)*r[dimension];
}

double SlaterDeterminant::phi_20(int particleNumber, int dimension){
    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    return 4*m_omega*r[dimension]*r[dimension]-2;
}

double SlaterDeterminant::phi_11(int particleNumber){
    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    return 4*m_omega*r[0]*r[1];
}

// Derivatives:

std::vector<double> SlaterDeterminant::phi_00_der(int particleNumber){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    std::vector<double> derivativeVector(numberOfDimensions);

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        derivativeVector[b0] = -m_alpha*m_omega*r[b0];
    }

    return derivativeVector;
}

std::vector<double> SlaterDeterminant::phi_10_der(int particleNumber, int dimension){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    std::vector<double> derivativeVector(numberOfDimensions);

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();
    int dim_other = 0;

    if (dimension == 0){
        dim_other = 1;
    }if(dimension == 1){
        dim_other = 0;}


    derivativeVector[dimension] = -2*sqrt(m_omega)*(m_alpha*m_omega*r[dimension]*r[dimension]-1);
    derivativeVector[dim_other] = -2*m_alpha*sqrt(m_omega)*m_omega*r[dimension]*r[dim_other];

    return derivativeVector;
}

std::vector<double> SlaterDeterminant::phi_20_der(int particleNumber, int dimension){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    std::vector<double> derivativeVector(numberOfDimensions);

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();
    int dim_other = 0;

    if (dimension == 0){
        dim_other = 1;
    }else{
        dim_other = 0;}

    derivativeVector[dimension] = -2*(2*m_alpha*m_omega*m_omega*r[dimension]*r[dimension]*r[dimension]
                                -m_alpha*m_omega*r[dimension]-4*m_omega*r[dimension]);
    derivativeVector[dim_other] = -2*(2*m_alpha*m_omega*m_omega*r[dimension]*r[dimension]*r[dim_other]
                                -m_alpha*m_omega*r[dim_other]);
    
    return derivativeVector;
}

std::vector<double> SlaterDeterminant::phi_11_der(int particleNumber){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    std::vector<double> derivativeVector(numberOfDimensions);

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    std::vector<double> other_dim = {1,0};

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        derivativeVector[b0] = -4*m_omega*r[other_dim[b0]]*(m_alpha*m_omega*r[b0]*r[b0]-1);
    }

    return derivativeVector;
}

// double derivatives

double SlaterDeterminant::phi_00_double_der(int particleNumber){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double rSquared = 0;

    auto m_alpha = m_system->getWaveFunction()->getParameters()[0];
    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        rSquared += r[b0]*r[b0];
    }


    return m_alpha*m_alpha*m_omega*m_omega*rSquared-2*m_alpha*m_omega;
}

double SlaterDeterminant::phi_10_double_der(int particleNumber, int dimension){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double rSquared = 0;

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        rSquared += r[b0]*r[b0];
    }

    return 2*m_alpha*sqrt(m_omega)*m_omega*r[dimension]*(m_alpha*m_omega*rSquared-4);// /(2*sqrt(m_omega)*r[dimension]);
}

double SlaterDeterminant::phi_20_double_der(int particleNumber, int dimension){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double rSquared = 0;

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        rSquared += r[b0]*r[b0];
    }

    return 2*m_omega*(m_alpha*m_alpha*m_omega*(2*m_omega*r[dimension]*r[dimension]-1)*rSquared + m_alpha*(2-12*m_omega*r[dimension]*r[dimension])+4);
}

double SlaterDeterminant::phi_11_double_der(int particleNumber){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double rSquared = 0;

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        rSquared += r[b0]*r[b0];
    }

    return 4*m_alpha*m_omega*m_omega*r[0]*r[1]*(m_alpha*m_omega*rSquared-6);
}

