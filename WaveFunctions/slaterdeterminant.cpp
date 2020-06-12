#include "slaterdeterminant.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include "InitialStates/initialstate.h"
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

SlaterDeterminant::SlaterDeterminant(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    m_alpha = getParameters()[0];
    m_omega = m_system->getFrequency();
    // std::cout << m_omega;
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

    for(int p0 = 0; p0 < numberOfParticles/2; p0++){ // Loop over the different particles
        
        int p0_ = p0 + numberOfParticles/2; // Particle number for the spin down particles
        
        // Calculate first the exponential term of the single-particle wavefunction 
        // that is the same for all of them
        auto phi_00_p0 = phi_00(p0); 
        auto phi_00_p0_ = phi_00(p0_);

        for (int p1 = 0; p1 < numberOfParticles/2; p1++){ // Loop over the different states
            if (p1 == 0){ // test to get the correct single-particle function
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
        }
    }

    return doubleDerivative;

}

std::vector<double> SlaterDeterminant::computeDerivative(int particleIndex){
    /* This function calculates the derivative of the wavefunction with 
    regards to one spesific particle. This is used to calculate the drift
    force used in importance sampling. */
    
    int numberOfDimensions = m_system->getNumberOfDimensions();
    int numberOfParticles = m_system->getNumberOfParticles();

    int cor = 0;
    int pI = particleIndex;

    std::vector<double> vectorSum(numberOfDimensions, 0);

    // Using test to see if it is a spin up or a spin down particle that is being moved.
    // It is nesessary because the Slater deterimant is seperated into a spin up 
    // determinant and a spin down determinant.
        if(particleIndex < numberOfParticles/2){
            for (int p3 = 0; p3<numberOfDimensions; p3++){
                for (int p1 = 0; p1 < numberOfParticles/2; p1++){ // Loop over the different states
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
            cor = numberOfParticles/2; // Using a correction number to get the correct particle number for the
                                       // calculation of the gradient of the single particle wave function. This
                                       // is necessary because the spin down electrons are the upper half of the 
                                       // particles, but they have their own Slater determinant.
            for (int p3 = 0; p3<numberOfDimensions; p3++){
                for (int p1 = 0; p1 < numberOfParticles/2; p1++){ // Loop over the different states
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

    return (-0.5)*m_omega*vectorSumSquared;

}

double SlaterDeterminant::getDistance(){
    /* This function calculates distance between particle one and two. This 
    function is used to evaluate the mean distance in a two-particle system.*/

    auto m_particles = m_system->getParticles();
    
    auto r0 = m_particles[0]->getPosition();
    auto r1 = m_particles[1]->getPosition();

    double rLength = sqrt((r0[0]-r1[0])*(r0[0]-r1[0]) + (r0[1]-r1[1])*(r0[1]-r1[1]));

    return rLength;
}

void SlaterDeterminant::setupSlaterMatrix(){
    /* This function sets up the Slater determinants. One for the spin up electrons
     (first half) and one for the spin down electrons (last half).*/

    int numberOfParticles = m_system->getNumberOfParticles();

    std::vector <std::vector <double> > slaterMatrixSpinUp(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > slaterMatrixSpinDown(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));


    for (int m1 = 0; m1 < numberOfParticles/2; m1++){ // going through the different states (single particle wavefunctions)
        for (int m2 = 0; m2 < numberOfParticles/2; m2++){ // going through the different particles
            
            int m3 = m2 + numberOfParticles/2; // The particle number of the spin down particles.
            
            // Calculate first the exponential term of the single-particle wavefunction 
            // that is the same for all of them
            auto phi_00_m2 = phi_00(m2); 
            auto phi_00_m3 = phi_00(m3);
            
            if (m1 == 0){
                slaterMatrixSpinUp[m1][m2] = phi_00_m2;
                slaterMatrixSpinDown[m1][m2] = phi_00_m3;
            }if(m1 == 1){ // test to get the correct Hermite polynomial
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

    m_slaterMatrixSpinDown = slaterMatrixSpinDown;
    m_slaterMatrixSpinUp = slaterMatrixSpinUp;
}

void SlaterDeterminant::updateSlaterMatrix(int particleNumber){
    /* This function updates the row in the Slater determinant 
    that represents the particle that was moved.*/

    int numberOfParticles = m_system->getNumberOfParticles();
    std::vector <std::vector <double> > m_slaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    int correction = 0;
    int pN = particleNumber;

    // A test to find out if it is the spin up or the spin down Slater determinant
    // that is being updated
    if(particleNumber < numberOfParticles/2){
        m_slaterMatrix = m_slaterMatrixSpinUp;
        correction = 0;
    }else{
        m_slaterMatrix = m_slaterMatrixSpinDown;
        pN = particleNumber - numberOfParticles/2;
        correction = numberOfParticles/2; // A correction term i used to get the 
                                          // correct particel number if it is a spin 
                                          // down particle that is moved.
    }


    for (int m3 = 0; m3 < numberOfParticles/2; m3++){ // going through the different states (for particle pN)
        auto phi_00_pN = phi_00(pN+correction);
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
    }

    if(particleNumber < numberOfParticles/2){
        m_slaterMatrixSpinUp = m_slaterMatrix;
    }else{
        m_slaterMatrixSpinDown = m_slaterMatrix;
    }

}

void SlaterDeterminant::updateSlaterRelatedThings(int particleNumber){
    updateSlaterMatrix(particleNumber);
    updateInverseSlaterMatrix();
}

void SlaterDeterminant::setupSlaterRelatedThings(){
    setupSlaterMatrix();
    calculateInverseSlaterMatrix();
}

void SlaterDeterminant::calculateInverseSlaterMatrix(){
    /* This function calculates the inverse of the Slater determinants and
    evaluate their determinant. */
    const int numberOfParticles = m_system->getNumberOfParticles();

    std::vector <std::vector <double> > m_slaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > m_inverseSlaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));

    MatrixXd m(numberOfParticles/2, numberOfParticles/2);
    MatrixXd m_inv(numberOfParticles/2, numberOfParticles/2);


    for (int s = 0; s < 2; s++){

        if (s == 0){
            m_slaterMatrix = m_slaterMatrixSpinUp;
        }if(s == 1){
            m_slaterMatrix = m_slaterMatrixSpinDown;
        }

        
        // Coverting into Eigen library matrix
        for (int i_i = 0; i_i < numberOfParticles/2; i_i++){
            for (int j_i = 0; j_i < numberOfParticles/2; j_i++){
                m(i_i,j_i) = m_slaterMatrix[i_i][j_i];
            }
        }

        // Inverse:
        m_inv = m.inverse();

        // Coverting back to std::vector
        for (int i_i = 0; i_i < numberOfParticles/2; i_i++){
            for (int j_i = 0; j_i < numberOfParticles/2; j_i++){
                m_inverseSlaterMatrix[i_i][j_i] = m_inv(i_i,j_i);
            }
        }


        if (s == 0){
            m_inverseSlaterMatrixSpinUp = m_inverseSlaterMatrix;
            m_slaterDeterminantSpinUp = m.determinant();
            m_oldInverseSlaterMatrixSpinUp = m_inverseSlaterMatrixSpinUp;
        }else{
            m_inverseSlaterMatrixSpinDown = m_inverseSlaterMatrix;
            m_slaterDeterminantSpinDown = m.determinant();
            m_oldInverseSlaterMatrixSpinDown = m_inverseSlaterMatrixSpinDown;
        }
    
    }


}

void SlaterDeterminant::updateInverseSlaterMatrix(){
    /* This function should ideally just update the inverse Slater determinant 
    using the old inverse Slater determinant, but for now it just calculates it. */

    calculateInverseSlaterMatrix();

}


double SlaterDeterminant::computeRatio(double oldWaveFunction, double newWaveFunction){
    m_metropolisRatio = newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction);
    return m_metropolisRatio;
}


std::vector<double> SlaterDeterminant::computeQuantumForce(int particleIndex){
    /* This function calculates the quantum force/drift force with is used for importance
        sampling. The quantum force is given by the derivative of the wavefunction. */
  
    auto derivative = computeDerivative(particleIndex);

    for (int m=0;m<m_system->getNumberOfDimensions();m++){
        derivative[m] *= 2;
    }

    return derivative;
}

// -----------------------------------
// Basis:
// -----------------------------------

double SlaterDeterminant::phi_00(int particleNumber){
    
    double factor = 0;
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double rSquared = 0;

    auto m_alpha = m_system->getWaveFunction()->getParameters()[0];
   
    auto r = m_system->getParticles()[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        rSquared += r[b0]*r[b0];
    }

    factor = -0.5*m_alpha*m_omega*rSquared;

    return exp(factor);
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

