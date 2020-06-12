#include "slaterdeterminantinteraction.h"
#include "interactionharmonicoscillator.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include <iostream>
#include "InitialStates/initialstate.h"
#include <Eigen/Dense>

using namespace Eigen;



SlaterDeterminantInteraction::SlaterDeterminantInteraction(System* system, double alpha, double beta) :
/* This is the class used to model interaction particles. The parts of the functions that deals with
the Slater determinant (not Jastrow factor) is explained more thoroughly in the SlaterDeterminant class.*/


        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_omega = m_system->getFrequency(); 
}

double SlaterDeterminantInteraction::evaluate() {
    /* This function calculated the trial wavefuntion in the case 
    with interaction. */

    double interactionPart = evaluateCorrelationPart();

    return m_slaterDeterminantSpinUp*m_slaterDeterminantSpinDown*interactionPart;

}

double SlaterDeterminantInteraction::evaluateCorrelationPart() {
    /* This function calculates the correlation/interaction part of the 
    trial wavefunction i.e. the Jastrow factor. */

    double correlationPart = 0;
  
    auto a = m_system->getSpinFactor(); // This is a matrix of spin factors which are 1 if the 
                                        // spins are antiparallel and 1/3 if the spins are parallel

    auto distances = calculateInterparticleDistances();

    int numberOfParticles = m_system->getNumberOfParticles();

    for (int j7 = 0; j7<numberOfParticles-1; j7++){ 
        for (int j8 = j7+1; j8<numberOfParticles; j8++){ // Double loop which counts only a pair of 
                                                         // electrons only once

            correlationPart += a[j7][j8]*distances[j7][j8]/(1+m_parameters[1]*distances[j7][j8]);
        }
    }
    
    return exp(correlationPart);

}

double SlaterDeterminantInteraction::computeDoubleDerivative() {
    /* This function calculates the double derivative of the Slater determinant part
    of the wave function and add the interaction part to that. */

    int numberOfParticles = m_system->getNumberOfParticles();
    double doubleDerPsi_SD = 0;
    double interactionPart = 0;

    interactionPart = computeInteractionPartOfDoubleDerivative();

    for(int p0 = 0; p0 < numberOfParticles/2; p0++){
        
        int p0_ = p0 + numberOfParticles/2;
        auto phi_00_p0 = phi_00(p0);
        auto phi_00_p0_ = phi_00(p0_);

        for (int p1 = 0; p1 < numberOfParticles/2; p1++){
            if (p1 == 0){
                doubleDerPsi_SD += phi_00_double_der(p0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerPsi_SD += phi_00_double_der(p0_)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 1){
                doubleDerPsi_SD += phi_10_double_der(p0,0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerPsi_SD += phi_10_double_der(p0_,0)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 2){
                doubleDerPsi_SD += phi_10_double_der(p0,1)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerPsi_SD += phi_10_double_der(p0_,1)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 3){
                doubleDerPsi_SD += phi_20_double_der(p0,0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerPsi_SD += phi_20_double_der(p0_,0)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 4){
                doubleDerPsi_SD += phi_11_double_der(p0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerPsi_SD += phi_11_double_der(p0_)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }if (p1 == 5){
                doubleDerPsi_SD += phi_20_double_der(p0, 1)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p0][p1];
                doubleDerPsi_SD += phi_20_double_der(p0_,1)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p0][p1];
            }
        }
    }

    return doubleDerPsi_SD + interactionPart;
}

std::vector<double> SlaterDeterminantInteraction::computeDerivative(int particleIndex){
    /* This function calculates the derivative of the trial wavefunction with 
    interaction with regards to one particle. */


    int                     numberOfDimensions          = m_system->getNumberOfDimensions();
    std::vector <double>    vectorWithInteraction       (numberOfDimensions,0);

    // The derivative of the interaction part is calculated in
    // another function, because the same expression is used to 
    // evaluate the double derivative.
    auto derivativePsi_C   = computeDerivativeOfu(particleIndex);

    auto derivativePsi_SD = computeDerivativePsi_SD(particleIndex);

    for (int n8=0; n8<numberOfDimensions; n8++){
        vectorWithInteraction[n8] = derivativePsi_SD[n8] + derivativePsi_C[n8];
    }

    return vectorWithInteraction;
}

double SlaterDeterminantInteraction::computeAlphaDerivative(){
    /* This function calculates the normalized derivative of the wavefunction 
    with regards to the parameter alpha. This is used to perform optimization
    by the use of gradient descent methods.*/

    auto vectorSumSquared =  calculatePositionSumSquared();

    return (-0.5)*m_omega*vectorSumSquared;

}

double SlaterDeterminantInteraction::computeBetaDerivative(){
    /* This function calculates the normalized derivative of the wavefunction 
    with regards to the parameter beta. This is used to perform optimization
    by the use of gradient descent methods.*/

    double factorSum = 0;
    int numberOfParticles = m_system->getNumberOfParticles();
    auto distances = calculateInterparticleDistances();
    auto a = m_system->getSpinFactor(); // This is a matrix of spin factors which are 1 if the 
                                        // spins are antiparallel and 1/3 if the spins are parallel



    for (int j9 = 0; j9<numberOfParticles-1; j9++){
        for (int j10 = j9+1; j10<numberOfParticles; j10++){

            auto rLength = distances[j9][j10];
            double factor = (1+m_parameters[1]*rLength);
            factorSum += -a[j9][j10]*rLength*rLength/(factor*factor);
        }
    }

    return factorSum; 

}

double SlaterDeterminantInteraction::computeInteractionPartOfDoubleDerivative(){
    /* This function calculates the interaction part of the double derivative
    when it is evaluated analytically. */

    // These terms represents the three terms the interaction is involved in
    // in the analytical expression of the double derivative
    double      firstTerm = 0; double secondTerm = 0; double thirdTerm = 0;

    int         numberOfParticles = m_system->getNumberOfParticles();
    int         numberOfDimensions = m_system->getNumberOfDimensions();
    auto        m_particles = m_system->getParticles();
    std::vector<double> uPart(numberOfDimensions+1);


    for(int l5 = 0; l5 < numberOfParticles; l5++){ 

        // The first term also includes the non-interating part of
        // the trial wavefunction. This is evaluated here.
        auto derivativePsi_SD = computeDerivativePsi_SD(l5);

        // uPart[0], uPart[1] and uPart[2] is the derivative of
        // the interaction part of the trial wavefunction which is
        // a vector. uPart[3] is the last term which is a part of the
        // double derivative of the interaction part of the trial
        // wavefunction. 
        // See computeDerivativeOfu for more details.

         uPart = computeDerivativeOfu(l5); // first comes derivativePsi_C and the last value is part of doubleDerPsi_C


        for (int l4 = 0; l4<numberOfDimensions; l4++){
            firstTerm += derivativePsi_SD[l4]*uPart[l4];
            secondTerm += uPart[l4]*uPart[l4];
        }
        thirdTerm += uPart[numberOfDimensions];
        
    }

    return 2*firstTerm + secondTerm + thirdTerm; // second and third term is together doubleDerPsi_C
}

std::vector <double> SlaterDeterminantInteraction::computeDerivativeOfu(int particleNumber){
    /* This function calcualtes the derivative and double derivative of the interation term 
    in the trial wavefunction which are used in calcualting the double derivative and the 
    derivative with regards to one particle. */

    int                     numberOfParticles       = m_system->getNumberOfParticles();
    int                     numberOfDimensions      = m_system->getNumberOfDimensions();
    double                  uDerivative             = 0;
    auto                    a                       = m_system->getSpinFactor();
    double                  uDoubleDerivative       = 0;
    double                  uTotalDoubleDerivative  = 0;
    std::vector <double>    uTotalDerivative        (numberOfDimensions);
    std::vector <double>    uAll                    (numberOfDimensions+1);


    auto m_particles = m_system->getParticles();
    auto ri          = m_particles[particleNumber]->getPosition();   

    auto difference = calculateInterparticleDistances();


    // The elements necessary to calculate doubleDerPsi_C/PsiC is calculated here    

    for (int l1 = 0; l1 < numberOfParticles; l1++){

        if (particleNumber != l1){
            auto rj         = m_particles[l1]->getPosition();        
            auto rLength    = difference[particleNumber][l1];                           

            auto factor = (1+m_parameters[1]*rLength);
            uDerivative = a[l1][particleNumber]/(factor*factor);
            uDoubleDerivative = -2*a[l1][particleNumber]*m_parameters[1]/(factor*factor*factor);
       

            for (int l3 = 0; l3<numberOfDimensions; l3++){
                uTotalDerivative[l3] += ((ri[l3]-rj[l3])/rLength)*uDerivative;
            }

            uTotalDoubleDerivative += uDoubleDerivative + ((numberOfDimensions-1)/rLength)*uDerivative;
        }
    }

    for (int l4 = 0; l4<numberOfDimensions; l4++){
        uAll[l4] = uTotalDerivative[l4];
    }
    uAll[numberOfDimensions] = uTotalDoubleDerivative;
                

    return uAll;
}

std::vector<double> SlaterDeterminantInteraction::computeDerivativePsi_SD(int particleIndex){
    /* This function calculates the derivative of the wavefunction with 
    regards to one spesific particle. This is used to calculate the drift
    force used in importance sampling. */
    int numberOfDimensions = m_system->getNumberOfDimensions();
    int numberOfParticles = m_system->getNumberOfParticles();

    int cor = 0;
    int pI = particleIndex;

    std::vector<double> vectorSum(numberOfDimensions, 0);

        if(particleIndex < numberOfParticles/2){
        
            for (int p3 = 0; p3<numberOfDimensions; p3++){
                for (int p1 = 0; p1 < numberOfParticles/2; p1++){
                    auto phi_00_p1 = phi_00(pI);
                    if (p1 == 0){
                        vectorSum[p3] += phi_00_der(pI)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 1){
                        vectorSum[p3] += phi_10_der(pI,0)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 2){
                        vectorSum[p3] += phi_10_der(pI,1)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 3){
                        vectorSum[p3] += phi_20_der(pI,0)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 4){
                        vectorSum[p3] += phi_11_der(pI)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }if (p1 == 5){
                        vectorSum[p3] += phi_20_der(pI,1)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinUp[pI][p1];
                    }
                }
            }

        }if (particleIndex >= numberOfParticles/2){
            pI = particleIndex - numberOfParticles/2;
            cor = numberOfParticles/2;

            for (int p3 = 0; p3<numberOfDimensions; p3++){
                for (int p1 = 0; p1 < numberOfParticles/2; p1++){
                    auto phi_00_p1 = phi_00(pI+cor);
            
                    if (p1 == 0){
                        vectorSum[p3] += phi_00_der(pI+cor)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 1){
                        vectorSum[p3] += phi_10_der(pI+cor,0)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 2){
                        vectorSum[p3] += phi_10_der(pI+cor,1)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 3){
                        vectorSum[p3] += phi_20_der(pI+cor,0)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 4){
                        vectorSum[p3] += phi_11_der(pI+cor)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }if (p1 == 5){
                        vectorSum[p3] += phi_20_der(pI+cor,1)[p3]*phi_00_p1*m_inverseSlaterMatrixSpinDown[pI][p1];
                    }
                }
            }
        }

    return vectorSum;

}

std::vector <std::vector <double> >  SlaterDeterminantInteraction::calculateInterparticleDistances(){
    /* This function calculates the distance between the particles, thus the output is a
    matrix which contains the distance between each and every particle. */

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    auto m_particles = m_system->getParticles();

    std::vector <std::vector <double> > distances(numberOfParticles, std::vector<double>(numberOfParticles, (double) 0)); // a matrix of the distance between all particles 
    
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

    m_distance = distances[0][1];

    return distances;
}

void SlaterDeterminantInteraction::setupSlaterMatrix(){
    /* This function sets up the Slater determinants. One for the spin up electrons
     (first half) and one for the spin down electrons (last half).*/
    int numberOfParticles = m_system->getNumberOfParticles();

    std::vector <std::vector <double> > slaterMatrixSpinUp(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > slaterMatrixSpinDown(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));


    for (int m1 = 0; m1 < numberOfParticles/2; m1++){
        for (int m2 = 0; m2 < numberOfParticles/2; m2++){
            
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

    m_slaterMatrixSpinDown = slaterMatrixSpinDown;
    m_slaterMatrixSpinUp = slaterMatrixSpinUp;
}

void SlaterDeterminantInteraction::updateSlaterMatrix(int particleNumber){
    /* This function updates the row in the Slater determinant 
    that represents the particle that was moved.*/

    int numberOfParticles = m_system->getNumberOfParticles();
    std::vector <std::vector <double> > m_slaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    int correction = 0;
    int pN = particleNumber;



    if(particleNumber < numberOfParticles/2){
        m_slaterMatrix = m_slaterMatrixSpinUp;
        correction = 0;
    }else{
        m_slaterMatrix = m_slaterMatrixSpinDown;
        pN = particleNumber - numberOfParticles/2;
        correction = numberOfParticles/2;
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

void SlaterDeterminantInteraction::updateSlaterRelatedThings(int particleNumber){
    updateSlaterMatrix(particleNumber);
    updateInverseSlaterMatrix();
}

void SlaterDeterminantInteraction::setupSlaterRelatedThings(){
    setupSlaterMatrix();
    calculateInverseSlaterMatrix();
}


void SlaterDeterminantInteraction::calculateInverseSlaterMatrix(){
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

        for (int i_i = 0; i_i < numberOfParticles/2; i_i++){
            for (int j_i = 0; j_i < numberOfParticles/2; j_i++){
                m(i_i,j_i) = m_slaterMatrix[i_i][j_i];
            }
        }

        // Inverse:
        m_inv = m.inverse();


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


void SlaterDeterminantInteraction::updateInverseSlaterMatrix(){
    /* This function should ideally just update the inverse Slater determinant 
    using the old inverse Slater determinant, but for now it just calculates it. */

    calculateInverseSlaterMatrix();


}

double SlaterDeterminantInteraction::computeRatio(double oldWaveFunction, double newWaveFunction){
    m_metropolisRatio = newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction);
    return m_metropolisRatio;
}


std::vector<double> SlaterDeterminantInteraction::computeQuantumForce(int particleIndex){
    /* This function calculates the quantum force/drift force with is used for importance
        sampling. The quantum force is given by the derivative of the wavefunction. */

    auto derivative = computeDerivativePsi_SD(particleIndex);

    for (int m=0;m<m_system->getNumberOfDimensions();m++){
        derivative[m] *= 2;
    }

    return derivative;
}

// -----------------------------------
// Basis:
// -----------------------------------

double SlaterDeterminantInteraction::phi_00(int particleNumber){
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

double SlaterDeterminantInteraction::phi_10(int particleNumber, int dimension){
    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    return 2*sqrt(m_omega)*r[dimension];
}

double SlaterDeterminantInteraction::phi_20(int particleNumber, int dimension){
    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    return 4*m_omega*r[dimension]*r[dimension]-2;
}

double SlaterDeterminantInteraction::phi_11(int particleNumber){
    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    return 4*m_omega*r[0]*r[1];
}

// Derivatives:

std::vector<double> SlaterDeterminantInteraction::phi_00_der(int particleNumber){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    std::vector<double> derivativeVector(numberOfDimensions);

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        derivativeVector[b0] = -m_alpha*m_omega*r[b0];
    }

    return derivativeVector;
}

std::vector<double> SlaterDeterminantInteraction::phi_10_der(int particleNumber, int dimension){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    std::vector<double> derivativeVector(numberOfDimensions);

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();
    int dim_other = 0;

    if (dimension == 0){
        dim_other = 1;
    }else{
        dim_other = 0;}


    derivativeVector[dimension] = -2*sqrt(m_omega)*(m_alpha*m_omega*r[dimension]*r[dimension]-1);
    derivativeVector[dim_other] = -2*m_alpha*sqrt(m_omega)*m_omega*r[dimension]*r[dim_other];

    return derivativeVector;
}

std::vector<double> SlaterDeterminantInteraction::phi_20_der(int particleNumber, int dimension){
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

std::vector<double> SlaterDeterminantInteraction::phi_11_der(int particleNumber){
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

double SlaterDeterminantInteraction::phi_00_double_der(int particleNumber){
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

double SlaterDeterminantInteraction::phi_10_double_der(int particleNumber, int dimension){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double rSquared = 0;

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        rSquared += r[b0]*r[b0];
    }

    return 2*m_alpha*sqrt(m_omega)*m_omega*r[dimension]*(m_alpha*m_omega*rSquared-4);
}

double SlaterDeterminantInteraction::phi_20_double_der(int particleNumber, int dimension){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double rSquared = 0;

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        rSquared += r[b0]*r[b0];
    }

    return 2*m_omega*(m_alpha*m_alpha*m_omega*(2*m_omega*r[dimension]*r[dimension]-1)*rSquared + m_alpha*(2-12*m_omega*r[dimension]*r[dimension])+4);
}

double SlaterDeterminantInteraction::phi_11_double_der(int particleNumber){
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double rSquared = 0;

    auto m_particles = m_system->getParticles();
    auto r = m_particles[particleNumber]->getPosition();

    for (int b0 = 0; b0 < numberOfDimensions; b0++){
        rSquared += r[b0]*r[b0];
    }

    return 4*m_alpha*m_omega*m_omega*r[0]*r[1]*(m_alpha*m_omega*rSquared-6);
}

