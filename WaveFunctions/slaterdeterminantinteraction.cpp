#include "slaterdeterminantinteraction.h"
#include "interactionharmonicoscillator.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include <iostream>
#include "InitialStates/initialstate.h"


SlaterDeterminantInteraction::SlaterDeterminantInteraction(System* system, double alpha, double beta, double spinFactor) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_system->setSpinFactor(spinFactor);
    m_omega = m_system->getFrequency(); // A little bit vonrable because now the hamiltonain has to be added before the wavefunction.
}

double SlaterDeterminantInteraction::evaluate() {
    /* This function calculated the trial wavefuntion in the case 
    with interaction. */

    // auto rSum = calculatePositionSumSquared();
    // double norm = 1;

    double interactionPart = evaluateCorrelationPart();

    return m_slaterDeterminantSpinUp*m_slaterDeterminantSpinDown*interactionPart;

}

double SlaterDeterminantInteraction::evaluateCorrelationPart() {
    /* This function calculates the correlation/interaction part of the 
    trial wavefunction. */

    double correlationPart = 0;
  
    auto a = m_system->getSpinFactor();

    auto distances = calculateInterparticleDistances();

    int numberOfParticles = m_system->getNumberOfParticles();

    for (int j7 = 0; j7<numberOfParticles-1; j7++){
        for (int j8 = j7+1; j8<numberOfParticles; j8++){
            correlationPart += a*distances[j7][j8]/(1+m_parameters[1]*distances[j7][j8]);
        }
    }
    
    return exp(correlationPart);

}

double SlaterDeterminantInteraction::computeDoubleDerivative() {

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double      doubleDerivative    = 0;

    double interactionPart = 0;

    auto rSum2 = calculatePositionSumSquared();

    interactionPart = computeInteractionPartOfDoubleDerivative();


    for(int p0 = 0; p0 < numberOfParticles/2; p0++){
        int p0_ = p0 + numberOfParticles/2;
        auto phi_00_p0 = phi_00(p0);
        auto phi_00_p0_ = phi_00(p0_);
        for (int p1 = 0; p1 < numberOfParticles/2; p1++){
            if (p1 == 0){
                doubleDerivative += phi_00_double_der(p0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p1][p0];
                doubleDerivative += phi_00_double_der(p0_)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p1][p0];
            }if (p1 == 1){
                doubleDerivative += phi_10_double_der(p0,0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p1][p0];
                doubleDerivative += phi_10_double_der(p0_,0)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p1][p0];
            }if (p1 == 2){
                doubleDerivative += phi_10_double_der(p0,1)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p1][p0];
                doubleDerivative += phi_10_double_der(p0_,1)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p1][p0];
            }if (p1 == 3){
                doubleDerivative += phi_20_double_der(p0,0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p1][p0];
                doubleDerivative += phi_20_double_der(p0_,0)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p1][p0];
            }if (p1 == 4){
                doubleDerivative += phi_11_double_der(p0)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p1][p0];
                doubleDerivative += phi_11_double_der(p0_)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p1][p0];
            }if (p1 == 5){
                doubleDerivative += phi_20_double_der(p0, 1)*phi_00_p0*m_inverseSlaterMatrixSpinUp[p1][p0];
                doubleDerivative += phi_20_double_der(p0_,1)*phi_00_p0_*m_inverseSlaterMatrixSpinDown[p1][p0];
            }
        }
    }
    return doubleDerivative + interactionPart;
}

std::vector<double> SlaterDeterminantInteraction::computeDerivative(int particleIndex){
    /* This function calculates the derivative of the trial wavefunction with 
    interaction with regards to one particle. */


    int                     numberOfDimensions          = m_system->getNumberOfDimensions();
    std::vector <double>    vectorWithInteraction       (numberOfDimensions);

    auto m_particles         = m_system->getParticles();
    auto ri                  = m_particles[particleIndex]->getPosition();


    // The derivative of the interaction part is calculated in
    // another function, because the same expression is used to 
    // evaluate the double derivative.
    auto derivative_psi_in   = computeDerivativeOfu(particleIndex);

    int numberOfParticles = m_system->getNumberOfParticles();

    std::vector <std::vector <double> > m_inverseSlaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    int cor = 0;
    int pI = 0;

    std::vector<double> derivative_psi_ob(numberOfDimensions);

    auto r = m_particles[particleIndex]->getPosition();


    // if (numberOfParticles == 2){
    //     for (int n8=0; n8<numberOfDimensions; n8++){
    //         vectorSum[n8] = -m_omega*getParameters()[0]*r[n8];
    //     }

    //     return vectorSum;
    // }else{

        double R_inv = 1/m_metropolisRatio;

        if(particleIndex < numberOfParticles/2){
            m_inverseSlaterMatrix = m_slaterMatrixSpinUp;
            cor = 0;
        }else{
            m_inverseSlaterMatrix = m_slaterMatrixSpinDown;
            pI = particleIndex - numberOfParticles/2;
            cor = numberOfParticles/2;
        }

        for (int p1 = 0; p1 < numberOfParticles/2; p1++){
            for (int p3 = 0; p3<numberOfDimensions; p3++){
                if (p1 == 0){
                    derivative_psi_ob[p3] += phi_00_der(pI+cor)[p3];
                }if (p1 == 1){
                    derivative_psi_ob[p3] += phi_10_der(pI+cor,0)[p3];
                }if (p1 == 2){
                    derivative_psi_ob[p3] += phi_10_der(pI+cor,1)[p3];
                }if (p1 == 3){
                    derivative_psi_ob[p3] += phi_20_der(pI+cor,0)[p3];
                }if (p1 == 4){
                    derivative_psi_ob[p3] += phi_11_der(pI+cor)[p3];
                }if (p1 == 5){
                    derivative_psi_ob[p3] += phi_20_der(pI+cor,1)[p3];
                }
            }
        }


    for (int n8=0; n8<numberOfDimensions; n8++){
        vectorWithInteraction[n8] = derivative_psi_ob[n8] + derivative_psi_in[n8];
    }

    return vectorWithInteraction;
}

double SlaterDeterminantInteraction::computeAlphaDerivative(){
    /* This function calculates the normalized derivative of the wavefunction 
    with regards to the parameter alpha. This is used to perform optimization
    by the use of gradient descent methods.*/

    auto m_particles = m_system->getParticles();

    auto vectorSumSquared =  calculatePositionSumSquared();

    // std::cout << "vector sum sq.: " << vectorSumSquared << std::endl;

    return (-0.5)*m_omega*vectorSumSquared;

}

double SlaterDeterminantInteraction::computeBetaDerivative(){
    /* This function calculates the normalized derivative of the wavefunction 
    with regards to the parameter beta. This is used to perform optimization
    by the use of gradient descent methods.*/

    auto m_particles = m_system->getParticles();
    double factorSum = 0;
    int numberOfParticles = m_system->getNumberOfParticles();
    auto distances = calculateInterparticleDistances();
    auto a = m_system->getSpinFactor();


    for (int j9 = 0; j9<numberOfParticles-1; j9++){
        for (int j10 = j9+1; j10<numberOfParticles; j10++){
            auto rLength = distances[j9][j10];
            double factor = (1+m_parameters[1]*rLength);
            factorSum += -a*rLength*rLength/(factor*factor);
        }
    }

    // std::cout << "factor sum: " << factorSum << std::endl;

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
        auto derivativePhi = computeDerivativePhiOB(l5);

        // uPart[0], uPart[1] and uPart[2] is the derivative of
        // the interaction part of the trial wavefunction which is
        // a vector. uPart[3] is the last term which is a part of the
        // double derivative of the interaction part of the trial
        // wavefunction. 
        // See computeDerivativeOfu for more details.

         uPart = computeDerivativeOfu(l5);


        for (int l4 = 0; l4<numberOfDimensions; l4++){
            firstTerm += derivativePhi[l4]*uPart[l4];
            secondTerm += uPart[l4]*uPart[l4];
        }
        thirdTerm += uPart[numberOfDimensions];
        
    }

    return 2*firstTerm + secondTerm + thirdTerm;
}

std::vector <double> SlaterDeterminantInteraction::computeDerivativeOfu(int particleNumber){
    /* This function calcualtes the derivative and double derivative of the interation term 
    in the trial wavefunction which are used in calcualting the double derivative and the 
    derivative with regards to one particle. */

    int                     numberOfParticles       = m_system->getNumberOfParticles();
    int                     numberOfDimensions      = m_system->getNumberOfDimensions();
    double                  uDerivative             = 0;
    double                  a                       = m_system->getSpinFactor();
    double                  uDoubleDerivative       = 0;
    double                  uTotalDoubleDerivative  = 0;
    std::vector <double>    uTotalDerivative        (numberOfDimensions);
    std::vector <double>    uAll                    (numberOfDimensions+1);


    auto m_particles = m_system->getParticles();
    auto ri          = m_particles[particleNumber]->getPosition();   

    auto difference = calculateInterparticleDistances();


    

    for (int l1 = 0; l1 < numberOfParticles; l1++){
        if (particleNumber != l1){
            auto rj         = m_particles[l1]->getPosition();        
            auto rLength    = difference[particleNumber][l1];                           

            // Here sum u'(r_ij) is determined based on the relationship 
            // between r_ij (distance between particles) and a (hard core diameter)

            auto factor = (1+m_parameters[1]*rLength);
            uDerivative = a/(factor*factor);
            uDoubleDerivative = -2*a*m_parameters[1]/(factor*factor*factor);
       

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

std::vector<double> SlaterDeterminantInteraction::computeDerivativePhiOB(int particleIndex){
    /* This function calculates the derivative of the wavefunction with 
    regards to one spesific particle. This is used to calcualte the drift
    force used in importance sampling. */
    
    int numberOfDimensions = m_system->getNumberOfDimensions();
    int numberOfParticles = m_system->getNumberOfParticles();
    auto m_particles = m_system->getParticles();

    std::vector <std::vector <double> > m_inverseSlaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    int cor = 0;
    int pI = 0;

    std::vector<double> vectorSum(numberOfDimensions);

    auto r = m_particles[particleIndex]->getPosition();

        double R_inv = 1/m_metropolisRatio;

        if(particleIndex < numberOfParticles/2){
            m_inverseSlaterMatrix = m_slaterMatrixSpinUp;
            cor = 0;
        }else{
            m_inverseSlaterMatrix = m_slaterMatrixSpinDown;
            pI = particleIndex - numberOfParticles/2;
            cor = numberOfParticles/2;
        }

        for (int p1 = 0; p1 < numberOfParticles/2; p1++){
            for (int p3 = 0; p3<numberOfDimensions; p3++){
                if (p1 == 0){
                    vectorSum[p3] += phi_00_der(pI+cor)[p3];
                }if (p1 == 1){
                    vectorSum[p3] += phi_10_der(pI+cor,0)[p3];
                }if (p1 == 2){
                    vectorSum[p3] += phi_10_der(pI+cor,1)[p3];
                }if (p1 == 3){
                    vectorSum[p3] += phi_20_der(pI+cor,0)[p3];
                }if (p1 == 4){
                    vectorSum[p3] += phi_11_der(pI+cor)[p3];
                }if (p1 == 5){
                    vectorSum[p3] += phi_20_der(pI+cor,1)[p3];
                }
            }
        }
        return vectorSum;

}

std::vector <std::vector <double> >  SlaterDeterminantInteraction::calculateInterparticleDistances(){
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
    int numberOfParticles = m_system->getNumberOfParticles();

    std::vector <std::vector <double> > slaterMatrixSpinUp(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > slaterMatrixSpinDown(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    // Test if number of particles are 2, 6 or 12?


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
    int numberOfParticles = m_system->getNumberOfParticles();

    std::vector <std::vector <double> > slaterMatrixSpinUp(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > slaterMatrixSpinDown(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    // Test if number of particles are 2, 6 or 12?

    // std::cout << "calc wf: " << phi_00(0) <<"\n";

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

    // std::cout << "ok here \n";

    m_slaterMatrixSpinDown = slaterMatrixSpinDown;
    m_slaterMatrixSpinUp = slaterMatrixSpinUp;
}

void SlaterDeterminantInteraction::updateSlaterRelatedThings(int particleNumber){
    updateSlaterMatrix(particleNumber);
    calculateInverseSlaterMatrix();
}

void SlaterDeterminantInteraction::setupSlaterRelatedThings(){
    setupSlaterMatrix();
    calculateInverseSlaterMatrix();
}



void SlaterDeterminantInteraction::calculateInverseSlaterMatrix(){
    // Need to find out how to calculate this.
    // Should I do this in the initiation of this class?
    // Should I have another class for the slater determinant parts?
    int numberOfParticles = m_system->getNumberOfParticles();
    int i_, j_;

    std::vector <std::vector <double> > mat(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > m_slaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
    std::vector <std::vector <double> > m_inverseSlaterMatrix(numberOfParticles/2, std::vector<double>(numberOfParticles/2, (double) 0));
	float determinant = 0;
	
  

    for (int s = 0; s < 2; s++){

        if (s == 0){
            m_slaterMatrix = m_slaterMatrixSpinUp;
            m_oldInverseSlaterMatrixSpinUp = m_inverseSlaterMatrixSpinUp; // saving the old one
        }if(s == 1){
            m_slaterMatrix = m_slaterMatrixSpinDown;
            m_oldInverseSlaterMatrixSpinDown = m_inverseSlaterMatrixSpinDown; // saving the old one
        }


        for(i_ = 0; i_ < numberOfParticles/2; i_++){
            // std::cout << "\n";
            for(j_ = 0; j_ < numberOfParticles/2; j_++){
                mat[i_][j_] = m_slaterMatrix[i_][j_];
                // std::cout << mat[i_][j_] << "\t";
            }
        }

        // finding determinant
        for(i_ = 0; i_ < numberOfParticles/2; i_++)
            if (numberOfParticles == 6){
                determinant += (mat[0][i_] * (mat[1][(i_+1)%3] * mat[2][(i_+2)%3] - mat[1][(i_+2)%3] * mat[2][(i_+1)%3]));
            }if (numberOfParticles == 2){
                determinant = mat[0][0];
            }
        // std::cout << "\n ------- \n determinant: " << determinant << "\n ---------- \n";

        for(i_ = 0; i_ < numberOfParticles/2; i_++){
            // std::cout << "\n";
            for(j_ = 0; j_ < numberOfParticles/2; j_++){
                if (numberOfParticles == 6){
                    m_inverseSlaterMatrix[i_][j_] = ((mat[(j_+1)%3][(i_+1)%3] * mat[(j_+2)%3][(i_+2)%3]) - (mat[(j_+1)%3][(i_+2)%3] * mat[(j_+2)%3][(i_+1)%3]))/ determinant;
                // std::cout << m_inverseSlaterMatrix[i_][j_] << "\t";
                }if (numberOfParticles == 2){
                    m_inverseSlaterMatrix[i_][j_] = 1/mat[i_][j_];
                }
            }
        }
    
        if (s == 0){
            m_slaterMatrixSpinUp = m_slaterMatrix;
            m_inverseSlaterMatrixSpinUp = m_inverseSlaterMatrix;
            m_slaterDeterminantSpinUp = determinant;
        }else{
            m_slaterMatrixSpinDown = m_slaterMatrix;
            m_inverseSlaterMatrixSpinDown = m_inverseSlaterMatrix;
            m_slaterDeterminantSpinDown = determinant;
        }
    

    }


}



// void SlaterDeterminantInteraction::updateInverseSlaterMatrix(int particleNumber){
//     int pN = particleNumber;
//     int numberOfParticles = m_system->getNumberOfParticles();
//     std::vector<double> S(numberOfParticles, 0);



//     for (int row = 0; row < numberOfParticles; row++){
//         for (int col = 0; col < numberOfParticles; col++){
//                 for (int l = 0; l <numberOfParticles; l++){
//                     S[col] += m_oldSlaterMatrix[pN][l];
//                 }
//             if (row != col){
//                 m_inverseSlaterMatrix[col][row] = m_oldInverseSlaterMatrix[col][row] - S[col]/R*m_oldInverseSlaterDeterminant[pN][row];
//             } // Jeg trenger både den gamle Slater-determinanten og den nye... Skal jeg lagre to stykker?
//         }
//     }

// }

double SlaterDeterminantInteraction::computeRatio(double oldWaveFunction, double newWaveFunction){
    return newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction);
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

