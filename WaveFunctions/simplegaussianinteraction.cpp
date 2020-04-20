#include "simplegaussianinteraction.h"
#include "ellipticalharmonicoscillator.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include "InitialStates/initialstate.h"
#include <iostream>

SimpleGaussianInteraction::SimpleGaussianInteraction(System* system, double alpha, double hardCoreDiameter, double beta) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    m_system->setHardCoreDiameter(hardCoreDiameter);
    m_beta = beta;
}

double SimpleGaussianInteraction::evaluate() {
    /* This function calculated the trial wavefuntion in the case 
    with interaction. */

    auto rSum = calculatePositionSumSquared();

    double interactionPart = evaluateCorrelationPart();

    return exp(-m_parameters[0]*rSum)*interactionPart;

}

double SimpleGaussianInteraction::evaluateCorrelationPart() {
    /* This function calculates the correlation/interaction part of the 
    trial wavefunction. */

    int     numberOfParticles       = m_system->getNumberOfParticles();
    int     uSumCheck               = 0;
    double  correlationPart         = 1;
    double  uSum                    = 0;
    double  distances_j1_j2         = 0;
    auto    a                       = m_system->getHardCoreDiameter();
    // The interparticle distances is calculated when the distances is needed
    auto    distances               = getDistances(m_system->getParticles());
    
    // Here the function u is evaluated
    for (int j1 = 0; j1 < numberOfParticles-1; j1++){

        for (int j2 = j1+1; j2 <numberOfParticles; j2++){
            distances_j1_j2 = distances[j1][j2];
            
            if ( distances_j1_j2 <= a ) {
                uSumCheck += 1;
            }else{
                uSum += log(1-a/distances_j1_j2);
            }
        }
    }
    
    // Here the interaction part is set to zero directly (instead of exp(-infty))
    // if one of the distances are less than a
    if (uSumCheck > 0){
        correlationPart = 0;
    }else{
        correlationPart = exp(uSum);
    }

    return correlationPart;

}

double SimpleGaussianInteraction::computeDoubleDerivative() {

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double interactionPart = 0;

    auto rSum2 = calculatePositionSumSquared();

    // if ( m_system->getSteps() == 4533+(int)1e5){
    //     std::cout << "step: " << m_system->getSteps()-1e5 << std::endl;
    //     std::cout << "rSum2: " << rSum2 << std::endl;
    // }

    interactionPart = computeInteractionPartOfDoubleDerivative();

    // if (m_system->getSteps() == 4533+(int)1e5){
    //     std::cout << "step: " << m_system->getSteps()-1e5 << std::endl;
    //     std::cout << "Interaction part: " << interactionPart << std::endl;
    // }
    
    return (-2*m_parameters[0]*numberOfParticles*numberOfDimensions + 4*m_parameters[0]*m_parameters[0]*rSum2) + interactionPart;
}

std::vector<double> SimpleGaussianInteraction::computeDerivative(int particleIndex){
    /* This function calculates the derivative of the trial wavefunction with 
    interaction with regards to one particle. */


    int                     numberOfDimensions          = m_system->getNumberOfDimensions();
    double                  derivative_psi_ob           = 0;
    std::vector <double>    vectorWithInteraction       (numberOfDimensions);

    auto m_particles         = m_system->getParticles();
    auto ri                  = m_particles[particleIndex]->getPosition();

    // The derivative of the interaction part is calculated in
    // another function, because the same expression is used to 
    // evaluate the double derivative.
    auto derivative_psi_in   = computeDerivativeOfu(particleIndex);

    for (int n8=0; n8<numberOfDimensions; n8++){
        if (n8 == 2){
            derivative_psi_ob = -2*getParameters()[0]*ri[n8]*m_beta;
        }else{
            derivative_psi_ob = -2*getParameters()[0]*ri[n8];
        }
        vectorWithInteraction[n8] = derivative_psi_ob + derivative_psi_in[n8];
    }

    return vectorWithInteraction;
}

double SimpleGaussianInteraction::computeAlphaDerivative(){
    /* This function calculates the normalized derivative of the wavefunction 
    with regards to the parameter alpha. This is used to perform optimization
    by the use of gradient descent methods.*/

    auto m_particles = m_system->getParticles();

    auto vectorSumSquared =  calculatePositionSumSquared();

    return (-1)*vectorSumSquared; // No interaction part because it is divided away.

}

double SimpleGaussianInteraction::computeInteractionPartOfDoubleDerivative(){
    /* This function calculates the interaction part of the double derivative
    when it is evaluated analytically. */

    // These terms represents the three terms the interaction is involved in
    // in the analytical expression of the double derivative
    double      firstTerm = 0; double secondTerm = 0; double thirdTerm = 0;

    int         numberOfParticles = m_system->getNumberOfParticles();
    int         numberOfDimensions = m_system->getNumberOfDimensions();
    auto        m_particles = m_system->getParticles();


    for(int l5 = 0; l5 < numberOfParticles; l5++){ 

        // The first term also includes the non-interating part of
        // the trial wavefunction. This is evaluated here.
        auto derivativePhi = computeDerivativeOneParticle(l5);

        // uPart[0], uPart[1] and uPart[2] is the derivative of
        // the interaction part of the trial wavefunction which is
        // a vector. uPart[3] is the last term which is a part of the
        // double derivative of the interaction part of the trial
        // wavefunction. 
        // See computeDerivativeOfu for more details.
        auto uPart = computeDerivativeOfu(l5);

        for (int l4 = 0; l4<numberOfDimensions; l4++){
            firstTerm += derivativePhi[l4]*uPart[l4];
            secondTerm += uPart[l4]*uPart[l4];
        }
        thirdTerm += uPart[numberOfDimensions+1];
        
    }

    return 2*firstTerm + secondTerm + thirdTerm;
}

std::vector <double> SimpleGaussianInteraction::computeDerivativeOfu(int particleNumber){
    /* This function calcualtes the derivative and double derivative of the interation term 
    in the trial wavefunction which are used in calcualting the double derivative and the 
    derivative with regards to one particle. */

    int                     numberOfParticles       = m_system->getNumberOfParticles();
    int                     numberOfDimentions      = m_system->getNumberOfDimensions();
    double                  uDerivative             = 0;
    double                  a                       = m_system->getHardCoreDiameter();
    double                  uDoubleDerivative       = 0;
    double                  uTotalDoubleDerivative  = 0;
    std::vector <double>    uTotalDerivative        (numberOfDimentions);
    std::vector <double>    uAllStuff               (numberOfDimentions+1);

    auto m_particles = m_system->getParticles();
    auto ri         = m_particles[particleNumber]->getPosition();        
    auto difference = getDistances(m_particles);
    

    for (int l1 = 0; l1 < numberOfParticles; l1++){
        if (particleNumber != l1){
            auto rj         = m_particles[l1]->getPosition();        
            auto rLength    = difference[particleNumber][l1];                           

            // Here sum u'(r_ij) is determined based on the relationship 
            // between r_ij (distance between particles) and a (hard core diameter)

            if (rLength <= a){ // this case should not happen because the wavefuntion is 0 if
                               // an interparticle distance is smaller than a. And then the energy
                               // is not sampled.

                uDerivative = -1e50;
                uDoubleDerivative = -1e50;
                std::cout << "r_ij < a in computeDerivativeOfu" << std::endl;
            }else{
                uDerivative = -a/(a*rLength-rLength*rLength);
                uDoubleDerivative = a*(a-2*rLength)/(rLength*rLength*(a-rLength)*(a-rLength));
            }

            for (int l3 = 0; l3<numberOfDimentions; l3++){
                uTotalDerivative[l3] += ((ri[l3]-rj[l3])/rLength)*uDerivative;
            }

            uTotalDoubleDerivative += uDoubleDerivative + ((numberOfDimentions-1)/rLength)*uDerivative;
        }
    }

    for (int l4 = 0; l4<numberOfDimentions; l4++){
        uAllStuff[l4] = uTotalDerivative[l4];
    }
    uAllStuff[numberOfDimentions+1] = uTotalDoubleDerivative;

    return uAllStuff;
}

std::vector<double> SimpleGaussianInteraction::computeDerivativeOneParticle(int particleIndex){
    /* This function calcualtes the derivative of the not interacting term of the 
    trial wavefunction with regards to one particle. This is used to calculate the double
    derivative. */

    int                 numberOfDimensions          = m_system->getNumberOfDimensions();
    std::vector<double> derivativeVector            (numberOfDimensions);
    auto m_particles = m_system->getParticles();

    auto r = m_particles[particleIndex]->getPosition();

    for (int j3=0; j3<numberOfDimensions; j3++){
        if (j3 == 2){
            derivativeVector[j3] = -2*getParameters()[0]*m_beta*r[j3];
        }else{
            derivativeVector[j3] = -2*getParameters()[0]*r[j3];
        }
    }

    return derivativeVector;
    
}

bool SimpleGaussianInteraction::calculateInterparticleDistances(std::vector <class Particle*> particles){
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector <std::vector <double> > distances(numberOfParticles, std::vector<double>(numberOfParticles, (double) 0)); // a matrix of the distance between all particles 
    
    double a = m_system->getHardCoreDiameter();
    
    std::vector <class Particle*> m_particles = particles;

    int distanceCheck = 0;
    std::vector <double> r1(numberOfDimensions), r2(numberOfDimensions);

    int times = 0;

    for (int k1 = 0; k1 < numberOfParticles; k1++){

        r1 = m_particles[k1]->getPosition();
        
        for (int k2 = 0; k2 <numberOfParticles; k2++){
            
             if (k1 != k2){
                r2 = m_particles[k2]->getPosition();
                for (int k3 = 0; k3<numberOfDimensions; k3++){
                    distances[k1][k2] +=  (r1[k3]-r2[k3])*(r1[k3]-r2[k3]);    
                }
                distances[k1][k2] = sqrt(distances[k1][k2]);      
               
                if (distances[k1][k2] < a){
                     distanceCheck +=1;
                }
            }else{
                distances[k1][k2] = 0;
            }
        }
    }

    setDistances(distances);

    // Here the function returns false if any of the particles are closer
    // to each other than a. This is used in the initializing of the particle
    // positions, so that we start with particles that are far enough apart.
    if (distanceCheck > 0){
        return false;
    }
 
    return true;
}

void SimpleGaussianInteraction::setDistances(std::vector<std::vector<double>> distances){
    m_distances = distances;
}

double SimpleGaussianInteraction::calculatePositionSumSquared(){
    /* This function calcualtes the positions squared and summed over all
    particles. This is used in several of the functions in this class. */

    double vectorSumSquared = 0.0;
    auto m_particles = m_system->getParticles();
    
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    for (int i10 = 0; i10<numberOfParticles;i10++){
        auto r = m_particles[i10]->getPosition();
        for (int n10=0; n10<numberOfDimensions; n10++){
            if (n10 == 2){
                vectorSumSquared += m_beta*r[n10]*r[n10];
            }else{
                vectorSumSquared += r[n10]*r[n10];
            }
        }
    }


    return vectorSumSquared;
}