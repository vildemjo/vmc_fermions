#include "simplegaussianinteraction.h"
#include "interactionharmonicoscillator.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include "InitialStates/initialstate.h"
#include <iostream>

SimpleGaussianInteraction::SimpleGaussianInteraction(System* system, double alpha, double beta, double spinFactor):
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_system->setSpinFactor(spinFactor);
    m_omega = m_system->getFrequency(); // A little bit vounrable because now the hamiltonain has to be added before the wavefunction.
}

double SimpleGaussianInteraction::evaluate() {
    /* This function calculated the trial wavefuntion in the case 
    with interaction. */

    auto rSum = calculatePositionSumSquared();
    double norm = 1;

    double interactionPart = evaluateCorrelationPart();


    return norm*exp(-0.5*m_parameters[0]*m_omega*rSum)*interactionPart;

}

double SimpleGaussianInteraction::evaluateCorrelationPart() {
    /* This function calculates the correlation/interaction part of the 
    trial wavefunction. */

    double correlationPart = 0;
  
    auto a = m_system->getSpinFactor();

    auto distances = calculateInterparticleDistances(m_system->getParticles());

    int numberOfParticles = m_system->getNumberOfParticles();

    for (int j7 = 0; j7<numberOfParticles-1; j7++){
        for (int j8 = j7+1; j8<numberOfParticles; j8++){
            correlationPart *= exp(a*distances[j7][j8]/(1+m_parameters[1]*distances[j7][j8]));
        }
    }
    
    return correlationPart;

}

double SimpleGaussianInteraction::computeDoubleDerivative() {

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double interactionPart = 0;

    auto rSum2 = calculatePositionSumSquared();

    interactionPart = computeInteractionPartOfDoubleDerivative();
    
    return (-m_omega*m_parameters[0]*numberOfParticles*numberOfDimensions 
                        + m_omega*m_omega*m_parameters[0]*m_parameters[0]*rSum2) + interactionPart;
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
        derivative_psi_ob = -getParameters()[0]*m_omega*ri[n8];
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

    return (-0.5)*m_omega*vectorSumSquared; // No interaction part because it is divided away.

}

double SimpleGaussianInteraction::computeBetaDerivative(){
    /* This function calculates the normalized derivative of the wavefunction 
    with regards to the parameter alpha. This is used to perform optimization
    by the use of gradient descent methods.*/

    auto m_particles = m_system->getParticles();
    double factorSum = 0;
    int numberOfParticles = m_system->getNumberOfParticles();
    auto distances = calculateInterparticleDistances(m_particles);
    auto a = m_system->getSpinFactor();


    for (int j9 = 0; j9<numberOfParticles-1; j9++){
        for (int j10 = j10+1; j10<numberOfParticles; j10++){
            auto rLength = distances[j9][j10];
            double factor = (1+m_parameters[1]*rLength);
            factorSum += -a*rLength*rLength/(factor*factor);
        }
    }


    return factorSum; 

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
    double                  a                       = m_system->getSpinFactor();
    double                  uDoubleDerivative       = 0;
    double                  uTotalDoubleDerivative  = 0;
    std::vector <double>    uTotalDerivative        (numberOfDimentions);
    std::vector <double>    uAllStuff               (numberOfDimentions+1);


    auto m_particles = m_system->getParticles();
    auto ri         = m_particles[particleNumber]->getPosition();   

    auto difference = calculateInterparticleDistances(m_particles);

    

    for (int l1 = 0; l1 < numberOfParticles; l1++){
        if (particleNumber != l1){
            auto rj         = m_particles[l1]->getPosition();        
            auto rLength    = difference[particleNumber][l1];                           

            // Here sum u'(r_ij) is determined based on the relationship 
            // between r_ij (distance between particles) and a (hard core diameter)

            auto factor = (1+m_parameters[1]*rLength);
            uDerivative = a/(factor*factor);
            uDoubleDerivative = -2*a*m_parameters[1]/(factor*factor*factor);
       

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
        derivativeVector[j3] = -getParameters()[0]*m_omega*r[j3];
    }

    return derivativeVector;
    
}

std::vector <std::vector <double> >  SimpleGaussianInteraction::calculateInterparticleDistances(std::vector <class Particle*> particles){
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector <std::vector <double> > distances(numberOfParticles, std::vector<double>(numberOfParticles, (double) 0)); // a matrix of the distance between all particles 
    
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
            }else{
                distances[k1][k2] = 0;
            }
        }
    }

    return distances;
}