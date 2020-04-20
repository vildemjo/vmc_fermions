#pragma once
#include "wavefunction.h"

class SimpleGaussianInteraction : public WaveFunction {
public:
    SimpleGaussianInteraction(class System* system, double alpha, double SpinFactor, double beta);
    double evaluate();
    double computeDoubleDerivative();
    std::vector<double> computeDerivative(int particleIndex);
    double computeAlphaDerivative();
    double computeBetaDerivative();
    bool getDistanceCheck(std::vector <class Particle*> particles) { return calculateInterparticleDistances(particles); };

    std::vector<std::vector<double>> getDistances(std::vector <class Particle*> particles) {
        calculateInterparticleDistances(particles);
        return m_distances; }


private:
    double evaluateCorrelationPart();
    double computeInteractionPartOfDoubleDerivative();
    std::vector <double> computeDerivativeOfu( int particleNumber);
    std::vector <double> computeDerivativeOneParticle( int particleIndex);
    bool calculateInterparticleDistances(std::vector <class Particle*> particles);
    void setDistances(std::vector<std::vector<double>> distances);
    void updateDistances(int particleNumber);
    
    std::vector<std::vector<double>> m_distances;
    double m_omega = 1.0;
};
