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


private:
    double evaluateCorrelationPart();
    double computeInteractionPartOfDoubleDerivative();
    std::vector <std::vector <double> >  calculateInterparticleDistances();
    std::vector <double> computeDerivativeOfu( int particleNumber);
    std::vector <double> computeDerivativeOneParticle( int particleIndex);
    double m_omega = 1.0;
};
