#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate();
    double computeDoubleDerivative();
    std::vector<double> computeDerivative(int particleIndex);
    double computeAlphaDerivative();
    double computeBetaDerivative(){ return 0;}
    double getDistance();
};
