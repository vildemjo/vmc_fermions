#pragma once
#include "initialstate.h"

class GaussianDistribution : public InitialState {
public:
    GaussianDistribution(System* system, int numberOfDimensions, int numberOfParticles, double stepLength);
    void setupInitialState();

private:
    double m_stepLength = 0.5;

};

