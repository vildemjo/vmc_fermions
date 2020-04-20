#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double stepLength);
    void setupInitialState();

private:
    double m_stepLength = 0.5;

};

