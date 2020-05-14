#pragma once
#include "hamiltonian.h"
#include <vector>

class InteractionHarmonicOscillator : public Hamiltonian {
public:
    InteractionHarmonicOscillator(System* system, double omega);

    double                  computeKineticEnergy();
    double                  computePotentialEnergy();
    double                  computeInteractionEnergy();
    double                  computeLocalEnergy();

private:
    double                  m_omega = 0;

    std::vector <std::vector <double> >  calculateInterparticleDistances();

};

