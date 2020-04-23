#pragma once
#include "hamiltonian.h"
#include <vector>

class InteractionHarmonicOscillator : public Hamiltonian {
public:
    InteractionHarmonicOscillator(System* system, double omega);

    double                  computeKineticEnergy();
    double                  computePotentialEnergy();
    double                  computeLocalEnergy();
    std::vector<double>     computeQuantumForce         (int particleIndex);
    
private:
    double                  m_omega = 0;
    double              computeInteractionEnergy();
    std::vector <std::vector <double> >  calculateInterparticleDistances();
};

