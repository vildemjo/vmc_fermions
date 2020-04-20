#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);

    double              computeKineticEnergy();
    double              computePotentialEnergy();
    double              computeLocalEnergy();
    std::vector<double> computeQuantumForce(int particleIndex);

private:
    double m_omega = 0;
};

