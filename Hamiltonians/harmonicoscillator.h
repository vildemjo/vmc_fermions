#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);

    double              computeKineticEnergy();
    double              computePotentialEnergy();
    double              computeInteractionEnergy(){return 0;};
    double              computeLocalEnergy();

private:
    double m_omega = 0;
};

