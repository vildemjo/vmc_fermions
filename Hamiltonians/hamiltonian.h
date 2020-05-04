#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);

    double                      computeDoubleDerivativeNumerically();
    virtual double              computeKineticEnergy() = 0;
    virtual double              computePotentialEnergy() = 0;
    virtual double              computeInteractionEnergy() = 0;
    virtual double              computeLocalEnergy() = 0;
    virtual std::vector<double> computeQuantumForce(int particleIndex) = 0;


protected:
    class System* m_system = nullptr;
};

