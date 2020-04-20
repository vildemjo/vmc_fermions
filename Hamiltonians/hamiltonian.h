#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);

    double                      computeDoubleDerivativeNumerically(std::vector<class Particle*> particles);
    virtual double              computeKineticEnergy(std::vector<Particle*> particles) = 0;
    virtual double              computePotentialEnergy(std::vector<Particle*> particles) = 0;
    virtual double              computeLocalEnergy(std::vector<Particle*> particles) = 0;
    virtual std::vector<double> computeQuantumForce(int particleIndex, std::vector<class Particle*> particles) = 0;


protected:
    class System* m_system = nullptr;
};

