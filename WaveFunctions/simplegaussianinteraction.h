#pragma once
#include "wavefunction.h"

class SimpleGaussianInteraction : public WaveFunction {
public:
    SimpleGaussianInteraction(class System* system, double alpha, double beta, double SpinFactor);
    double evaluate();
    double computeDoubleDerivative();
    std::vector<double> computeDerivative(int particleIndex);
    double computeAlphaDerivative();
    double computeBetaDerivative();
    double getDistance()  { return m_distance; }


private:
    double evaluateCorrelationPart();
    double computeInteractionPartOfDoubleDerivative();
    std::vector <std::vector <double> >  calculateInterparticleDistances();
    std::vector <double> computeDerivativeOfu( int particleNumber);
    std::vector <double> computeDerivativeOneParticle( int particleIndex);
    double m_omega = 1.0;
    double m_distance = 0;

    std::vector<std::vector<double>> m_slaterDeterminantSpinUp = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_slaterDeterminantSpinDown = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_oldSlaterDeterminantSpinUp = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_oldSlaterDeterminantSpinDown = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_inverseSlaterDeterminantSpinUp = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_inverseSlaterDeterminantSpinDown = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_oldInverseSlaterDeterminantSpinUp = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_oldInverseSlaterDeterminantSpinDown = std::vector<std::vector<double>>();

    void setupSlaterDeterminant();
    void updateSlaterDeterminant(int particleNumber);
    double phi_00(int particleNumber);
    double phi_10(int particleNumber, int dimension);
    double phi_20(int particleNumber, int dimension);
    double phi_11(int particleNumber);

    void calculateInverseSlaterDeterminant();
    void updateInverseSlaterDeterminant(int particleNumber);
};
