#pragma once
#include "wavefunction.h"

class SlaterDeterminantInteraction : public WaveFunction {
public:
    SlaterDeterminantInteraction(class System* system, double alpha, double beta, double SpinFactor);
    double evaluate();
    double computeDoubleDerivative();
    std::vector<double> computeDerivative(int particleIndex);
    double computeAlphaDerivative();
    double computeBetaDerivative();
    double getDistance()  { return m_distance; }
    double computeRatio(double oldWaveFunction, double newWaveFunction);
    void updateSlaterRelatedThings(int particleNumber);
    void setupSlaterRelatedThings();


private:
    double evaluateCorrelationPart();
    double computeInteractionPartOfDoubleDerivative();
    std::vector <std::vector <double> >  calculateInterparticleDistances();
    std::vector <double> computeDerivativeOfu( int particleNumber);
    std::vector <double> computeDerivativeOneParticle( int particleIndex);
    double m_omega = 1.0;
    double m_distance = 0;

    double m_slaterDeterminantSpinUp = 0;
    double m_slaterDeterminantSpinDown = 0;

    std::vector<std::vector<double>> m_slaterMatrixSpinUp = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_slaterMatrixSpinDown = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_oldSlaterMatrixSpinUp = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_oldSlaterMatrixSpinDown = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_inverseSlaterMatrixSpinUp = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_inverseSlaterMatrixSpinDown = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_oldInverseSlaterMatrixSpinUp = std::vector<std::vector<double>>();
    std::vector<std::vector<double>> m_oldInverseSlaterMatrixSpinDown = std::vector<std::vector<double>>();

    void setupSlaterMatrix();
    void updateSlaterMatrix(int particleNumber);
    double phi_00(int particleNumber);
    double phi_10(int particleNumber, int dimension);
    double phi_20(int particleNumber, int dimension);
    double phi_11(int particleNumber);

    void calculateInverseSlaterMatrix();
    void updateInverseSlaterMatrix(int particleNumber);
};
