#pragma once
#include "wavefunction.h"

class SlaterDeterminant : public WaveFunction {
public:
    SlaterDeterminant(class System* system, double alpha);
    double evaluate();
    double computeDoubleDerivative();
    std::vector<double> computeDerivative(int particleIndex);
    double computeAlphaDerivative();
    double computeBetaDerivative(){ return 0;}
    double getDistance();
    double computeRatio(int particleNumber);
    void updateSlaterRelatedThings(int particleNumber);
    void updateSlaterMatrix(int particleNumber);
    void updateInverseSlaterMatrix(int particleNumber);
    void setupSlaterRelatedThings();
    std::vector<double> computeQuantumForce(int particleIndex, bool oldOrNew);

private:
    double m_omega = 1;
    double m_alpha = 1;

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
    double phi_00(int particleNumber);
    double phi_10(int particleNumber, int dimension);
    double phi_20(int particleNumber, int dimension);
    double phi_11(int particleNumber);
    
    std::vector<double> phi_00_der(int particleNumber);
    std::vector<double> phi_10_der(int particleNumber, int dimension);
    std::vector<double> phi_20_der(int particleNumber, int dimension);
    std::vector<double> phi_11_der(int particleNumber);

    double phi_00_double_der(int particleNumber);
    double phi_10_double_der(int particleNumber, int dimension);
    double phi_20_double_der(int particleNumber, int dimension);
    double phi_11_double_der(int particleNumber);

    void calculateInverseSlaterMatrix();

    double m_metropolisRatio = 1;
};
