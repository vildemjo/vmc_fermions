#pragma once
#include <vector>
#include "../system.h"


class WaveFunction {
public:
    WaveFunction(class System* system);
    void updateOneBodyDensity();
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate() = 0;
    virtual double computeDoubleDerivative() = 0;
    virtual std::vector<double> computeDerivative(int particleIndex) = 0;
    virtual double computeAlphaDerivative() = 0;
    virtual double computeBetaDerivative() = 0;
    std::vector <std::vector<double> > getOneBodyDensity(){ return m_oneBodyDensity; };
    std::vector <std::vector<double> > getOneBodyDensityRadial(){ return m_oneBodyDensityRadial; };
    void setOneBodyDensityBins(int numberOfBins, double densityLength);
    virtual double getDistance() = 0;
    virtual double computeRatio(double newWaveFunction, double oldWaveFunction) = 0;
    virtual void updateSlaterRelatedThings(int particleNumber) = 0;
    virtual void updateSlaterMatrix(int particleNumber) = 0;
    virtual void updateInverseSlaterMatrix(int particleNumber) = 0;
    virtual void setupSlaterRelatedThings() = 0;
    
    virtual std::vector<double> computeQuantumForce(int particleIndex, bool oldOrNew) = 0;
    
protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
    double calculatePositionSumSquared();

private:
    std::vector<std::vector<double> > m_oneBodyDensity = std::vector<std::vector<double> >(); 
    std::vector<std::vector<double> > m_oneBodyDensityRadial = std::vector<std::vector<double> >(); 
    int m_numberOfBins = 50;
    double m_densityLength = 1.0;
};

