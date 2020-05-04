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
    void setOneBodyDensityBins(int numberOfBins, double densityLength);
    
protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
    double calculatePositionSumSquared();

private:
    std::vector<std::vector<double> > m_oneBodyDensity = std::vector<std::vector<double> >(); 
    int m_numberOfBins = 50;
    double m_densityLength = 1.0;
};

