#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void sampleAllEnergies(bool acceptedStep);
    void printOutputToTerminal();
    void printOutputToEnergyAlphaFile();
    void setFileOutput(int firstCriteria);
    void computeAverages();
    void evaluateNumberOfCyclesIncluded();
    void printOutputToEnergyFile();
    void printOneBodyDensityToFile();
    double getEnergy()          { return m_energy; }
    double getDerivative()      { return m_derivative; }
    double getAcceptance()      { return 100*(double)m_numberOfAcceptedSteps/(double)m_numberOfCyclesIncluded; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_energySquared = 0;
    double  m_derivative = 0;
    double  m_distance = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySquared = 0;
    double  m_cumulativeEnergyDerivative = 0;
    double  m_cumulativeAlphaDerivative = 0;
    double  m_cumulativeDistance = 0;
    int     m_numberOfAcceptedSteps = 0;
    int     m_numberOfCyclesIncluded = 0;
    int     m_firstCriteria = 1;
    std::vector <double> m_localEnergyVector = std::vector <double>();
    class System* m_system = nullptr;
};
