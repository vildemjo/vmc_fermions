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
    std::vector<double> getDerivative()      { return m_derivative; }
    double getKineticEnergy()     { return m_kineticEnergy; }
    double getPotentialEnergy()     { return m_potentialEnergy; }
    double getInteractionEnergy()     { return m_interactionEnergy; }
    double getAcceptance()      { return 100*(double)m_numberOfAcceptedSteps/(double)m_numberOfCyclesIncluded; }
    double getMeanDistance()     { return m_distance; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_energySquared = 0;
    double  m_derivativeAlpha = 0;
    double  m_derivativeBeta = 0;
    double  m_kineticEnergy = 0;
    double  m_potentialEnergy = 0;
    double  m_interactionEnergy = 0;
    std::vector<double>  m_derivative = std::vector<double>(2,1);
    double  m_distance = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySquared = 0;
    double  m_cumulativeEnergyDerivativeAlpha = 0;
    double  m_cumulativeEnergyDerivativeBeta = 0;
    double  m_cumulativeAlphaDerivative = 0;
    double  m_cumulativeBetaDerivative = 0;
    double  m_cumulativeDistance = 0;
    double  m_cumulativeKineticEnergy = 0;
    double  m_cumulativePotentialEnergy = 0;
    double  m_cumulativeInteractionEnergy = 0;
    int     m_numberOfAcceptedSteps = 0;
    int     m_numberOfCyclesIncluded = 0;
    int     m_firstCriteria = 1;
    std::vector <double> m_localEnergyVector = std::vector <double>();
    class System* m_system = nullptr;
    std::vector <double> m_densityVector = std::vector<double>();
};
