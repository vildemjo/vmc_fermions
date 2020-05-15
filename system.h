#pragma once
#include <vector>
#include <string>
#include <Math/random.h>

class System {
public:
    System();
    System(int seed);

    bool metropolisStep             ();
    bool metropolisStepImportance   ();
    void runMetropolisSteps         (int numberOfMetropolisSteps, int firstCriteria, bool importanceOrNot, bool allEnergiesOrNot, double stepLength);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibration           (double equilibration);
    void setFrequency               (double omega);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void setAnalytical              (bool statement);
    void setImportance              (bool statement);
    void setAllEnergies             (bool statement);
    void setFileName                (std::string filename) {m_filename = filename;}
    void setSpinFactor              ();
    class Random*                    getRandomEngine()   { return m_random; }
    class WaveFunction*              getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*               getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                   getSampler()        { return m_sampler; }
    class InitialState*              getInitialState()   { return m_initialState; }
    std::vector<class Particle*>     getParticles()      { return m_particles; }
    std::string getFileName()                            { return m_filename; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    int getSteps()                      { return m_steps; }
    bool getAnalytical()                { return m_analytical; }
    bool getImportance()                { return m_importance; }
    bool getAllEnergies()               { return m_allEnergies;}
    double getStepLength()              { return m_stepLength; }
    double getEquilibration()           { return m_equilibration; }
    std::vector<std::vector<double>> getSpinFactor()              { return m_spinFactor; }
    double getFrequency()               { return m_omega; }


    double greensFunctionFraction(std::vector<double> posNew, std::vector<double> posOld, std::vector<double> forceNew, std::vector<double> forceOld);

private:
    int                              m_numberOfParticles = 0;
    int                              m_numberOfDimensions = 0;
    int                              m_numberOfMetropolisSteps = 0;
    int                              m_steps = 0;
    bool                             m_analytical = false;
    bool                             m_importance = false;
    bool                             m_allEnergies = false;
    double                           m_diffConstant = 0.5;
    int                              m_equilibration = 0;
    double                           m_stepLength = 0.1;
    double                           m_timeStep   = 0.1;
    std::vector<std::vector<double>> m_spinFactor = std::vector<std::vector<double>>();
    double                           m_omega = 1.0;
    std::string                      m_filename = "Output/no_name_specified";
    class WaveFunction*              m_waveFunction = nullptr;
    class Hamiltonian*               m_hamiltonian = nullptr;
    class InitialState*              m_initialState = nullptr;
    class Sampler*                   m_sampler = nullptr;
    class Random*                    m_random = nullptr;
    std::vector<class Particle*>     m_particles = std::vector<class Particle*>(); 
};

