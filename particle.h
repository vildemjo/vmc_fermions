#pragma once
#include <vector>

class Particle {
public:
    Particle();
    void setPosition(const std::vector<double> &position);
    void adjustPosition(double change, int dimension);
    void setNumberOfDimensions(int numberOfDimensions);
    std::vector<double> getPosition() { return m_position; }
    void setParticleIndex(int index) { m_particleIndex = index; }

private:
    int     m_numberOfDimensions = 0;
    std::vector<double> m_position = std::vector<double>();
    int m_particleIndex = 0;
    class System* m_system = nullptr;
};

