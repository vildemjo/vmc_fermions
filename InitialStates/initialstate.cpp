#include "initialstate.h"
#include <cmath>
#include "../particle.h"
#include "../system.h"
#include <cassert>
#include "../Hamiltonians/hamiltonian.h"
#include <iostream>

InitialState::InitialState(System* system) {
    m_system = system;
}
