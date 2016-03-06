#pragma once

// Contains settings for the simulation.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../grid/Grid.h"

class Simulator
{ 
public:
  Simulator(std::string filename);
  void iterate(); // performs one time step advancement.
  Grid grid;
private:
  // Run time control.
  int timesteps;
  bool refinement; // If true, enables solution-adaptive refinement.

  // Physical parameters.
  double viscosity_physical;
  double velocity_physical;
  double length_physical; 
  double dt_physical; // physical time scale derived from length and velocity.

  // Non-dimensional parameters.
  double Re; // Reynolds number.

  // Discrete parameters (for the coarsest cells!).
  double dt_lattice; // discrete, Lattice Boltzmann time step.
  double coarse_cell_size; // dimension of the square coarsest cell.
  double viscosity_lattice; // lattice viscosity of the coarsest cells.
  double tau; // relaxation time of the coarsest cells.
  double omega; // relaxation frequency of the coarsest cells.
  double velocity_lattice; // lattice velocity of the coarsest cell.
  double buffer_viscosity_factor; // For the viscosity-counteracting approach. The multiple of the discrete viscosity to add as a buffer viscosity.

  void read_settings(std::string filename);
  void process_settings();
};