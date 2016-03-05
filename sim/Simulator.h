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
  ~Simulator();
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

  // Discrete parameters.
  int cell_count[2]; // coarse cells in the x and y direction.
  double coarse_cell_dim; // dimension of the square coarsest cell.
  double dt_lattice; // discrete, Lattice Boltzmann time step.
  double viscosity_lattice; // lattice viscosity of the coarsest cells.
  double tau_coarse; // relaxation time of the coarsest cells.
  double omega_coarse; // relaxation frequency of the coarsest cells.
  double velocity_lattice; // lattice velocity of the coarsest cell.
  double buffer_viscosity_factor; // For the viscosity-counteracting approach. The multiple of the discrete viscosity to add as a buffer viscosity.

  // Initial values.
  // Initial distribution function values are set the equilibrium distribution values for the given macroscopic initial values.
  double rho0,u0,v0; // initial density and velocity.

  void read_settings(std::string filename);
  void process_settings();
  void create_grid();
};