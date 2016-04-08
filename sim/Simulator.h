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
  const double get_Re() const { return Re; }
  const double get_nu() const { return nu; }
  const double get_tau() const { return tau; }
  const double get_U() const { return U; }
  const double get_timesteps() const { return timesteps; }
  void output_coarse_field(std::string output_file_name);
private:
  // Run time control.
  bool refinement = false; // If true, enables solution-adaptive refinement.

  // Non-dimensional parameters.
  double Re = 0; // Reynolds number.

  // Discrete parameters (for the coarsest cells!).
  double U = 0; // characteristic velocity, used for velocity BC values.
  double nu = 0; // lattice viscosity of the coarsest cells.
  double nuc = 0; // For the viscosity-counteracting approach. The buffer viscosity.
  double tau = 0; // relaxation time of the coarsest cells.
  double omega = 0; // relaxation frequency of the coarsest cells.

  // Grid temporary variables.
  int cell_count[2] = { 0, 0 }; // coarse cells in the x and y direction.
  // Initial values.
  // Initial distribution function values are set the equilibrium distribution values for the given macroscopic initial values.
  double rho0 = 0, u0 = 0, v0 = 0; // initial density and velocity.
  // Boundary conditions
  char bc[4] = { 'w', 'w', 'w', 'w' }; // bottom, right, top, left
  // double bcv[4] = {  }; // values; bottom, right, top, left

  std::size_t timesteps;

  void read_settings(std::string filename);
  void process_settings();
};