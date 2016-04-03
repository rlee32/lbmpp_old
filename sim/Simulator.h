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
  int timesteps;
  const double get_Re(){ return Re; }
  const double get_nu(){ return nu; }
  const double get_tau(){ return tau; }
  const double get_uc(){ return uc; }
  void output_coarse_field(std::string output_file_name);
private:
  // Run time control.
  bool refinement; // If true, enables solution-adaptive refinement.

  // Non-dimensional parameters.
  double Re; // Reynolds number.

  // Discrete parameters (for the coarsest cells!).
  double uc; // lattice velocity.
  double dt; // discrete, Lattice Boltzmann time step.
  double coarse_cell_size; // dimension of the square coarsest cell.
  double nu; // lattice viscosity of the coarsest cells.
  double nuc; // For the viscosity-counteracting approach. The buffer viscosity.
  double tau; // relaxation time of the coarsest cells.
  double omega; // relaxation frequency of the coarsest cells.

  // Grid temporary variables.
  int cell_count[2]; // coarse cells in the x and y direction.
  // Initial values.
  // Initial distribution function values are set the equilibrium distribution values for the given macroscopic initial values.
  double rho0,u0,v0; // initial density and velocity.
  // Boundary conditions
  char bc[4]; // bottom, right, top, left
  double bcv[4]; // values; bottom, right, top, left

  void read_settings(std::string filename);
  void process_settings();
};