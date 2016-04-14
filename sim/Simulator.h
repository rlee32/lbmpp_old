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
  void iteration(); // performs one time step advancement.
  Grid grid;
  std::size_t get_cell_count_0() const { return cell_count[0]; }
  std::size_t get_cell_count_1() const { return cell_count[1]; }
  double get_Re() const { return Re; }
  double get_nu() const { return nu; }
  double get_tau() const { return tau; }
  double get_M() const { return M; }
  double get_U() const { return U; }
  double get_timesteps() const { return timesteps; }
  double get_relax_model() const { return relax_model; }
  double get_vc_model() const { return vc_model; }
  double get_nucf() const { return nuc / nu; }
  void output_coarse_field( std::string output_file_name );
private:
  // Run time control.
  bool refinement = false; // If true, enables solution-adaptive refinement.
  std::size_t timesteps = 0;
  // Solver parameters.
  std::size_t relax_model = 1;
  std::size_t vc_model = 0;
  // Physical parameters.
  double Re = 0; // Reynolds number.
  double M = 0; // Mach number for characteristic velocity.
  double U = 0; // characteristic velocity, used for velocity BC values.
  double nu = 0; // lattice viscosity of the coarsest cells.
  double L = 0;
  double tau = 0; // relaxation time of the coarsest cells.
  double omega = 0; // relaxation frequency of the coarsest cells.
  // Initial values.
  double rho0 = 1;
  double u0 = 0;
  double v0 = 0;
  // For the viscosity-counteracting approach. The buffer viscosity.
  double nuc = 0;
  double nucf = 0;
  // Boundary conditions
  char bc[4] = { 'w', 'w', 'w', 'w' }; // bottom, right, top, left
  std::string face_order[4] = { "bottom", "right", "top", "left" };
  char face_order_char[4] = { 'b', 'r', 't', 'l' };
  // coarse grid dimension
  std::size_t cell_count[2] = { 0, 0 };

  void read_settings(std::string filename);
  void process_settings();
};