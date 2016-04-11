#pragma once

// Holds all of the coarsest cells in the domain in row-major order.

#include <vector>
#include <string>
#include <iostream>
#include <omp.h>

#include "../cell/Cell.h"

class Grid 
{
public:
  const static std::size_t MAX_LEVELS = 32;
  std::vector<Cell> grid_levels[MAX_LEVELS]; //Holds the cells associated with each level. 0: coarsest, MAX_LEVEL-1: finest.
  double scale_factors[MAX_LEVELS] = { 0 };
  std::size_t cell_count[2]; // The number of cells in x and y direction on the coarsest level.
  // void initialize(int cell_count_x, int cell_count_y, 
  //   double rho0, double u0, double v0, 
  //   double tau, double omega, double nu, double nuc,
  //   char bc_[4], double bcv_[4], double uc_);
  // void iterate(size_t level);
  void initialize(int cell_count_x, int cell_count_y, 
    double rho0, double u0, double v0, 
    double tau, double omega, double nu, double nuc,
    char bc_[4], double U_, std::size_t relax_model_, std::size_t vc_model_);
  void iterate(std::size_t level);
  double get_max_velocity_magnitude() const; // Mainly for post-processing purposes.
  double get_min_velocity_magnitude() const; // Mainly for post-processing purposes.
private:
  double U = 0; // characteristic velocity.
  char bc[4] = {'w', 'w', 'w', 'w'}; // Boundary conditions for bottom, right, top, left, respectively.
  // double bcv[4] = {0, 0, 0, 0}; // Boundary condition values for bottom, right, top, left, respectively.
  std::size_t relax_model = 1;
  std::size_t vc_model = 0;
  void assign_coarse_neighbours();
  void enforce_bc_side(int side, char type, double value);
  void enforce_bc();
  void reconstruct_macro(int level);
  void bufferize_parallel(int level);
  void stream_parallel(int level);
  void calculate_scale_factors();
};