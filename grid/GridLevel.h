#pragma once

// Holds all cells of same dimension, regardless of spatial proximity.

#include <vector>
#include <string>
#include <iostream>
#include <omp.h>

#include "../cell/Cell.h"

class Grid 
{
public:
  std::vector<Cell> cells;
  double scale = 1; // 2^(tree depth)
  double scale_inv = 1; // 2^-(tree depth)
  std::size_t tree_depth = 0; // 0: coarsest cells
  double tau = 0; // relaxation time.
  double omega = 0; // inverse of tau. relaxation frequency.
  double nu = 0; // viscosity.
  double nuc = 0; // viscosity buffer.
  void initialize(int cell_count_x, int cell_count_y, 
    double rho0, double u0, double v0, 
    double tau, double omega, double nu, double nuc,
    char bc_[4], double U_, std::size_t relax_model_, std::size_t vc_model_);
  void iterate(std::size_t level);
  double get_max_velocity_magnitude() const; // Mainly for post-processing purposes.
  double get_min_velocity_magnitude() const; // Mainly for post-processing purposes.
private:
  void assign_coarse_neighbours();
  void enforce_bc_side(int side, char type, double value);
  void enforce_bc();
  void reconstruct_macro(int level);
  void bufferize_parallel(int level);
  void stream_parallel(int level);
  void calculate_scale_factors();
};