#pragma once

// Holds all cells of same dimension, regardless of spatial proximity.

#include <vector>
#include <string>
#include <iostream>
#include <omp.h>

#include "../cell/Cell.h"

class GridLevel
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
  std::size_t total_cells() { return cells.size(); }
  // Mainly for post-processing purposes.
  double get_max_velocity_magnitude() const;
  double get_min_velocity_magnitude() const;
private:
  void assign_coarse_neighbours();
  void enforce_bc_side(int side, char type, double value);
  void enforce_bc();
  void reconstruct_macro(int level);
  void bufferize_parallel(int level);
  void stream_parallel(int level);
  void calculate_scale_factors();
};