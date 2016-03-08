#pragma once

// Holds all of the coarsest cells in the domain in row-major order.

#include <vector>
#include <string>
#include <iostream>

#include "../cell/Cell.h"

#define MAX_LEVELS 32

class Grid 
{
public:
  std::vector<Cell> grid_levels[MAX_LEVELS]; //Holds the cells associated with each level. 0: coarsest, MAX_LEVEL-1: finest.
  int cell_count[2]; // The number of cells in x and y direction on the coarsest level.
  void initialize(int cell_count_x, int cell_count_y, 
    double rho0, double u0, double v0, 
    double tau, double omega, double nu, double nuc);
  void enforce_coarse_bc(int side, char type, double value);
  void iterate(int level);
  double get_max_velocity_magnitude(); // Mainly for post-processing purposes.
  double get_min_velocity_magnitude(); // Mainly for post-processing purposes.
private:
  void assign_coarse_neighbours();
};