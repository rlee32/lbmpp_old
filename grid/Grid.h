#pragma once

// Holds all of the coarsest cells in the domain in row-major order.

#include <vector>
#include <string>

#include "../cell/Cell.h"


class Grid 
{
public:
  void initialize(int cell_count_x, int cell_count_y, 
    double rho0, double u0, double v0);
  std::vector<Cell> cells; // row-major order of cells.
  double get_max_velocity_magnitude(); // Mainly for post-processing purposes.
  double get_min_velocity_magnitude(); // Mainly for post-processing purposes.
  void enforce_bc(int side, char type, double value);
  int cell_count[2]; // The number of cells in x and y direction.
private:
  void assign_neighbours();
};