#pragma once

// Holds all of the coarsest cells in the domain in row-major order.

#include <vector>
#include <string>
#include <iostream>
#include <omp.h>

#include "../cell/Cell.h"
#include "GridLevel.h"

class Grid 
{
public:
  void iteration( std::size_t level );
  std::size_t active_cells() const;
  std::vector<Cell>& get_cells(std::size_t index) 
    { return levels[index].get_cells(); }
  //
  void initialize( std::size_t cell_count_x, std::size_t cell_count_y, 
    double rho0, double u0, double v0, 
    double nu, double nuc,
    char sides[4], char bc[4], double U,
    double relax_model, double vc_model );
  void set_coarse_solution(
    double rho0, std::vector<double>& u, std::vector<double>& v );
  // 
  std::size_t cell_count_x() { return cell_count[0]; }
  std::size_t cell_count_y() { return cell_count[1]; }
  // Post-processing functions.
  double max_mag() const;
  double min_mag() const;
  double max_rho() const;
  double min_rho() const;
  double mag( std::size_t level, std::size_t cell_index ) const;
private:
  const static std::size_t MAX_LEVELS = 32;
  GridLevel levels[MAX_LEVELS];
  double relax_model = 1;
  double vc_model = 0;
  std::size_t cell_count[2] = { 0, 0 };// Coarsest level cell dimension.
};