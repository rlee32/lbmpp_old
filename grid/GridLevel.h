#pragma once

// Holds all cells of same dimension, 
// regardless of connectivity or spatial proximity.

#include <vector>
#include <string>
#include <iostream>
#include <omp.h>

#include "../cell/Cell.h"
#include "BoundaryConditions.h"

class GridLevel
{
public:
  // Performs one whole iteration at this grid level.
  void iteration( std::size_t relax_model, std::size_t vc_model );
  // Initialization
  void initialize( double scale_increase, double nu, double nuc, 
    char sides[4], char bc[4], double U );
  void create_coarse_grid( 
    std::size_t cell_count_x, std::size_t cell_count_y, Cell& default_cell );
  // Mainly for post-processing purposes.
  const Cell* cell( std::size_t index ) const { return &cells[index]; }
  std::size_t get_active_cells() const { return active_cells; }
  double max_mag() const;
  double min_mag() const;
  double mag(std::size_t cell_index) const;
  const std::vector<Cell>& get_cells() const { return cells; }
private:
  std::vector<Cell> cells;
  // Quantities scaled for this grid level.
  double scale_decrease = 1; // 2^-(tree depth)
  double scale_increase = 1; // 2^(tree depth)
  double tau = 0; // relaxation time.
  double omega = 0; // inverse of tau. relaxation frequency.
  double nu = 0; // viscosity.
  double nuc = 0; // viscosity buffer.
  //
  std::size_t active_cells = 0;
  BoundaryConditions bcs;
  void reconstruct_macro();
  void bufferize_parallel();
  void stream_parallel();
};
