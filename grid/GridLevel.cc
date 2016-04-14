#include "GridLevel.h"

using namespace std;

// Collide, explode and stream all cells on this grid level.
void GridLevel::iteration( std::size_t relax_model, std::size_t vc_model )
{
  reconstruct_macro();
  
  // Collide
  #pragma omp parallel for
  for(size_t i = 0; i < cells.size(); ++i)
  {
    cells[i].collide(
      relax_model, vc_model, omega, scale_increase, scale_decrease, nuc);
    cells[i].explode();
  }
  
  bcs.apply_bc();
  
  // Stream
  stream_parallel();
  bufferize_parallel();
  
}

void GridLevel::stream_parallel()
{
  #pragma omp parallel for
  for(size_t i = 0; i < cells.size(); ++i)
  {
    cells[i].stream_parallel();
  }
}

void GridLevel::bufferize_parallel()
{
  #pragma omp parallel for
  for(size_t i = 0; i < cells.size(); ++i)
  {
    cells[i].bufferize_parallel();
  }
}

void GridLevel::reconstruct_macro()
{
  #pragma omp parallel for
  for(size_t i = 0; i < cells.size(); ++i)
  {
    cells[i].reconstruct_macro();
  }
}

void GridLevel::initialize( double scale_increase, 
  double nu0, double nuc0, 
  char sides[4], char bc[4], double U )
{
  bcs.initialize( sides, bc, U );
  // compute scale factors.
  scale_decrease = 1.0 / scale_increase;
  nu = scale_decrease * nu0;
  nuc = scale_decrease * nuc0;
  tau = 3 * (nu + nuc) + 0.5;
  omega = 1 / tau;
}

void GridLevel::create_coarse_grid( size_t cell_count_x, size_t cell_count_y,
  Cell& default_cell )
{
  active_cells = cell_count_x * cell_count_y;
  cells.resize( active_cells, default_cell );
  
  // East neighbours
  for (size_t i = 0; i < cell_count_x-1; ++i)
  {
    for (size_t j = 0; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].tree.neighbours[0] = &cells[ii+1];
    }
  }
  // West neighbours
  for (size_t i = 1; i < cell_count_x; ++i)
  {
    for (size_t j = 0; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].tree.neighbours[4] = &cells[ii-1];
    }
  }
  // North neighbours
  for (size_t i = 0; i < cell_count_x; ++i)
  {
    for (size_t j = 0; j < cell_count_y-1; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].tree.neighbours[2] = &cells[ii + cell_count_x];
    }
  }
  // South neighbours
  for (size_t i = 0; i < cell_count_x; ++i)
  {
    for (size_t j = 1; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].tree.neighbours[6] = &cells[ii - cell_count_x];
    }
  }
  // Northeast neighbours
  for (size_t i = 0; i < cell_count_x-1; ++i)
  {
    for (size_t j = 0; j < cell_count_y-1; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].tree.neighbours[1] = &cells[ii + 1 + cell_count_x];
    }
  }
  // Northwest neighbours
  for (size_t i = 1; i < cell_count_x; ++i)
  {
    for (size_t j = 0; j < cell_count_y-1; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].tree.neighbours[3] = &cells[ii - 1 + cell_count_x];
    }
  }
  // Southwest neighbours
  for (size_t i = 1; i < cell_count_x; ++i)
  {
    for (size_t j = 1; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].tree.neighbours[5] = &cells[ii - 1 - cell_count_x];
    }
  }
  // Southeast neighbours
  for (size_t i = 0; i < cell_count_x-1; ++i)
  {
    for (size_t j = 1; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].tree.neighbours[7] = &cells[ii + 1 - cell_count_x];
    }
  }

  // Top and bottom near-wall nodes.
  for (size_t i = 0; i < cell_count_x; ++i)
  {
    // bottom adjacent
    size_t ii = i;
    cells[ii].tree.nn[3] = 0;
    cells[ii].tree.fully_interior_cell = false;
    // bottom one away
    size_t j = 1;
    ii = i + j*cell_count_x;
    cells[ii].tree.nn[3] = 1;
    cells[ii].tree.fully_interior_cell = false;
    // top one away
    j = cell_count_y-2;
    ii = i + j*cell_count_x;
    cells[ii].tree.nn[1] = 1;
    cells[ii].tree.fully_interior_cell = false;
    // top adjacent
    j = cell_count_y-1;
    ii = i + j*cell_count_x;
    cells[ii].tree.nn[1] = 0;
    cells[ii].tree.fully_interior_cell = false;
  }
  // Left and right near-wall nodes.
  for (size_t j = 0; j < cell_count_y; ++j)
  {
    // left adjacent
    size_t i = 0;
    size_t ii = i + j*cell_count_x;
    cells[ii].tree.nn[2] = 0;
    cells[ii].tree.fully_interior_cell = false;
    // left one away
    i = 1;
    ii = i + j*cell_count_x;
    cells[ii].tree.nn[2] = 1;
    cells[ii].tree.fully_interior_cell = false;
    // right one away
    i = cell_count_x-2;
    ii = i + j*cell_count_x;
    cells[ii].tree.nn[0] = 1;
    cells[ii].tree.fully_interior_cell = false;
    // right adjacent
    i = cell_count_x-1;
    ii = i + j*cell_count_x;
    cells[ii].tree.nn[0] = 0;
    cells[ii].tree.fully_interior_cell = false;
  }
  // Now to identify boundaries.
  for(size_t i = 0; i < cell_count_x; ++i)
  {
    size_t ii = i;
    bcs.add_cell( &cells[ii], 'b' );
    size_t j = cell_count_y - 1;
    ii += j * cell_count_x;
    bcs.add_cell( &cells[ii], 't' );
  }
  for(size_t j = 0; j < cell_count_y; ++j)
  {
    size_t ii = j * cell_count_x;
    bcs.add_cell( &cells[ii], 'l' );
    ii += cell_count_x - 1;
    bcs.add_cell( &cells[ii], 'r' );
  }
}

double GridLevel::max_mag() const
{
  double max = 0;
  if (cells.size() > 0)
  {
    max = cells[0].get_mag();
    for(size_t i = 1; i < cells.size(); ++i)
    {
      double test = cells[i].get_mag();
      if (test > max) max = test;
    }
  }
  return max;
}

double GridLevel::min_mag() const
{
  double min = 0;
  if (cells.size() > 0)
  {
    min = cells[0].get_mag();
    for(size_t i = 1; i < cells.size(); ++i)
    {
      double test = cells[i].get_mag();
      if (test < min) min = test;
    }
  }
  return min;
}

double GridLevel::max_rho() const
{
  double min = -1;
  if (cells.size() > 0)
  {
    min = cells[0].rho();
    for(size_t i = 1; i < cells.size(); ++i)
    {
      double test = cells[i].rho();
      if (test > min) min = test;
    }
  }
  return min;
}
double GridLevel::min_rho() const
{
  double min = -1;
  if (cells.size() > 0)
  {
    min = cells[0].rho();
    for(size_t i = 1; i < cells.size(); ++i)
    {
      double test = cells[i].rho();
      if (test < min) min = test;
    }
  }
  return min;
}

double GridLevel::mag( size_t cell_index ) const
{
  return cells[cell_index].get_mag();
}


