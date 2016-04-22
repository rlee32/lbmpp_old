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
  }

  // cout << "BC" << endl;
  bcs.apply_bc();

  // cout << "Stream" << endl;
  // Stream
  stream_parallel();
  bufferize_parallel();
  
  // cout << "Finish stream" << endl;
  refresh_active_cells();
  // cout << active_cells << endl;
}

void GridLevel::refresh_active_cells()
{
  active_cells = compute_active_cells();
}


size_t GridLevel::compute_active_cells() const
{
  size_t active = 0;
  // #pragma omp parallel for reduction(+:active)
  for(size_t i = 0; i < cells.size(); ++i)
  {
    // cout << cells[i].state.active << endl;
    active += cells[i].state.active;
  }
  return active;
}

void GridLevel::stream_parallel()
{
  #pragma omp parallel for
  for(size_t i = 0; i < cells.size(); ++i)
  {
    cells[i].stream_parallel( cells );
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
  char sides[4], char bc[4], double U,
  GridLevel* next_grid_level_ )
{
  bcs.initialize( sides, bc, U, next_grid_level_->get_bcs(), &cells );
  // compute scale factors.
  scale_decrease = 1.0 / scale_increase;
  nu = scale_decrease * nu0;
  nuc = scale_decrease * nuc0;
  tau = 3 * (nu + nuc) + 0.5;
  omega = 1 / tau;
  // cout << nu << " " << nuc << " " << tau << " " << omega << endl;
  next_grid_level = next_grid_level_;
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
      cells[ii].local.neighbours[0] = ii+1;
    }
  }
  // West neighbours
  for (size_t i = 1; i < cell_count_x; ++i)
  {
    for (size_t j = 0; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].local.neighbours[4] = ii-1;
    }
  }
  // North neighbours
  for (size_t i = 0; i < cell_count_x; ++i)
  {
    for (size_t j = 0; j < cell_count_y-1; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].local.neighbours[2] = ii + cell_count_x;
    }
  }
  // South neighbours
  for (size_t i = 0; i < cell_count_x; ++i)
  {
    for (size_t j = 1; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].local.neighbours[6] = ii - cell_count_x;
    }
  }
  // Northeast neighbours
  for (size_t i = 0; i < cell_count_x-1; ++i)
  {
    for (size_t j = 0; j < cell_count_y-1; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].local.neighbours[1] = ii + 1 + cell_count_x;
    }
  }
  // Northwest neighbours
  for (size_t i = 1; i < cell_count_x; ++i)
  {
    for (size_t j = 0; j < cell_count_y-1; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].local.neighbours[3] = ii - 1 + cell_count_x;
    }
  }
  // Southwest neighbours
  for (size_t i = 1; i < cell_count_x; ++i)
  {
    for (size_t j = 1; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].local.neighbours[5] = ii - 1 - cell_count_x;
    }
  }
  // Southeast neighbours
  for (size_t i = 0; i < cell_count_x-1; ++i)
  {
    for (size_t j = 1; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      cells[ii].local.neighbours[7] = ii + 1 - cell_count_x;
    }
  }

  // Top and bottom near-wall nodes.
  for (size_t i = 0; i < cell_count_x; ++i)
  {
    // bottom adjacent
    size_t ii = i;
    cells[ii].local.nn[3] = 0;
    cells[ii].local.fully_interior_cell = false;
    // bottom one away
    size_t j = 1;
    ii = i + j*cell_count_x;
    cells[ii].local.nn[3] = 1;
    cells[ii].local.fully_interior_cell = false;
    // top one away
    j = cell_count_y-2;
    ii = i + j*cell_count_x;
    cells[ii].local.nn[1] = 1;
    cells[ii].local.fully_interior_cell = false;
    // top adjacent
    j = cell_count_y-1;
    ii = i + j*cell_count_x;
    cells[ii].local.nn[1] = 0;
    cells[ii].local.fully_interior_cell = false;
  }
  // Left and right near-wall nodes.
  for (size_t j = 0; j < cell_count_y; ++j)
  {
    // left adjacent
    size_t i = 0;
    size_t ii = i + j*cell_count_x;
    cells[ii].local.nn[2] = 0;
    cells[ii].local.fully_interior_cell = false;
    // left one away
    i = 1;
    ii = i + j*cell_count_x;
    cells[ii].local.nn[2] = 1;
    cells[ii].local.fully_interior_cell = false;
    // right one away
    i = cell_count_x-2;
    ii = i + j*cell_count_x;
    cells[ii].local.nn[0] = 1;
    cells[ii].local.fully_interior_cell = false;
    // right adjacent
    i = cell_count_x-1;
    ii = i + j*cell_count_x;
    cells[ii].local.nn[0] = 0;
    cells[ii].local.fully_interior_cell = false;
  }
  // Now to identify boundaries.
  for(size_t i = 0; i < cell_count_x; ++i)
  {
    size_t ii = i;
    bcs.add_cell( ii, 'b' );
    size_t j = cell_count_y - 1;
    ii += j * cell_count_x;
    bcs.add_cell( ii, 't' );
  }
  for(size_t j = 0; j < cell_count_y; ++j)
  {
    size_t ii = j * cell_count_x;
    bcs.add_cell( ii, 'l' );
    ii += cell_count_x - 1;
    bcs.add_cell( ii, 'r' );
  }
  refresh_active_cells();
}

// Checks all cells to see if they have been marked for linking 
//  and links them (with neighbours).
// Meant to be called after entire-grid refinement.
void GridLevel::link_marked()
{
  // Can parallelize
  for ( size_t i = 0; i < cells.size(); ++i )
  {
    if ( cells[i].action.link_children )
    {
      // cout << "Linking dawg" << endl;
      cells[i].link_children( cells, next_grid_level->get_cells() );
    }
  }
}
// Checks all cells to see if they have been marked for refinement 
//  and refines them.
void GridLevel::refine_marked()
{
  // Do not parallelize without using locks!!!
  for (size_t i = 0; i < cells.size(); ++i)
  {
    if ( cells[i].action.refine ) 
    {
      cells[i].refine( next_grid_level->get_cells() );
    }
  }
}
// This is a test function to be called during initialization.
void GridLevel::refine_all()
{
  for(size_t i = 0; i < cells.size(); ++i) cells[i].action.refine = true;
  refine_marked();
  link_marked();
  refined_cell_bc();
  active_cells = compute_active_cells();
  next_grid_level->refresh_active_cells();
  for(size_t i = 0; i < cells.size(); ++i) cells[i].action.refine = false;
}
// This is a test function to be called during initialization.
void GridLevel::refine_half( size_t i_cells, size_t j_cells )
{
  for(size_t j = 0; j < j_cells; ++j)
  {
    for(size_t i = 0; i < i_cells/2; ++i)
    {
      size_t ii = i + j*i_cells;
      cells[ii].action.refine = true;
    }
  }
  // These are end-of-whole-grid-iteration operations.
  refine_marked();
  link_marked();
  refined_cell_bc();
  active_cells = compute_active_cells();
  next_grid_level->refresh_active_cells();
  for(size_t i = 0; i < cells.size(); ++i) cells[i].action.refine = false;
}

// Call AFTER refinement (via action.refine), and BEFORE linking children.
// Goes through all cells and identifies newly-created cells that 
//  have parent neighbours but no corresponding same-level neighbours, 
//  and creates the neighbours.
void GridLevel::identify_interfaces()
{
  for( size_t i = 0; i < cells.size(); ++i )
  {
    if( cells[i].action.refine )
    {
      for ( size_t n = 0; n < 8; ++n )
      {
        if ( cells[i].has_neighbour(n) )
        {
          Cell& nc = get_neighbour(i,n);
          if ( not nc.action.refine )
          {
            if ( nc.has_children() )
            {
              if (not nc.has_interface_children(next_grid_level->get_cells()))
              {
                // Set children as interface.
                next_grid_level->set_interface(nc.local.children[0]);
                next_grid_level->set_interface(nc.local.children[1]);
                next_grid_level->set_interface(nc.local.children[2]);
                next_grid_level->set_interface(nc.local.children[3]);
              }
            }
            else
            {
              // Create children as interface.
              nc.create_interface_children( next_grid_level->get_cells() );
            }
          }
        }
      }
    }
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


