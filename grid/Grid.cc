#include "Grid.h"

using namespace std;

void Grid::iteration( size_t level )
{
  if (level < MAX_LEVELS-1)
  {
    if ( levels[level+1].get_active_cells() > 0 )
    {
      // cout << "Calling iteration from level " << level << endl; 
      iteration( level+1 );
      // cout << "Calling iteration from level " << level << endl; 
      iteration( level+1 );
    }
  }
  // cout << "Iterating Level "<< level << endl;
  levels[level].iteration( relax_model, vc_model );
}

void Grid::initialize(size_t cell_count_x, size_t cell_count_y, 
  double rho0, double u0, double v0, 
  double nu0, double nuc0,
  char sides[4], char bc[4], double U,
  double relax_model_, double vc_model_ )
{
  // Assign parameters.
  relax_model = relax_model_;
  vc_model = vc_model_;
  cell_count[0] = cell_count_x;
  cell_count[1] = cell_count_y;

  // Compute numerical parameters for each grid level.
  double scale_increase = 1;
  for (size_t i = 0; i < MAX_LEVELS; ++i)
  {
    GridLevel* next = ( i < MAX_LEVELS - 1 ) ? &levels[i+1] : nullptr;
    GridLevel* parent = ( i > 0 ) ? &levels[i-1] : nullptr;
    levels[i].initialize( scale_increase, nu0, nuc0, 
      sides, bc, U, 
      next, parent );
    scale_increase *= 2.0;
  }
  
  // Create and initialize coarse cells.
  // Boundary conditions are also defined here.
  Cell default_cell( rho0, u0, v0, 
    &levels[0].get_cells(), &levels[1].get_cells() );
  levels[0].create_coarse_grid( cell_count_x, cell_count_y, default_cell );

  // Testing refine operation.
  // levels[0].refine_all();
}

void Grid::set_coarse_solution(
  double rho0, vector<double>& u, vector<double>& v )
{
  vector<Cell>& c = levels[0].get_cells();
  // cout << c.size() << endl;
  for(size_t i = 0; i < c.size(); ++i)
  {
    c[i].reconstruct_distribution( rho0, u[i], v[i] );
  }
}

double Grid::max_mag() const
{
  double max = levels[0].max_mag();
  // cout << "max: " << max << endl;  
  for(size_t i = 1; i < MAX_LEVELS; ++i)
  {
    double test = levels[i].max_mag();
    if (test > max) max = test;
  }
  return max;
}

double Grid::min_mag() const
{
  double min = levels[0].min_mag();
  // cout << "min: " << min << endl;
  for(size_t i = 1; i < MAX_LEVELS; ++i)
  {
    double test = levels[i].min_mag();
    if (test < min) min = test;
  }
  return min;
}
double Grid::max_rho() const
{
  double max = levels[0].max_rho();
  for(size_t i = 1; i < MAX_LEVELS; ++i)
  {
    double test = levels[i].max_rho();
    if ( test > max and test > 0 ) max = test;
  }
  return max;
}

double Grid::min_rho() const
{
  double min = levels[0].min_rho();
  for(size_t i = 1; i < MAX_LEVELS; ++i)
  {
    double test = levels[i].min_rho();
    if ( test < min and test > 0 ) min = test;
  }
  return min;
}

double Grid::mag( size_t level, size_t cell_index) const
{
  return levels[level].mag(cell_index);
}

size_t Grid::active_cells() const
{
  size_t total = 0;
  for(size_t i = 0; i < MAX_LEVELS; ++i)
  {
    total += levels[i].get_active_cells();
  }
  // cout << total << endl;
  return total;
}

