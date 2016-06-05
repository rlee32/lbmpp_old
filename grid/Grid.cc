#include "Grid.h"

using namespace std;


void Grid::experimental_initialize()
{
  // Testing refine operation.
  // levels[0].refine_all();
  // levels[0].refine_half( cell_count_x, cell_count_y );
  // levels[0].refine_three_parts( cell_count_x, cell_count_y );
  // levels[0].refine_three_parts_rotated( cell_count_x, cell_count_y );
  // levels[0].refine_three_parts_rotated_flipped( cell_count_x, cell_count_y );
  // cout << "Successful refinement" << endl;
  // levels[0].print_cell_status( cell_count_x, cell_count_y );
}


// void printdist(Cell& cell)
// {
//   // cout << "grid f: ";
//   // for(int i = 0; i < 8; ++i ) cout << cell.state.f[i] << " ";
//   // cout <<endl;
//   // cout << "grid b: ";
//   // for(int i = 0; i < 8; ++i ) cout << cell.state.b[i] << " ";
//   // cout <<endl;
// }
// void printneighbours(Cell& cell)
// {
//   cout << "n: ";
//   for(int i = 0; i < 8; ++i ) cout << cell.local.neighbours[i] << " ";
//   cout <<endl;
// }
// // static void printchildren(Cell& cell)
// // {
// //   cout << "children: ";
// //   for(int i = 0; i < 4; ++i ) cout << cell(i).local.me << " ";
// //   cout <<endl;
// // }
// static void printinterface(Cell& c)
// {
//   cout << "coal interface: ";
//   for(size_t i = 0; i < 8; ++i ) if(c.bc.coalesce[i]) cout << c.bc.coalesce[i] << " ";
//   cout <<endl;
// }
// const int ccc = 10;

void Grid::iteration( size_t level )
{
  // printchildren(levels[0].get_cell(84));
  // printchildren(levels[0].get_cell(85));
  // printchildren(levels[0].get_cell(94));
  // printchildren(levels[0].get_cell(95));
  // // printneighbours(levels[0].get_cell(94));
  // printneighbours(levels[1].get_cell(197));
  // printchildren(levels[0].get_cell(5));
  // printchildren(levels[0].get_cell(15));
  // printchildren(levels[0].get_cell(25));
  // printchildren(levels[0].get_cell(35));
  // printchildren(levels[0].get_cell(45));
  // printchildren(levels[0].get_cell(55));
  // printchildren(levels[0].get_cell(65));
  // printchildren(levels[0].get_cell(75));
  // printchildren(levels[0].get_cell(85));
  // printchildren(levels[0].get_cell(95));
  // printchildren(levels[0].get_cell(9));
  // printchildren(levels[0].get_cell(19));
  // printchildren(levels[0].get_cell(29));
  // printchildren(levels[0].get_cell(39));
  // printchildren(levels[0].get_cell(49));
  // printchildren(levels[0].get_cell(59));
  // printchildren(levels[0].get_cell(69));
  // printchildren(levels[0].get_cell(79));
  // printchildren(levels[0].get_cell(89));
  // printchildren(levels[0].get_cell(99));
  // printinterface(levels[0].get_cell(4));
  // printinterface(levels[0].get_cell(14));
  // printinterface(levels[0].get_cell(24));
  // printinterface(levels[0].get_cell(34));
  // printinterface(levels[0].get_cell(44));
  // printinterface(levels[0].get_cell(54));
  // printinterface(levels[0].get_cell(64));
  // printinterface(levels[0].get_cell(74));
  // printinterface(levels[0].get_cell(84));
  // printinterface(levels[0].get_cell(94));

  // cout << "child: " << (levels[0].get_cell(95))(1).local.me << endl; 
  // cout << "child: " << (levels[0].get_cell(95))(3).local.me << endl; 
  // printneighbours(levels[1].get_cell(1));
  // cout << "child: " << (levels[0].get_cell(5))(1).local.me << endl; 
  // cout << "Start level " << level << endl;
  // First determine if next level exists.
  bool go_to_next_level = false;
  if (level < MAX_LEVELS-1)
  {
    if ( levels[level+1].get_active_cells() > 0 ) go_to_next_level = true;
  }
  // Main recursive iteration.
  // cout << "Collide level " << level << endl;
  // printdist(levels[1].get_cell(ccc));
  levels[level].collide( relax_model, vc_model );
  // printdist(levels[1].get_cell(ccc));
  if ( go_to_next_level )
  {
    // cout << "Explode level " << level << endl;
    // printdist(levels[1].get_cell(ccc));
    levels[level].explode();
    // printdist(levels[1].get_cell(ccc));
    iteration( level+1 );
    iteration( level+1 );
  }
  // cout << "Stream level " << level << endl;
  // printdist(levels[1].get_cell(ccc));
  levels[level].stream();
  // printdist(levels[1].get_cell(ccc));
  // cout << "Coalesce level " << level << endl;
  // printdist(levels[1].get_cell(ccc));
  if ( go_to_next_level ) levels[level].coalesce();
  // printdist(levels[1].get_cell(ccc));
}

void Grid::initialize(size_t cell_count_x, size_t cell_count_y, 
  double rho0, double u0, double v0, 
  double nu0, double nuc0,
  char sides[4], char bc[4], double U,
  double relax_model_, double vc_model_, bool experimental_ )
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

  experimental = experimental_;
  if (experimental) experimental_initialize();
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

void Grid::reconstruct_macro()
{
  for(size_t i = 0; i < MAX_LEVELS; ++i) levels[i].reconstruct_macro();
}
