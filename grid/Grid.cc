#include "Grid.h"

using namespace std;

// void Grid::initialize(int cell_count_x, int cell_count_y, 
//   double rho0, double u0, double v0, 
//   double tau, double omega, double nu, double nuc,
//   char bc_[4], double bcv_[4], double uc_)
void Grid::initialize(int cell_count_x, int cell_count_y, 
  double rho0, double u0, double v0, 
  double tau, double omega, double nu, double nuc,
  char bc_[4], double U_)
{
  U = U_;
  vector<Cell>& cells = grid_levels[0];
  cell_count[0] = cell_count_x;
  cell_count[1] = cell_count_y;
  Cell default_cell( rho0, u0, v0,
    tau, omega, nu, nuc, 
    grid_levels );
  cells.resize(cell_count_x*cell_count_y, default_cell);
  assign_coarse_neighbours();
  for (size_t i = 0; i < 4; ++i) bc[i] = bc_[i];
  // for (size_t i = 0; i < 4; ++i) bcv[i] = bcv_[i];
  // cout << bcv[0] << endl;
  // cout << bcv[1] << endl;
  // cout << bcv[2] << endl;
  // cout << bcv[3] << endl;
}


void Grid::iterate(size_t level)
{
  // Check to see if this level needs to iterate.
  if (level >= MAX_LEVELS) return;
  // Collide, explode and stream all cells on current level.
  vector<Cell>& cells = grid_levels[level];
  if (cells.size() == 0) return;

  // 
  reconstruct_macro(level);
  #pragma omp parallel for
  for(uint i = 0; i < cells.size(); ++i)
  {
    cells[i].collide();
    cells[i].explode();
  }
  // enforce_bc();
  reconstruct_macro(level);
  stream_parallel(level);
  bufferize_parallel(level);
  iterate(level+1);
  iterate(level+1);
  reconstruct_macro(level);
  enforce_bc();
}

void Grid::stream_parallel(int level)
{
  #pragma omp parallel for
  for(uint i = 0; i < grid_levels[level].size(); ++i)
  {
    grid_levels[level][i].stream_parallel();
  }
}
void Grid::bufferize_parallel(int level)
{
  #pragma omp parallel for
  for(uint i = 0; i < grid_levels[level].size(); ++i)
  {
    grid_levels[level][i].bufferize_parallel();
  }
}
void Grid::reconstruct_macro(int level)
{
  #pragma omp parallel for
  for(uint i = 0; i < grid_levels[level].size(); ++i)
  {
    grid_levels[level][i].reconstruct_macro();
  }
}

void Grid::enforce_bc()
{
  enforce_bc_side('b',bc[0], U);
  enforce_bc_side('r',bc[1], U);
  enforce_bc_side('t',bc[2], U);
  enforce_bc_side('l',bc[3], U);
}

// Currently only works for coarsest level.
double Grid::get_max_velocity_magnitude() const
{
  const vector<Cell>& cells = grid_levels[0];
  double max = cells[0].get_velocity_magnitude();
  #pragma omp parallel for
  for(uint i = 1; i < cells.size(); ++i)
  {
    double test = cells[i].get_velocity_magnitude();
    if (test > max) max = test;
  }
  return max;
}
// Currently only works for coarsest level.
double Grid::get_min_velocity_magnitude() const
{
  const vector<Cell>& cells = grid_levels[0];
  double min = cells[0].get_velocity_magnitude();
  #pragma omp parallel for
  for(uint i = 1; i < cells.size(); ++i)
  {
    double test = cells[i].get_velocity_magnitude();
    if (test < min) min = test;
  }
  return min;
}

void Grid::assign_coarse_neighbours()
{
  vector<Cell>& cells = grid_levels[0];
  #pragma omp parallel for
  for (size_t i = 0; i < cell_count[0]; ++i)
  {
    for (size_t j = 0; j < cell_count[1]; ++j)
    {
      int ii = i + j*cell_count[0];
      bool right = i == cell_count[0]-1;
      bool top = j == 0;
      bool left = i == 0;
      bool bottom = j == cell_count[1]-1;
      int dj = -cell_count[0];
      // East
      if ( not right ) cells[ii].tree.neighbours[0] = &cells[ii+1];
      // Northeast
      if ( not right and not top ) cells[ii].tree.neighbours[1] = &cells[ii+1+dj];
      // North
      if ( not top ) cells[ii].tree.neighbours[2] = &cells[ii+dj];
      // Northwest
      if ( not top and not left ) cells[ii].tree.neighbours[3] = &cells[ii-1+dj];
      // West
      if ( not left ) cells[ii].tree.neighbours[4] = &cells[ii-1];
      // Southwest
      if ( not left and not bottom ) cells[ii].tree.neighbours[5] = &cells[ii-1-dj];
      // South
      if ( not bottom ) cells[ii].tree.neighbours[6] = &cells[ii-dj];
      // Southeast
      if ( not bottom and not right ) cells[ii].tree.neighbours[7] = &cells[ii-dj+1];
    }
  }
}

// side: b (bottom), r (right), t (top), l (left)
// type: w (wall), o (outlet), i (inlet), m (moving wall)
// value: applicable to m, i
// Currently only works for coarsest level.
void Grid::enforce_bc_side(int side, char type, double value)
{
  vector<Cell>& cells = grid_levels[0];
  int ii = -1;
  int dii = -1, djj = -1;
  int imax = -1;
  int ortho = -1; // the index pointing inwards, (unknown) orthogonal to boundary.
  int before = -1; // the index before ortho. In the case of 0, it is 7.
  double ubc = 0, vbc = 0;
  switch (side)
  {
    case 'b': // bottom
      // ii = 0;
      ii = (cell_count[1]-1)*cell_count[0];
      dii = 1;
      djj = -cell_count[0];
      // imax = cell_count[0] - 1;
      imax = cell_count[1]*cell_count[0] - 1;
      ortho = 2;
      break;
    case 'r': // right
      ii = cell_count[0]-1;
      dii = cell_count[0];
      djj = -1;
      imax = cell_count[1]*cell_count[0] - 1;
      ortho = 4;
      break;
    case 't': // top
      // ii = (cell_count[1]-1)*cell_count[0];
      ii = 0;
      dii = 1;
      djj = cell_count[0];
      // imax = cell_count[1]*cell_count[0] - 1;
      imax = cell_count[0] - 1;
      ortho = 6;
      break;
    case 'l': // left
      ii = 0;
      dii = cell_count[0];
      djj = 1;
      imax = (cell_count[1]-1)*cell_count[0];
      ortho = 0;
      break;
    default:
      break;
  }
  before = (ortho != 0) ? ortho-1 : 7;
  switch (type)
  {
    case 'w': // wall
      for (; ii <= imax; ii+=dii)
      {
        // meso
        cells[ii].state.f[ortho] = cells[ii].state.f[OPPOSITE(ortho)];
        cells[ii].state.f[before] = cells[ii].state.f[OPPOSITE(before)];
        cells[ii].state.f[ortho+1] = cells[ii].state.f[OPPOSITE(ortho+1)];
        // macro
        cells[ii].state.u = 0;
        cells[ii].state.v = 0;
      }
      break;
    case 'i': // inlet
      ubc = (side == 'b' or side == 't') ? 0 : value;
      vbc = (side == 'r' or side == 'l') ? 0 : value;
      for (ii+=dii; ii <= imax-dii; ii+=dii)
      {
        //meso
        int parallel = ( ortho + OPPOSITE(ortho) ) / 2.0;
        double coeff = 1.0 / ( 1.0 - value );
        double tangent = cells[ii].state.fc
          + cells[ii].state.f[parallel] 
          + cells[ii].state.f[OPPOSITE(parallel)];
        double out = cells[ii].state.f[OPPOSITE(ortho)]
          + cells[ii].state.f[OPPOSITE(ortho+1)]
          + cells[ii].state.f[OPPOSITE(before)];
        double density = coeff * ( tangent + 2.0 * out );
        double normal = 2.0 / 3.0 * value * density;
        double diagonal = value / 6.0 * density;
        cells[ii].state.f[ortho] = cells[ii].state.f[OPPOSITE(ortho)] + normal;
        cells[ii].state.f[ortho+1] = cells[ii].state.f[OPPOSITE(ortho+1)] + diagonal;
        cells[ii].state.f[before] = cells[ii].state.f[OPPOSITE(before)] + diagonal;
        // macro
        // cells[ii].state.rho = density;
        cells[ii].state.u = ubc;
        cells[ii].state.v = vbc;
      }
      break;
    case 'o': // outlet (zero-gradient extrapolation)
      for (ii+=dii; ii <= imax-dii; ii+=dii)
      {
        cells[ii].state.f[ortho] = cells[ii+djj].state.f[ortho];
        cells[ii].state.f[ortho+1] = cells[ii+djj].state.f[ortho+1];
        cells[ii].state.f[before] = cells[ii+djj].state.f[before];
      }
      break;
    case 'm': // moving wall
      ubc = (side == 'b' or side == 't') ? value : 0;
      vbc = (side == 'r' or side == 'l') ? value : 0;
      for (ii+=dii; ii <= imax-dii; ii+=dii)
      {
        // meso
        int sign = (side == 't' or side == 'l') ? 1.0: -1.0;
        int parallel = ( ortho + OPPOSITE(ortho) ) / 2.0;
        double coeff = 1.0 / ( 1.0 - value );
        double tangent = cells[ii].state.fc
          + cells[ii].state.f[parallel] 
          + cells[ii].state.f[OPPOSITE(parallel)];
        double out = cells[ii].state.f[OPPOSITE(ortho)]
          + cells[ii].state.f[OPPOSITE(ortho+1)]
          + cells[ii].state.f[OPPOSITE(before)];
        double density = coeff * ( tangent + 2.0 * out );
        double diagonal = value / 6.0 * density;
        cells[ii].state.f[ortho] = cells[ii].state.f[OPPOSITE(ortho)];
        cells[ii].state.f[ortho+1] = cells[ii].state.f[OPPOSITE(ortho+1)] 
          + sign * diagonal;
        cells[ii].state.f[before] = cells[ii].state.f[OPPOSITE(before)] 
          - sign * diagonal;
        // macro
        // cells[ii].state.rho = density;
        cells[ii].state.u = ubc;
        cells[ii].state.v = vbc;
      }
      break;
    default:
      break;
  }
}

void Grid::calculate_scale_factors()
{
  double mult = 1;
  for (size_t i = 0; i < MAX_LEVELS; ++i)
  {
    scale_factors[i] = 1.0 / mult;
    mult *= 2;
  }
}