#include "Grid.h"

using namespace std;

void Grid::initialize(int cell_count_x, int cell_count_y, 
  double rho0, double u0, double v0)
{
  vector<Cell>& cells = grid_levels[0];
  cell_count[0] = cell_count_x;
  cell_count[1] = cell_count_y;
  Cell default_cell(rho0,u0,v0);
  cells.resize(cell_count_x*cell_count_y, default_cell);
  assign_coarse_neighbours();
}


void Grid::iterate(int level)
{
  if (level >= MAX_LEVELS) return;
  // Collide, explode and stream all cells on current level.
  vector<Cell>& cells = grid_levels[level];
  if (cells.size() == 0) return;
  for(vector<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    it->ces();
  }
  iterate(level+1);
  iterate(level+1);
  for(vector<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    it->coalesce();
  }
}


// Currently only works for coarsest level.
double Grid::get_max_velocity_magnitude()
{
  vector<Cell>& cells = grid_levels[0];
  double max = cells[0].get_velocity_magnitude();
  for(vector<Cell>::iterator it = cells.begin()+1; it != cells.end(); ++it)
  {
    double test = it->get_velocity_magnitude();
    if (test > max) max = test;
  }
  return max;
}
// Currently only works for coarsest level.
double Grid::get_min_velocity_magnitude()
{
  vector<Cell>& cells = grid_levels[0];
  double min = cells[0].get_velocity_magnitude();
  for(vector<Cell>::iterator it = cells.begin()+1; it != cells.end(); ++it)
  {
    double test = it->get_velocity_magnitude();
    if (test < min) min = test;
  }
  return min;
}

void Grid::assign_coarse_neighbours()
{
  vector<Cell>& cells = grid_levels[0];
  for (int i = 0; i < cell_count[0]; ++i)
  {
    for (int j = 0; j < cell_count[1]; ++j)
    {
      int ii = i + j*cell_count[0];
      bool right = i == cell_count[0]-1;
      bool top = j == 0;
      bool left = i == 0;
      bool bottom = j == cell_count[1]-1;
      // East
      if ( not right ) cells[ii].tree.neighbours[0] = &cells[ii+1];
      // Northeast
      if ( not right and not top ) cells[ii].tree.neighbours[1] = &cells[ii+1-j];
      // North
      if ( not top ) cells[ii].tree.neighbours[2] = &cells[ii-j];
      // Northwest
      if ( not top and not left ) cells[ii].tree.neighbours[3] = &cells[ii-1-j];
      // West
      if ( not left ) cells[ii].tree.neighbours[4] = &cells[ii-1];
      // Southwest
      if ( not left and not bottom ) cells[ii].tree.neighbours[5] = &cells[ii-1+j];
      // South
      if ( not bottom ) cells[ii].tree.neighbours[6] = &cells[ii+j];
      // Southeast
      if ( not bottom and not right ) cells[ii].tree.neighbours[7] = &cells[ii+j+1];
    }
  }
}

// side: b (bottom), r (right), t (top), l (left)
// type: w (wall), o (outlet), i (inlet), m (moving wall)
// value: applicable to m, i
// Currently only works for coarsest level.
void Grid::enforce_coarse_bc(int side, char type, double value)
{
  vector<Cell>& cells = grid_levels[0];
  int ii = -1;
  int dii = -1, djj = -1;
  int imax = -1;
  int ortho = -1; // the index pointing inwards, (unknown) orthogonal to boundary.
  int before = -1; // the index before ortho. In the case of 0, it is 7.
  switch (side)
  {
    case 'b': // bottom
      ii = 0;
      dii = 1;
      djj = -cell_count[0];
      imax = cell_count[0] - 1;
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
      ii = (cell_count[1]-1)*cell_count[0];
      dii = 1;
      djj = cell_count[0];
      imax = cell_count[1]*cell_count[0] - 1;
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
        cells[ii].state.f[ortho] = cells[ii].state.f[OPPOSITE(ortho)];
        cells[ii].state.f[before] = cells[ii].state.f[OPPOSITE(before)];
        cells[ii].state.f[ortho+1] = cells[ii].state.f[OPPOSITE(ortho+1)];
      }
      break;
    case 'i': // inlet
      for (; ii <= imax; ii+=dii)
      {
        int parallel = ( ortho + OPPOSITE(ortho) ) / 2;
        double coeff = 1.0 / ( 1.0 - value );
        double tangent = cells[ii].state.fc
          + cells[ii].state.f[parallel] 
          + cells[ii].state.f[OPPOSITE(parallel)];
        double out = cells[ii].state.f[OPPOSITE(ortho)]
          + cells[ii].state.f[OPPOSITE(ortho+1)]
          + cells[ii].state.f[OPPOSITE(before)];
        double density = coeff * ( tangent + 2 * out );
        double normal = 2 / 3 * value * density;
        double diagonal = value / 6 * density;
        cells[ii].state.f[ortho] = cells[ii].state.f[OPPOSITE(ortho)] + normal;
        cells[ii].state.f[ortho+1] = cells[ii].state.f[OPPOSITE(ortho+1)] + diagonal;
        cells[ii].state.f[before] = cells[ii].state.f[OPPOSITE(before)] + diagonal;
      }
      break;
    case 'o': // outlet (zero-gradient extrapolation)
      for (; ii <= imax; ii+=dii)
      {
        cells[ii].state.f[ortho] = cells[ii+djj].state.f[ortho];
        cells[ii].state.f[ortho+1] = cells[ii+djj].state.f[ortho+1];
        cells[ii].state.f[before] = cells[ii+djj].state.f[before];
      }
      break;
    case 'm': // moving wall
      for (; ii <= imax; ii+=dii)
      {
        int sign = (side == 't' or side == 'l') ? 1: -1;
        int parallel = ( ortho + OPPOSITE(ortho) ) / 2;
        double coeff = 1.0 / ( 1.0 - value );
        double tangent = cells[ii].state.fc
          + cells[ii].state.f[parallel] 
          + cells[ii].state.f[OPPOSITE(parallel)];
        double out = cells[ii].state.f[OPPOSITE(ortho)]
          + cells[ii].state.f[OPPOSITE(ortho+1)]
          + cells[ii].state.f[OPPOSITE(before)];
        double density = coeff * ( tangent + 2 * out );
        double diagonal = value / 6 * density;
        cells[ii].state.f[ortho] = cells[ii].state.f[OPPOSITE(ortho)];
        cells[ii].state.f[ortho+1] = cells[ii].state.f[OPPOSITE(ortho+1)] 
          + sign * diagonal;
        cells[ii].state.f[before] = cells[ii].state.f[OPPOSITE(before)] 
          - sign * diagonal;
      }
      break;
    default:
      break;
  }
}