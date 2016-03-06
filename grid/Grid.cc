#include "Grid.h"

using namespace std;

void Grid::initialize(int cell_count_x, int cell_count_y, 
  double rho0, double u0, double v0)
{
  vector<Cell>& cells = grid_levels[0];
  cell_count[0] = cell_count_x;
  cell_count[1] = cell_count_y;
  Cell default_cell(0,rho0,u0,v0);
  cells.resize(cell_count_x*cell_count_y, default_cell);
  assign_neighbours();
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

void Grid::assign_neighbours()
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
      if ( not right ) cells[ii].e = &cells[ii+1];
      // Northeast
      if ( not right and not top ) cells[ii].ne = &cells[ii+1-j];
      // North
      if ( not top ) cells[ii].n = &cells[ii-j];
      // Northwest
      if ( not top and not left ) cells[ii].nw = &cells[ii-1-j];
      // West
      if ( not left ) cells[ii].w = &cells[ii-1];
      // Southwest
      if ( not left and not bottom ) cells[ii].sw = &cells[ii-1+j];
      // South
      if ( not bottom ) cells[ii].s = &cells[ii+j];
      // Southeast
      if ( not bottom and not right ) cells[ii].se = &cells[ii+j+1];
    }
  }
}

// void set_bc(char sides[4])
// {
//   for(int i = 0; i < 4; ++i) bc[i] = sides[i];
// }

// side: b (bottom), r (right), t (top), l (left)
// type: w (wall), o (outlet), i (inlet), m (moving wall)
// value: applicable to m, i
// Currently only works for coarsest level.
void Grid::enforce_bc(int side, char type, double value)
{
  vector<Cell>& cells = grid_levels[0];
  int ii = 0;
  int dii = 1;
  int imax = 1;
  switch (side)
  {
    case 'b': // bottom
      ii = 0;
      dii = 1;
      imax = cell_count[0] - 1;
      break;
    case 'r': // right
      ii = cell_count[0]-1;
      dii = cell_count[0];
      imax = cell_count[1]*cell_count[0] - 1;
      break;
    case 't': // top
      ii = (cell_count[1]-1)*cell_count[0];
      dii = 1;
      imax = cell_count[1]*cell_count[0] - 1;
      break;
    case 'l': // left
      ii = 0;
      dii = cell_count[0];
      imax = (cell_count[1]-1)*cell_count[0];
      break;
    default:
      break;
  }
  switch (type)
  {
    case 'w': // wall
      switch (side)
      {
        case 'b':
          for (; ii <= imax; ii+=dii)
          {
            cells[ii].fn = cells[ii].fs;
            cells[ii].fne = cells[ii].fsw;
            cells[ii].fnw = cells[ii].fse;
          }
          break;
        case 'r':
          for (; ii <= imax; ii+=dii)
          {
            cells[ii].fw = cells[ii].fe;
            cells[ii].fnw = cells[ii].fse;
            cells[ii].fsw = cells[ii].fne;
          }
          break;
        case 't':
          for (; ii <= imax; ii+=dii)
          {
            cells[ii].fs = cells[ii].fn;
            cells[ii].fse = cells[ii].fnw;
            cells[ii].fsw = cells[ii].fne;
          }
          break;
        case 'l':
          for (; ii <= imax; ii+=dii)
          {
            cells[ii].fe = cells[ii].fw;
            cells[ii].fse = cells[ii].fnw;
            cells[ii].fne = cells[ii].fsw;
          }
          break;
        default:
          break;
      }
    case 'i': // inlet
      switch (side)
      {
        case 'b':
          for (; ii <= imax; ii+=dii)
          {
            double coeff = 1.0 / ( 1.0 - value );
            double tangent = cells[ii].f + cells[ii].fe + cells[ii].fw;
            double out = cells[ii].fs + cells[ii].fsw + cells[ii].fse;
            double density = coeff * ( tangent + 2 * out );
            double normal = 2 / 3 * value * density;
            double diagonal = value / 6 * density;
            cells[ii].fn = cells[ii].fs + normal;
            cells[ii].fne = cells[ii].fsw + diagonal;
            cells[ii].fnw = cells[ii].fse + diagonal;
          }
          break;
        case 'r':
          for (; ii <= imax; ii+=dii)
          {
            double coeff = 1.0 / ( 1.0 - value );
            double tangent = cells[ii].f + cells[ii].fn + cells[ii].fs;
            double out = cells[ii].fe + cells[ii].fse + cells[ii].fne;
            double density = coeff * ( tangent + 2 * out );
            double normal = 2 / 3 * value * density;
            double diagonal = value / 6 * density;
            cells[ii].fw = cells[ii].fe + normal;
            cells[ii].fnw = cells[ii].fse + diagonal;
            cells[ii].fsw = cells[ii].fne + diagonal;
          }
          break;
        case 't':
          for (; ii <= imax; ii+=dii)
          {
            double coeff = 1.0 / ( 1.0 - value );
            double tangent = cells[ii].f + cells[ii].fe + cells[ii].fw;
            double out = cells[ii].fn + cells[ii].fnw + cells[ii].fne;
            double density = coeff * ( tangent + 2 * out );
            double normal = 2 / 3 * value * density;
            double diagonal = value / 6 * density;
            cells[ii].fs = cells[ii].fn + normal;
            cells[ii].fse = cells[ii].fnw + diagonal;
            cells[ii].fsw = cells[ii].fne + diagonal;
          }
          break;
        case 'l':
          for (; ii <= imax; ii+=dii)
          {
            double coeff = 1.0 / ( 1.0 - value );
            double tangent = cells[ii].f + cells[ii].fn + cells[ii].fs;
            double out = cells[ii].fw + cells[ii].fsw + cells[ii].fnw;
            double density = coeff * ( tangent + 2 * out );
            double normal = 2 / 3 * value * density;
            double diagonal = value / 6 * density;
            cells[ii].fe = cells[ii].fw + normal;
            cells[ii].fse = cells[ii].fnw + diagonal;
            cells[ii].fne = cells[ii].fsw + diagonal;
          }
          break;
        default:
          break;
      }
      break;
    case 'o': // outlet (zero-gradient extrapolation)
      switch (side)
      {
        case 'b':
          for (; ii <= imax; ii+=dii)
          {
            cells[ii].fn = cells[ii-cell_count[0]].fn;
            cells[ii].fne = cells[ii-cell_count[0]].fne;
            cells[ii].fnw = cells[ii-cell_count[0]].fnw;
          }
          break;
        case 'r':
          for (; ii <= imax; ii+=dii)
          {
            cells[ii].fw = cells[ii-1].fw;
            cells[ii].fnw = cells[ii-1].fnw;
            cells[ii].fsw = cells[ii-1].fsw;
          }
          break;
        case 't':
          for (; ii <= imax; ii+=dii)
          {
            cells[ii].fs = cells[ii+cell_count[0]].fs;
            cells[ii].fse = cells[ii+cell_count[0]].fse;
            cells[ii].fsw = cells[ii+cell_count[0]].fsw;
          }
          break;
        case 'l':
          for (; ii <= imax; ii+=dii)
          {
            cells[ii].fe = cells[ii+1].fe;
            cells[ii].fse = cells[ii+1].fse;
            cells[ii].fne = cells[ii+1].fne;
          }
          break;
        default:
          break;
      }
      break;
    case 'm': // inlet
      switch (side)
      {
        case 'b':
          for (; ii <= imax; ii+=dii)
          {
            double coeff = 1.0 / ( 1.0 - value );
            double tangent = cells[ii].f + cells[ii].fe + cells[ii].fw;
            double out = cells[ii].fs + cells[ii].fsw + cells[ii].fse;
            double density = coeff * ( tangent + 2 * out );
            double diagonal = value / 6 * density;
            cells[ii].fn = cells[ii].fs;
            cells[ii].fne = cells[ii].fsw + diagonal;
            cells[ii].fnw = cells[ii].fse - diagonal;
          }
          break;
        case 'r':
          for (; ii <= imax; ii+=dii)
          {
            double coeff = 1.0 / ( 1.0 - value );
            double tangent = cells[ii].f + cells[ii].fn + cells[ii].fs;
            double out = cells[ii].fe + cells[ii].fse + cells[ii].fne;
            double density = coeff * ( tangent + 2 * out );
            double diagonal = value / 6 * density;
            cells[ii].fw = cells[ii].fe;
            cells[ii].fnw = cells[ii].fse + diagonal;
            cells[ii].fsw = cells[ii].fne - diagonal;
          }
          break;
        case 't':
          for (; ii <= imax; ii+=dii)
          {
            double coeff = 1.0 / ( 1.0 - value );
            double tangent = cells[ii].f + cells[ii].fe + cells[ii].fw;
            double out = cells[ii].fn + cells[ii].fnw + cells[ii].fne;
            double density = coeff * ( tangent + 2 * out );
            double diagonal = value / 6 * density;
            cells[ii].fs = cells[ii].fn;
            cells[ii].fse = cells[ii].fnw + diagonal;
            cells[ii].fsw = cells[ii].fne - diagonal;
          }
          break;
        case 'l':
          for (; ii <= imax; ii+=dii)
          {
            double coeff = 1.0 / ( 1.0 - value );
            double tangent = cells[ii].f + cells[ii].fn + cells[ii].fs;
            double out = cells[ii].fw + cells[ii].fsw + cells[ii].fnw;
            double density = coeff * ( tangent + 2 * out );
            double diagonal = value / 6 * density;
            cells[ii].fe = cells[ii].fw;
            cells[ii].fse = cells[ii].fnw - diagonal;
            cells[ii].fne = cells[ii].fsw + diagonal;
          }
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
}