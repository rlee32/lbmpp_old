#include "Grid.h"

void Grid::initialize(int cell_count_x, int cell_count_y, 
  double rho0, double u0, double v0)
{
  cell_count[0] = cell_count_x;
  cell_count[1] = cell_count_y;
  Cell default_cell(0,rho0,u0,v0);
  cells.resize(cell_count_x*cell_count_y, default_cell);
}
