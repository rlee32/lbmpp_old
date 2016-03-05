#include "Grid.h"

using namespace std;

void Grid::initialize(int cell_count_x, int cell_count_y, 
  double rho0, double u0, double v0)
{
  cell_count[0] = cell_count_x;
  cell_count[1] = cell_count_y;
  Cell default_cell(0,rho0,u0,v0);
  cells.resize(cell_count_x*cell_count_y, default_cell);
}

double Grid::get_max_velocity_magnitude()
{
  double max = cells[0].get_velocity_magnitude();
  for(vector<Cell>::iterator it = cells.begin()+1; it != cells.end(); ++it)
  {
    double test = it->get_velocity_magnitude();
    if (test > max) max = test;
  }
  return max;
}
double Grid::get_min_velocity_magnitude()
{
  double min = cells[0].get_velocity_magnitude();
  for(vector<Cell>::iterator it = cells.begin()+1; it != cells.end(); ++it)
  {
    double test = it->get_velocity_magnitude();
    if (test < min) min = test;
  }
  return min;
}