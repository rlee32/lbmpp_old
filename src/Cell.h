#pragma once

// A D2Q9 Lattice Boltzmann cell.

class Cell
{
  public:
    Cell(int level_, double f_, double u_, double v_, double rho_);
  private:
    Cell* parent;
    Cell* children[4];
    Cell* neighbours[8];
    int level; // level of refinement. 0: root (coarsest).
    double f[9];
    double u,v,rho;
};
