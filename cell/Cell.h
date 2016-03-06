#pragma once

// A D2Q9 Lattice Boltzmann cell.

#include <cstdint>
#include <cmath>

class Cell
{
  public:
    Cell(int level_, double u_, double v_, double rho_);
    double dim; // dimension of this square cell.
    // neighbour cells using north, east, south, west notation.
    Cell* e;
    Cell* ne;
    Cell* n;
    Cell* nw;
    Cell* w;
    Cell* sw;
    Cell* s;
    Cell* se;
    double f, fe, fne, fn, fnw, fw, fsw, fs, fse; // distribution function with north, east, south, west notations.
    double fbuffer; // Holds a streamed value to facilitate parallelization.
    double get_velocity_magnitude();
    void advect(int level);
    void stream();
  private:
    Cell* parent;
    Cell* children[4]; 
    int level; // level of refinement. 0: root (coarsest).
    double u,v,rho;
    uint64_t key; // Morton key for relative positioning.
    bool refined; // a flag for presence of active children.
    bool interface; // a flag for interface status. An interface cell will not collide, only advect. Being an interface and being a 'real' cell are mutually exclusive. Real cells must ALWAYS have neighbours (whether real or interface), but interface cells do not need (to create new) neighbours.
};
