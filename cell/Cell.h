#pragma once

// A D2Q9 Lattice Boltzmann cell.

#include <cstdint>
#include <cmath>

#define WCENTER (4/9)
#define WORTHO (1/9)
#define WDIAG (1/36)

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
    // Buffers for each direction; more memory, but more locality. Do not have to access the tree (which is very unstructured many times).
    double b, be, bne, bn, bnw, bw, bsw, bs, bse; // Holds streamed values to facilitate parallelization.
    double get_velocity_magnitude();
    void stream();
    void collide();
    void explode();
  private:
    std::vector<Cell>* grid_levels; // array of vectors, where each vector corresponds to a refinement level.
    Cell* parent;
    Cell* children[4]; 
    int level; // level of refinement. 0: root (coarsest).
    double u,v,rho;
    uint64_t key; // Morton key for relative positioning.
    bool physical; // distinguishes cells that are physically active and that need to collide, explode, or stream (interfaces are physical).
    bool refined; // a flag for presence of active children.
    bool interface; // a flag for interface status. An interface cell will not collide, only advect. Being an interface and being a 'real' cell are mutually exclusive. Real cells must ALWAYS have neighbours (whether real or interface), but interface cells do not need (to create new) neighbours.
    double tau; // the relaxation time specific to this cell.
    double omega; // the relaxation frequency specific to this cell.
};
