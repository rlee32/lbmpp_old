#pragma once

// A D2Q9 Lattice Boltzmann cell.

#include <cstdint>
#include <cmath>
#include <vector>

#define WCENTER (4/9)
#define WORTHO (1/9)
#define WDIAG (1/36)

class Cell
{
public:
  Cell(int level_, double u_, double v_, double rho_);
  Cell(Cell* parent); // meant to be called in a refine operation.
  typedef struct State
  {
    double fc; // center distribution.
    double f[8]; // advected distributions ( to index-correspond to the neighbours ).
    double rho,u,v;
  } State;
  typedef struct TreeData
  {
    int level; // tree level.
    double dim; // dimension of this cell.
    std::vector<Cell>* grid_levels;
    Cell* parent;
    Cell* children[4];
    Cell* neighbours[8];
    uint64_t mk; // Morton N-order key (interleave x and y, x first).
    bool refined; // a flag for presence of active children.
  } TreeData;
  typedef struct PhysicalData
  {
    bool physical;
    bool interface; // a flag for interface status. An interface cell will not collide, only advect. Being an interface and being a 'real' cell are mutually exclusive. Real cells must ALWAYS have neighbours (whether real or interface), but interface cells do not need (to create new) neighbours.
    bool cut; // true if physical surface resides in this cell.
    double tau, omega, lattice_viscosity;
  } PhysicalData; 
  double dim; // dimension of this square cell.
  int level; // level of refinement. 0: root (coarsest).
  double u,v,rho;
  double f, fe, fne, fn, fnw, fw, fsw, fs, fse; // distribution function with north, east, south, west notations.
  // Buffers for each direction; more memory, but more locality. Do not have to access the tree (which is very unstructured many times).
  double b, be, bne, bn, bnw, bw, bsw, bs, bse; // Holds streamed values to facilitate parallelization.
  // neighbour cells using north, east, south, west notation.
  Cell* e;
  Cell* ne;
  Cell* n;
  Cell* nw;
  Cell* w;
  Cell* sw;
  Cell* s;
  Cell* se; 
  double lattice_viscosity;
  double get_velocity_magnitude();
  void copy_state(Cell* parent);
  void ces();
  void coalesce();
private:
  std::vector<Cell>* grid_levels; // array of vectors, where each vector corresponds to a refinement level.
  Cell* parent;
  Cell* children[4];// Morton N-order (interleave x and y bits, x first.)
  uint64_t key; // Morton key for relative positioning.
  bool physical; // distinguishes cells that are physically active and that need to collide, explode, or stream (interfaces are physical).
  bool refined; 
  bool interface; // a flag for interface status. An interface cell will not collide, only advect. Being an interface and being a 'real' cell are mutually exclusive. Real cells must ALWAYS have neighbours (whether real or interface), but interface cells do not need (to create new) neighbours.
  bool cut; // true if physical surface resides in this cell.
  double tau; // the relaxation time specific to this cell.
  double omega; // the relaxation frequency specific to this cell.
  void collide();
  void explode();
  void stream_parallel();
  void bufferize_parallel();
};
