#pragma once

// A D2Q9 Lattice Boltzmann cell.

#include <cstdint>
#include <cmath>
#include <vector>

#define WCENTER (4/9)
#define WORTHO (1/9)
#define WDIAG (1/36)
#define OPPOSITE(index) (((index)<4)?((index)+4):((index)-4))
#define WEIGHT(index) (((index)%2==0)?(WORTHO):(WDIAG))

// Indexing of distribution functions and neighbours correspond to the same direction. 
// Index goes from 0 to 7 in this order (starting from E, go CCW): E, NE, N, NW, W,SW, S, SE.
// Result is that even-numbered directions are othogonal to coordinate axes, and odd-numbered indices are diagonal.

class Cell
{
public:
  Cell(int level_, double u_, double v_, double rho_);
  Cell(Cell* parent); // meant to be called in a refine operation.
  struct
  {
    double fc = 0; // center distribution.
    double f[8] = {}; // advected distributions ( to index-correspond to the neighbours ).
    double b[8] = {}; // buffers for advected distributions (for parallel advection).
    double rho = 0;
    double u = 0;
    double v = 0;
  } state;
  struct
  {
    int level; // tree level.
    double dim; // dimension of this cell. used for adaptive-gridding.
    std::vector<Cell>* grid_levels;
    Cell* parent;
    Cell* children[4];
    Cell* neighbours[8]; // cell neighbours for every lattice direction (have same level)
    uint64_t mk; // Morton N-order key (interleave x and y, x first).
    bool refined; // a flag for presence of active children.
  } tree;
  struct
  {
    bool physical = true; // flag whether this cell is participating in the flow solution.
    bool interface = false; // a flag for interface status. An interface cell will not collide, only advect. Being an interface and being a 'real' cell are mutually exclusive. Real cells must ALWAYS have neighbours (whether real or interface), but interface cells do not need (to create new) neighbours.
    bool cut = false; // true if physical surface resides in this cell.
    double tau = 0;
    double omega = 0;
    double lattice_viscosity = 0;;
  } numerics;
  double get_velocity_magnitude();
  void ces();
  void coalesce();
private:
  void copy_state(Cell* parent);
  void collide();
  void explode();
  void stream_parallel();
  void bufferize_parallel();
};
