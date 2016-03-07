#pragma once

// A D2Q9 Lattice Boltzmann cell.

#include <cstdint>
#include <cmath>
#include <vector>

#define OPPOSITE(index) (((index)<4)?((index)+4):((index)-4))
#define WEIGHT(index) (((index)%2==0)?((1/9)):((1/36)))
#define FEQ(W,RHO,UVC,MSQ) (W*RHO*( 1 + 3*(UVC) + 4.5*(UVC)*(UVC) - 1.5*(MSQ) ))
const int CX[8] = {1,1,0,-1,-1,-1,0,1};
const int CY[8] = {0,1,1,1,0,-1,-1,-1};

// Indexing of distribution functions and neighbours correspond to the same direction. 
// Index goes from 0 to 7 in this order (starting from E, go CCW): E, NE, N, NW, W,SW, S, SE.
// Result is that even-numbered directions are othogonal to coordinate axes, and odd-numbered indices are diagonal.

class Cell
{
public:
  Cell(double u_, double v_, double rho_, std::vector<Cell>* grid_levels_); // Meant to make coarsest cells.
  Cell(Cell* parent); // meant to be called in a refine operation.
  double get_velocity_magnitude();
  void ces();
  void coalesce();
  struct
  {
    double fc = 0; // center distribution.
    double f[8] = {}; // advected distributions ( to index-correspond to the neighbours ).
    double rho = 0;
    double u = 0;
    double v = 0;
  } state;
  struct
  {
    int level = 0; // tree level. default coarse cell.
    double dim = 0; // dimension of this cell. used for adaptive-gridding.
    std::vector<Cell>* grid_levels;
    Cell* parent = nullptr;
    Cell* children[4] = { nullptr };
    Cell* neighbours[8] = { nullptr }; // cell neighbours for every lattice direction (have same level)
    uint64_t mk = (uint64_t)0; // Morton N-order key (interleave x and y, x first).
    bool refined = false; // a flag for presence of active children.
  } tree;
  struct
  {
    bool physical = true; // flag whether this cell is participating in the flow solution.
    bool interface = false; // a flag for interface status. An interface cell will not collide, only advect. Being an interface and being a 'real' cell are not mutually exclusive. Real cells must ALWAYS have neighbours (whether real or interface), but interface cells do not need (to create new) neighbours.
    bool cut = false; // true if physical surface resides in this cell.
    double tau = 0;
    double omega = 0;
    double lattice_viscosity = 0;
    double b[8] = {}; // buffers for advected distributions (for parallel advection).
  } numerics;
private:
  void collide();
  void explode();
  void stream_parallel();
  void bufferize_parallel();
};
