#include "Cell.h"

using namespace std;

Cell::Cell(double rho, double u, double v, 
  double tau, double omega, double nu, double nuc, 
  vector<Cell>* grid_levels)
{
  state.rho = rho;
  state.u = u;
  state.v = v;
  tree.grid_levels = grid_levels;
  // Iniitalize f via equilibrium distribution function.
  double msq = state.u*state.u + state.v*state.v;
  state.fc = FEQ(4/9,state.rho,0,msq);
  for (int i = 0; i < 8; ++i) 
  {
    state.f[i] = FEQ(
      WEIGHT(i),
      state.rho,
      CX[i] * state.u + CY[i] * state.v,
      msq );
  }
}

// Currently a homogeneous copy of state.
Cell::Cell(Cell* parent)
{
  tree.parent = parent;
  // TODO: neighbour assignment.
  state = parent->state;
  tree.grid_levels = parent->tree.grid_levels;
  tree.level = parent->tree.level+1;
  numerics.nu = parent->numerics.nu / 2.0;
}

// Stands for collide, explode, stream.
void Cell::ces()
{
  if ( numerics.physical )
  {
    if ( not numerics.interface ) collide();
    if ( numerics.interface ) explode();
    stream_parallel();
    bufferize_parallel();
  }
}
// Currently just averages children.
// Need to account for cut cells, with different-volume cells.
void Cell::coalesce()
{
  if ( numerics.interface )
  {
    double fcavg = 0;
    double favg[8] = {};
    double uavg = 0;
    double vavg = 0;
    double rhoavg = 0;
    int nc = 0;
    // Sum over children.
    for(int i = 0; i < 4; ++i)
    {
      Cell* c = tree.children[i];
      if (c != nullptr)
      {
        uavg += c->state.u;
        vavg += c->state.v;
        rhoavg += c->state.rho;
        fcavg += c->state.fc;
        for(int j = 0; j < 8; ++j) favg[j] += c->state.f[j];
        ++nc;
      }
    }
    // Average.
    if (nc > 0)
    {
      state.u = uavg / nc;
      state.v = vavg / nc;
      state.rho = rhoavg / nc;
      state.fc = fcavg / nc;
      for(int j = 0; j < 8; ++j) state.f[j] = favg[j] / nc;
    }
  }
}

// For cut cells, the distribution function needs to be paid close attention to for MME conservation.
// For homogeneous explosion, no problem. 
void Cell::explode()
{
  for (int c = 0; c < 4; ++c)
  {
    Cell* ch = tree.children[c];
    if (ch != nullptr)
    {
      ch->state = state;
    }
    else
    {
      //Create cell
      Cell new_child(this);
      tree.grid_levels[tree.level+1].push_back(new_child);
      tree.children[c] = &tree.grid_levels[tree.level+1].back();
    }
  }
}

void Cell::collide()
{
  double msq = state.u*state.u + state.v*state.v;
  state.fc = numerics.omega * FEQ( 4/9, state.rho, 0, msq ) + ( 1 - numerics.omega ) * state.fc;
  for (int i = 0; i < 8; ++i)
  {
    state.f[i] = numerics.omega * FEQ(WEIGHT(i), state.rho, CX[i]*state.u + CY[i]*state.v, msq) 
      + ( 1-numerics.omega ) * state.f[i];
  }
}

void Cell::stream_parallel()
{
  for(int i = 0; i < 8; ++i)
  {
    numerics.b[i] = (tree.neighbours[OPPOSITE(i)] != nullptr) ? tree.neighbours[OPPOSITE(i)]->state.f[i] : state.f[i];
  }
}

void Cell::bufferize_parallel()
{
  for(int i = 0; i < 8; ++i) state.f[i] = numerics.b[i];
}

double Cell::get_velocity_magnitude()
{
  return sqrt(state.u*state.u + state.v*state.v);
}

void Cell::reconstruct_macro()
{
  state.rho = state.fc;
  for (int i = 0; i < 8; ++i) state.rho += state.f[i];
  state.u = 0;
  state.v = 0;
  for (int i = 0; i < 8; ++i) state.u += state.f[i]*CX[i];
  for (int i = 0; i < 8; ++i) state.v += state.f[i]*CX[i];
  state.u /= state.rho;
  state.v /= state.rho;
}


void Cell::recompute_relaxation()
{
  // After a refine operation, the lattice viscosity is updated.
  tau = 3*(nu + nuc) + 0.5;
  omega = 1.0 / tau;
}






