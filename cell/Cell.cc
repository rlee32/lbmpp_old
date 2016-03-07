#include "Cell.h"


Cell::Cell(double u, double v, double rho)
{
  state.u = u;
  state.v = v;
  state.rho = rho;
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

Cell::Cell(Cell* parent_)
{
  tree.parent = parent_;
  // TODO: neighbour assignment.
  copy_state(parent_);
  tree.level = parent_->tree.level+1;
  numerics.lattice_viscosity = parent_->numerics.lattice_viscosity / 2.0;
}

// Stands for collide, explode, stream.
void Cell::ces()
{
  if ( numerics.physical )
  {
    if ( not numerics.interface ) collide();
    explode();
    stream_parallel();
    bufferize_parallel();
  }
}
void Cell::coalesce()
{
  // Call upon and average children values.
}

// Simple copy of state.
void Cell::copy_state(Cell* other)
{
  state = other->state;
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
      ch->copy_state(this);
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










