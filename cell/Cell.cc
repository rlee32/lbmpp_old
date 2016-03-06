#include "Cell.h"


Cell::Cell(int level_, double u_, double v_, double rho_) : 
  level(level_), u(u_), v(v_), rho(rho_)
{
  parent = nullptr;
  for(int i = 0; i < 4; ++i) children[i] = nullptr;
  e = nullptr;
  ne = nullptr;
  n = nullptr;
  nw = nullptr;
  w = nullptr;
  sw = nullptr;
  s = nullptr;
  se = nullptr;
  f = 0;
  fe = 0;
  fne = 0;
  fn = 0;
  fnw = 0;
  fw = 0;
  fsw = 0;
  fs = 0;
  fse = 0;
  b=0; be=0; bne=0; bn=0; bnw=0; bw=0; bsw=0; bs=0; bse=0;
  // Iniitalize f via equilibrium distribution function.
}

Cell::Cell(Cell* parent_)
{
  parent = parent_;
  for(int i = 0; i < 4; ++i) children[i] = nullptr;
  // TODO: neighbour assignment.
  e = nullptr;
  ne = nullptr;
  n = nullptr;
  nw = nullptr;
  w = nullptr;
  sw = nullptr;
  s = nullptr;
  se = nullptr;
  copy_state(parent_);
  b=0; be=0; bne=0; bn=0; bnw=0; bw=0; bsw=0; bs=0; bse=0;
  level = parent_->level+1;
  lattice_viscosity = parent_->lattice_viscosity / 2;
}

// void Cell::recurse()
// {
//   if ( physical )
//   {
//     if ( not interface ) collide();
//     explode();
//     stream_all();
//   }
//   if ( children[0] != nullptr ) children[0]->recurse();
//   if ( children[1] != nullptr ) children[1]->recurse();
//   if ( children[2] != nullptr ) children[2]->recurse();
//   if ( children[3] != nullptr ) children[3]->recurse();
// }

// Stands for collide, explode, stream.
void Cell::ces()
{
  if ( physical )
  {
    if ( not interface ) collide();
    explode();
    stream_parallel();
    bufferize_parallel();
  }
}
void Cell::coalesce()
{
}

// Simple copy of state.
void Cell::copy_state(Cell* other)
{
  f = other->f;
  fe = other->fe;
  fn = other->fn;
  fs = other->fs;
  fw = other->fw;
  fne = other->fne;
  fnw = other->fnw;
  fsw = other->fsw;
  fse = other->fse;
  rho = other->rho;
  u = other->u;
  v = other->v;
}

// For cut cells, the distribution function needs to be paid close attention to for MME conservation.
// For homogeneous explosion, no problem. 
void Cell::explode()
{
  for (int c = 0; c < 4; ++c)
  {
    Cell* ch = children[c];
    if (ch != nullptr)
    {
      ch->copy_state(this);
    }
    else
    {
      //Create cell
      Cell new_child(this);
      grid_levels[level+1].push_back(new_child);
      children[c] = &grid_levels[level+1].back();
    }
  }
}

void Cell::collide()
{
  double usq = u*u;
  double vsq = v*v;
  double uv = u*v;
  double msq = usq + vsq;
  double feq;
  // Center
  feq = WCENTER*rho*( 1 - 1.5*msq );
  f = omega*feq + (1-omega)*f;
  // East
  feq = WORTHO*rho*( 1 + 3*u + 4.5*usq - 1.5*msq );
  fe = omega*feq + (1-omega)*fe;
  // North
  feq = WORTHO*rho*( 1 + 3*v + 4.5*vsq - 1.5*msq );
  fn = omega*feq + (1-omega)*fn;
  // West
  feq = WORTHO*rho*( 1 + 3*-u + 4.5*usq - 1.5*msq );
  fw = omega*feq + (1-omega)*fw;
  // South
  feq = WORTHO*rho*( 1 + 3*-v + 4.5*vsq - 1.5*msq );
  fs = omega*feq + (1-omega)*fs;
  // Northeast
  feq = WDIAG*rho*( 1 + 3*(u+v) + 4.5*(msq + 2*uv) - 1.5*msq );
  fne = omega*feq + (1-omega)*fne;
  // Northwest
  feq = WDIAG*rho*( 1 + 3*(-u+v) + 4.5*(msq - 2*uv) - 1.5*msq );
  fnw = omega*feq + (1-omega)*fnw;
  // Southwest
  feq = WDIAG*rho*( 1 + 3*(-u-v) + 4.5*(msq + 2*uv) - 1.5*msq );
  fsw = omega*feq + (1-omega)*fsw;
  // Southeast
  feq = WDIAG*rho*( 1 + 3*(u-v) + 4.5*(msq - 2*uv) - 1.5*msq );
  fse = omega*feq + (1-omega)*fse;
}

void Cell::stream_parallel()
{
  be = (w != nullptr) ? w->fe : 0;
  bn = (s != nullptr) ? s->fn : 0;
  bw = (e != nullptr) ? e->fw : 0;
  bs = (n != nullptr) ? n->fs : 0;
  bne = (sw != nullptr) ? sw->fne : 0;
  bnw = (se != nullptr) ? se->fnw : 0;
  bsw = (ne != nullptr) ? ne->fsw : 0;
  bse = (nw != nullptr) ? nw->fse : 0;
}

void Cell::bufferize_parallel()
{
  fe = be;
  fn = bn;
  fw = bw;
  fs = bs;
  fne = bne;
  fnw = bnw;
  fsw = bsw;
  fse = bse;
}

double Cell::get_velocity_magnitude()
{
  return sqrt(u*u + v*v);
}










