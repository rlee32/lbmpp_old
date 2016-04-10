#include "Cell.h"

using namespace std;

Cell::Cell(double rho, double u, double v, 
  double tau, double omega, double nu, double nuc, 
  vector<Cell>* grid_levels)
{
  state.rho = rho;
  // cout << state.rho << endl;
  state.u = u;
  state.v = v;
  numerics.tau = tau;
  numerics.omega = omega;
  numerics.nu = nu;
  numerics.nuc = nuc;
  tree.grid_levels = grid_levels;
  // Iniitalize f via equilibrium distribution function.
  double msq = state.u*state.u + state.v*state.v;
  // cout << state.rho << endl;
  state.fc = FEQ(4.0/9.0,state.rho,0,msq);
  // cout << state.fc << endl;
  for (int i = 0; i < 8; ++i) 
  {
    state.f[i] = FEQ(
      WEIGHT(i),
      state.rho,
      CX[i] * state.u + CY[i] * state.v,
      msq );
    // cout << state.f[i] << endl;
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
  if ( numerics.physical and numerics.interface ) 
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
}

void Cell::collide(size_t relax_model)
{
  if ( numerics.physical and not numerics.interface )
  {
    double msq = state.u*state.u + state.v*state.v;
    switch(relax_model)
    {
      case 1:
        state.fc = next_fc_srt(msq);
        for (size_t i = 0; i < 8; ++i) state.f[i] = next_fi_srt(i, msq);
        break;
      case 2:

        break;
      case 3:
        {
        double m[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        compute_moment( m );
        next_f_mrt( m );
        // state.fc = next_fc_mrt( m );
        // for (size_t i = 0; i < 8; ++i) state.f[i] = next_fi_mrt( i, m );
        }
        break;
      default:
        break;
    }
  }
}

// Single-relaxation time advance of fc.
inline double Cell::next_fc_srt(double msq) const
{
  return ( numerics.omega * FEQ( 4.0/9.0, state.rho, 0, msq ) 
    + ( 1 - numerics.omega ) * state.fc );
}
// Single-relaxation time advance of f[i].
inline double Cell::next_fi_srt(size_t i, double msq) const
{
  return ( numerics.omega 
    * FEQ(WEIGHT(i), state.rho, CX[i]*state.u + CY[i]*state.v, msq) 
    + ( 1 - numerics.omega ) * state.f[i] );
}

// For MRT, computes M * x = b, for the ith entry of b.
inline double Cell::premultiply_M(size_t i) const
{
  const double* f = state.f;
  double fc = state.fc;
  switch( i )
  {
    case 0:
      return ( fc + f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] );
      break;
    case 1:
      return ( -4*fc 
        - ( f[0] + f[2] + f[4] + f[6] ) 
        + 2*( f[1] + f[3] + f[5] + f[7] ) );
      break;
    case 2:
      return ( 4*fc 
        - 2*( f[0] + f[2] + f[4] + f[6] ) 
        + ( f[1] + f[3] + f[5] + f[7] ) );
      break;
    case 3:
      return ( f[0] - f[4] 
        + f[1] - f[3] - f[5] + f[7] );
      break;
    case 4:
      return ( 2*( -f[0] + f[4] ) 
        + f[1] - f[3] - f[5] + f[7] );
      break;
    case 5:
      return ( f[1] + f[2] + f[3] 
        - f[5] - f[6] - f[7] );
      break;
    case 6:
      return ( f[1] - 2*f[2] + f[3] 
        - f[5] + 2*f[6] - f[7] );
      break;
    case 7:
      return ( f[0] - f[2] + f[4] - f[6] );
      break;
    case 8:
      return ( f[1] - f[3] + f[5] - f[7] );
      break;
    default:
      return 0;
      break;
  }
}
// For MRT, computes change in moment-space variables.
inline void Cell::compute_moment( double m[9] ) const
{
  double rhou = state.rho*state.u;
  double rhov = state.rho*state.v;
  double uu = state.u*state.u;
  double vv = state.v*state.v;
  double magmag = uu + vv;
  m[0] = premultiply_M(0) - state.rho;
  m[1] = premultiply_M(1) - state.rho * ( -2 + 3 * magmag );
  m[2] = premultiply_M(2) - state.rho * ( 1 - 3 * magmag );
  m[3] = premultiply_M(3) - rhou;
  m[4] = premultiply_M(4) + rhou;
  m[5] = premultiply_M(5) - rhov;
  m[6] = premultiply_M(6) + rhov;
  m[7] = premultiply_M(7) - state.rho * ( uu - vv );
  m[8] = premultiply_M(8) - state.rho * state.u * state.v;
}
// For MRT, computes ( (M^-1) * S ) * x = b, for the ith entry of b.
inline double Cell::premultiply_MinvS( size_t i, const double m[9] ) const
{
  switch( i )
  {
    case 0:
      return ( ( 4.0*m[0] - 4.8*m[1] + 4.0*m[2] ) / 36.0 );
      break;
    case 1:
      return ( 
        ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] + 6.0*m[3] - 7.2*m[4] 
          + 9*numerics.omega*m[7] ) 
        / 36.0 );
      break;
    case 2:
      return ( 
        ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] + 6.0*m[5] - 7.2*m[6] 
          - 9*numerics.omega*m[7] ) 
        / 36.0 );
      // return ( 
      //   ( 4.0*m[0] + 2.4*m[1] + m[2] + 6.0*m[3] + 3.6*m[4] + 6.0*m[5] + 3.6*m[6] 
      //     + 9*numerics.omega*m[8] ) 
      //   / 36.0 );
      break;
    case 3:
      return ( 
        ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] - 6.0*m[3] + 7.2*m[4] 
          + 9*numerics.omega*m[7] ) 
        / 36.0 );
      // return ( 
      //   ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] + 6.0*m[5] - 7.2*m[6] 
      //     - 9*numerics.omega*m[7] ) 
      //   / 36.0 );
      break;
    case 4:
      return ( 
        ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] - 6.0*m[5] + 7.2*m[6] 
          - 9*numerics.omega*m[7] ) 
        / 36.0 );
      // return ( 
      //   ( 4.0*m[0] + 2.4*m[1] + m[2] - 6.0*m[3] - 3.6*m[4] + 6.0*m[5] + 3.6*m[6] 
      //     - 9*numerics.omega*m[8] ) 
      //   / 36.0 );
      break;
    case 5:
      return ( 
        ( 4.0*m[0] + 2.4*m[1] + 1.0*m[2] + 6.0*m[3] + 3.6*m[4] + 6.0*m[5] + 3.6*m[6] 
          + 9*numerics.omega*m[8] ) 
        / 36.0 );
      // return ( 
      //   ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] - 6.0*m[3] + 7.2*m[4] 
      //     + 9*numerics.omega*m[7] ) 
      //   / 36.0 );
      break;
    case 6:
      return ( 
        ( 4.0*m[0] + 2.4*m[1] + 1.0*m[2] - 6.0*m[3] - 3.6*m[4] + 6.0*m[5] + 3.6*m[6] 
          - 9*numerics.omega*m[8] ) 
        / 36.0 );
      // return ( 
      //   ( 4.0*m[0] + 2.4*m[1] + m[2] - 6.0*m[3] - 3.6*m[4] - 6.0*m[5] - 3.6*m[6] 
      //     + 9*numerics.omega*m[8] ) 
      //   / 36.0 );
      break;
    case 7:
      return ( 
        ( 4.0*m[0] + 2.4*m[1] + 1.0*m[2] - 6.0*m[3] - 3.6*m[4] - 6.0*m[5] - 3.6*m[6] 
          + 9*numerics.omega*m[8] ) 
        / 36.0 );
      // return ( 
      //   ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] - 6.0*m[5] + 7.2*m[7] 
      //     - 9*numerics.omega*m[7] ) 
      //   / 36.0 );
      break;
    case 8:
      return ( 
        ( 4.0*m[0] + 2.4*m[1] + 1.0*m[2] + 6.0*m[3] + 3.6*m[4] - 6.0*m[5] - 3.6*m[6] 
          - 9*numerics.omega*m[8] ) 
        / 36.0 );
      // return ( 
      //   ( 4.0*m[0] + 2.4*m[1] + m[2] + 6.0*m[3] + 3.6*m[4] - 6.0*m[5] - 3.6*m[6] 
      //     - 9*numerics.omega*m[8] ) 
      //   / 36.0 );
      break;
    default:
      return 0;
      break;
  }
}
// MRT advance of fc.
inline double Cell::next_fc_mrt( const double m[9] ) const
{
  return state.fc - premultiply_MinvS(0, m);
}
// MRT advance of f[i].
void Cell::next_f_mrt( const double m[9] )
{
  state.fc -= premultiply_MinvS(0, m);
  state.f[0] -= premultiply_MinvS(1, m);
  state.f[2] -= premultiply_MinvS(2, m);
  state.f[4] -= premultiply_MinvS(3, m);
  state.f[6] -= premultiply_MinvS(4, m);
  state.f[1] -= premultiply_MinvS(5, m);
  state.f[3] -= premultiply_MinvS(6, m);
  state.f[5] -= premultiply_MinvS(7, m);
  state.f[7] -= premultiply_MinvS(8, m);
}
void Cell::test_mrt()
{
  double m[9] = {1,2,3,4,5,6,7,8,9};
  cout << "Testing MinvS mult: " << endl;
  for(int i = 0; i < 9; ++i) cout << premultiply_MinvS(i,m) << endl;
  state.fc = 1;
  for(int i = 0; i < 8; ++i) state.f[i] = i+2;
  cout << "Testing M mult: " << endl;
  for(int i = 0; i < 9; ++i) cout << premultiply_M(i) << endl;
}

void Cell::stream_parallel()
{
  if( numerics.physical )
  {
    for(int i = 0; i < 8; ++i)
    {
      numerics.b[i] = (tree.neighbours[OPPOSITE(i)] != nullptr) ? tree.neighbours[OPPOSITE(i)]->state.f[i] : state.f[i];
    }
  }
}

void Cell::bufferize_parallel()
{
  for(int i = 0; i < 8; ++i) state.f[i] = numerics.b[i];
}

double Cell::get_velocity_magnitude() const
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
  for (int i = 0; i < 8; ++i) state.v += state.f[i]*CY[i];
  state.u /= state.rho;
  state.v /= state.rho;
  // cout << state.rho << endl;
}


void Cell::recompute_relaxation()
{
  // After a refine operation, the lattice viscosity is updated.
  numerics.tau = 3*(numerics.nu + numerics.nuc) + 0.5;
  numerics.omega = 1.0 / numerics.tau;
}




