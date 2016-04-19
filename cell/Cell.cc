#include "Cell.h"

using namespace std;

Cell::Cell( double rho, double u, double v )
{
  reconstruct_distribution( rho, u, v );
}

// Constructor for cells generated from refinement.
// Currently a homogeneous copy of state.
Cell::Cell( Cell* parent )
{
  local.parent = parent->local.me;
  state = parent->state;
  parent->action.link_children = true;
}

void Cell::reconstruct_distribution( double rho, double u, double v )
{
  state.rho = rho;
  state.u = u;
  state.v = v;
  // Iniitalize f to equilibrium distribution function.
  double msq = state.u*state.u + state.v*state.v;
  state.fc = FEQ( 4.0/9.0, state.rho, 0, msq );
  for (int i = 0; i < 8; ++i) 
  {
    state.f[i] = FEQ(
      WEIGHT(i),
      state.rho,
      CX[i] * state.u + CY[i] * state.v,
      msq );
  }
}

// Links newly created children and reciprocates linking for those cells who 
//  are not flagged to be linked in this iteration (older cells).
// pg: parent grid, next-level / finer grid that the children reside on.
// cg: child grid, next-level / finer grid that the children reside on.
void Cell::link_children( vector<Cell>& pg, vector<Cell>& cg )
{
  // A reference to the this (parent) cell's immediate children and neighbours.
  int (&plc)[4] = local.children;
  int (&pln)[8] = local.neighbours;

  // 0th position (bottom-left)
  {
    int (&cln)[8] = cg[ plc[0] ].local.neighbours; // this child's neighbours.
    int cme = plc[0]; // this child's id.
    cln[0] = plc[2]; // Sibling
    cln[1] = plc[3]; // Sibling
    cln[2] = plc[1]; // Sibling
    int pn = pln[4]; // parent neighbour
    if ( pn >= 0 )
    {
      cln[3] = pg[ pn ].local.children[3]; // Cousin
      cln[4] = pg[ pn ].local.children[2]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[3] >= 0 ) cg[ cln[3] ].local.neighbours[7] = cme;
        if ( cln[4] >= 0 ) cg[ cln[4] ].local.neighbours[0] = cme;
      }
    }
    // Cousins
    pn = pln[5];
    if ( pn >= 0 )
    {
      cln[5] = pg[ pn ].local.children[3]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[5] >= 0 ) cg[ cln[5] ].local.neighbours[1] = cme;
      }
    }
    // Cousins
    pn = pln[6];
    if ( pn >= 0 )
    {
      cln[6] = pg[ pn ].local.children[1]; // Cousin
      cln[7] = pg[ pn ].local.children[3]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[6] >= 0 ) cg[ cln[6] ].local.neighbours[2] = cme;
        if ( cln[7] >= 0 ) cg[ cln[7] ].local.neighbours[3] = cme;
      }
    }
  }

  // 1st position (top-left)
  {
    int (&cln)[8] = cg[ plc[1] ].local.neighbours; // this child's neighbours.
    int cme = plc[1]; // this child's id.
    cln[0] = plc[3]; // sibling
    int pn = pln[2]; // parent neighbour
    if ( pn >= 0 )
    {
      cln[1] = pg[ pn ].local.children[2]; // Cousin
      cln[2] = pg[ pn ].local.children[0]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[1] >= 0 ) cg[ cln[1] ].local.neighbours[5] = cme;
        if ( cln[2] >= 0 ) cg[ cln[2] ].local.neighbours[6] = cme;
      }
    }
    pn = pln[3];
    if ( pn >= 0 )
    {
      cln[3] = pg[ pn ].local.children[2]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[3] >= 0 ) cg[ cln[3] ].local.neighbours[7] = cme;
      }
    }
    pn = pln[4];
    if ( pn >= 0 )
    {
      cln[4] = pg[ pn ].local.children[3]; // Cousin
      cln[5] = pg[ pn ].local.children[2]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[4] >= 0 ) cg[ cln[4] ].local.neighbours[0] = cme;
        if ( cln[5] >= 0 ) cg[ cln[5] ].local.neighbours[1] = cme;
      }
    }
    cln[6] = plc[0]; // sibling
    cln[7] = plc[2]; // sibling
  }

  // 2nd position (bottom-right)
  {
    int (&cln)[8] = cg[ plc[2] ].local.neighbours; // this child's neighbours.
    int cme = plc[2]; // this child's id.
    int pn = pln[0];
    if ( pn >= 0 )
    {
      cln[0] = pg[ pn ].local.children[0]; // Cousin
      cln[1] = pg[ pn ].local.children[1]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[0] >= 0 ) cg[ cln[0] ].local.neighbours[4] = cme;
        if ( cln[1] >= 0 ) cg[ cln[1] ].local.neighbours[5] = cme;
      }
    }
    cln[2] = plc[3]; // sibling
    cln[3] = plc[1]; // sibling
    cln[4] = plc[0]; // sibling
    pn = pln[6];
    if ( pn >= 0 )
    {
      cln[5] = pg[ pn ].local.children[1]; // Cousin
      cln[6] = pg[ pn ].local.children[3]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[5] >= 0 ) cg[ cln[5] ].local.neighbours[1] = cme;
        if ( cln[6] >= 0 ) cg[ cln[6] ].local.neighbours[2] = cme;
      }
    }
    pn = pln[7];
    if ( pn >= 0 )
    {
      cln[7] = pg[ pn ].local.children[1]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[7] >= 0 ) cg[ cln[7] ].local.neighbours[3] = cme;
      }
    }
  }
    
  // 3rd position (top-right)
  {
    int (&cln)[8] = cg[ plc[3] ].local.neighbours; // this child's neighbours.
    int cme = plc[3]; // this child's id.
    int pn = pln[0];
    if ( pn >= 0 )
    {
      cln[0] = pg[ pn ].local.children[1]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[0] >= 0 ) cg[ cln[0] ].local.neighbours[4] = cme;
      }
    }
    pn = pln[1];
    if ( pn >= 0 )
    {
      cln[1] = pg[ pn ].local.children[0]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[1] >= 0 ) cg[ cln[1] ].local.neighbours[5] = cme;
      }
    }
    pn = pln[2];
    if ( pn >= 0 )
    {
      cln[2] = pg[ pn ].local.children[2]; // Cousin
      cln[3] = pg[ pn ].local.children[0]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[2] >= 0 ) cg[ cln[2] ].local.neighbours[6] = cme;
        if ( cln[3] >= 0 ) cg[ cln[3] ].local.neighbours[7] = cme;
      }
    }
    cln[4] = plc[1]; // sibling
    cln[5] = plc[0]; // sibling
    cln[6] = plc[2]; // sibling
    pn = pln[0];
    if ( pn >= 0 )
    {
      cln[7] = pg[ pn ].local.children[0]; // Cousin
      if ( not pg[ pn ].action.link_children ) // reciprocate
      {
        if ( cln[7] >= 0 ) cg[ cln[7] ].local.neighbours[3] = cme;
      }
    }
  }
}

// Currently just averages children.
// Need to account for cut cells, with different-volume cells.
void Cell::coalesce()
{
  // if ( state.interface )
  // {
  //   double fcavg = 0;
  //   double favg[8] = {};
  //   double uavg = 0;
  //   double vavg = 0;
  //   double rhoavg = 0;
  //   int nc = 0;
  //   // Sum over children.
  //   for(int i = 0; i < 4; ++i)
  //   {
  //     Cell* c = tree.children[i];
  //     if (c != nullptr)
  //     {
  //       uavg += c->state.u;
  //       vavg += c->state.v;
  //       rhoavg += c->state.rho;
  //       fcavg += c->state.fc;
  //       for(int j = 0; j < 8; ++j) favg[j] += c->state.f[j];
  //       ++nc;
  //     }
  //   }
  //   // Average.
  //   if (nc > 0)
  //   {
  //     state.u = uavg / nc;
  //     state.v = vavg / nc;
  //     state.rho = rhoavg / nc;
  //     state.fc = fcavg / nc;
  //     for(int j = 0; j < 8; ++j) state.f[j] = favg[j] / nc;
  //   }
  // }
}

// For cut cells, the distribution of particles should ensure MME conservation.
// For homogeneous explosion, no problem.
// cg: child grid
void Cell::explode_homogeneous( vector<Cell>& cg )
{
  for ( size_t c = 0; c < 4; ++c )
  {
    cout << local.children[c] << endl;
    cg[ local.children[c] ].state = state;
  }
}
void Cell::activate_children( vector<Cell>& cg )
{
  explode_homogeneous( cg );
  for ( size_t i = 0; i < 4; ++i ) cg[ local.children[i] ].state.active = true;
}

// Creates activated children.
void Cell::create_children( vector<Cell>& next_level_cells )
{
  // create the children.
  for ( size_t i = 0; i < 4; ++i )
  {
    Cell child( this );
    child.local.me = next_level_cells.size();
    local.children[i] = next_level_cells.size();
    next_level_cells.push_back(child);
  }
}

// If children do not already exist, create them.
// Then, activate children and deactivate current (parent) cell.
void Cell::refine( vector<Cell>& next_level_cells )
{
  // We assume all children made or no children made.
  bool no_children = local.children[0] < 0;
  if( no_children )
  {
    create_children( next_level_cells );
  }
  else
  {
    // cout << "ACTIVATION" << endl;
    activate_children( next_level_cells );
  }
  state.active = false;
}

void Cell::collide( size_t relax_model, size_t vc_model, double omega, 
  double scale_decrease, double scale_increase, double nuc )
{
  if ( state.active and not state.interface )
  {
    if ( local.me == 0 )
    {
      cout << scale_increase << " ";
      for (int i = 0; i < 8; ++i) cout << local.neighbours[i] << " ";
      cout << endl;
      cin.ignore();
    }
    // Relaxation
    switch( relax_model )
    {
      case 1:
        {
        double msq = state.u*state.u + state.v*state.v;
        state.fc = next_fc_srt( msq, omega );
        for (size_t i=0;i < 8;++i) state.f[i] = next_fi_srt( i, msq, omega );
        }
        break;
      case 2:
        break;
      case 3:
        {
        double m[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        compute_moment( m );
        next_f_mrt( m, omega );
        }
        break;
      default:
        break;
    }
    // VC
    switch( vc_model )
    {
      case 0:
        break;
      case 1:
        apply_steady_vc_body_force( 
          omega, scale_decrease, scale_increase, nuc );
        break;
      case 2:
        break;
      case 3:
        break;
      default:
        break;
    }
  }
}

// For viscosity counteraction method.
// Compute equilibrium distribution function.
inline void Cell::compute_feq( double feq[9] ) const
{
  const double& rho = state.rho;
  const double& u = state.u;
  const double& v = state.v;
  double uu = u * u;
  double vv = v * v;
  double upv = u + v;
  double upvupv = upv * upv;
  double magmag = uu + vv;
  double mag_term = 1.5 * magmag;
  double center_term = 4.0 / 9.0 * rho;
  double orth_term = 1.0 / 9.0 * rho;
  double diag_term = 1.0 / 36.0 * rho;
  feq[0] = center_term * ( 1.0 - mag_term );
  feq[1] = orth_term * ( 1.0 + 3.0*u + 4.5*uu - mag_term );
  feq[2] = diag_term * ( 1.0 + 3.0*upv + 4.5*upvupv - mag_term );
  feq[3] = orth_term * ( 1.0 + 3.0*v + 4.5*vv - mag_term );
  feq[4] = diag_term * ( 1.0 + 3.0*(-u+v) + 4.5*(-u+v)*(-u+v) - mag_term );
  feq[5] = orth_term * ( 1.0 - 3.0*u + 4.5*uu - mag_term );
  feq[6] = diag_term * ( 1.0 - 3.0*upv + 4.5*upvupv - mag_term );
  feq[7] = orth_term * ( 1.0 - 3.0*v + 4.5*vv - mag_term );
  feq[8] = diag_term * ( 1.0 + 3.0*(u-v) + 4.5*(u-v)*(u-v) - mag_term );
}
// This computes the strain terms S11, S12, and S22, 
// but without 2*rho (it cancels in the spatial differencing). 
inline void Cell::compute_strain_terms( double& s11, double& s12, double& s22, 
  const double feq[9], double omega ) const
{
  double coeff = 3.0*omega;
  double df1 = state.f[1] - feq[2] + state.f[5] - feq[6];
  double df2 = state.f[3] - feq[4] + state.f[7] - feq[8];
  double df3 = df1 + df2;
  s11 = coeff * ( df1 - df2 );
  s12 = coeff * ( df3 + state.f[0] - feq[1] + state.f[4] - feq[5] );
  s22 = coeff * ( df3 + state.f[2] - feq[3] + state.f[6] - feq[7] );
}
inline void Cell::fill_strain_terms( double omega )
{
  double feq[9] = {0,0,0,0,0,0,0,0,0};
  compute_feq(feq);
  compute_strain_terms(vc.s11, vc.s12, vc.s22, feq, omega);
}
// differencing.
// we assume for each x and y direction, 
//  there can only be a one-way deficiency of neighbours (if any).
// assumes at least 5 cells present in domain in each x and y direction.
// dh_inv: inverse of grid spacing.
inline void Cell::compute_strain_differences( 
  double& s11x, double& s12x, double& s12y, double& s22y, double dh_inv ) const
{
  // switch ( tree.fully_interior_cell )
  // {
  //   case true: // 4th order, 4-point stencil
  //     {
  //     const Cell* const& e1 = tree.neighbours[0];
  //     const Cell* const& e2 = tree.neighbours[0]->tree.neighbours[0];
  //     const Cell* const& w1 = tree.neighbours[4];
  //     const Cell* const& w2 = tree.neighbours[4]->tree.neighbours[4];
  //     s11x = dh_inv * ( 
  //       -e2->vc.s11 + 8.0*e1->vc.s11 - 8.0*w1->vc.s11 + w2->vc.s11 ) / 12.0;
  //     s12x = dh_inv * ( 
  //       -e2->vc.s12 + 8.0*e1->vc.s12 - 8.0*w1->vc.s12 + w2->vc.s12 ) / 12.0;
  //     const Cell* const& n1 = tree.neighbours[2];
  //     const Cell* const& n2 = tree.neighbours[2]->tree.neighbours[2];
  //     const Cell* const& s1 = tree.neighbours[6];
  //     const Cell* const& s2 = tree.neighbours[6]->tree.neighbours[6];
  //     s12y = dh_inv * ( 
  //       -n2->vc.s12 + 8.0*n1->vc.s12 - 8.0*s1->vc.s12 + s2->vc.s12 ) / 12.0;
  //     s22y = dh_inv * ( 
  //       -n2->vc.s22 + 8.0*n1->vc.s22 - 8.0*s1->vc.s22 + s2->vc.s22 ) / 12.0;
  //     }
  //     break;
  //   case false:
  //     {
  //     // x-direction
  //     switch( tree.nn[0] )
  //     {
  //       case 0: // backward 2nd order, 3-point stencil
  //         {
  //         const Cell* const& w1 = tree.neighbours[4];
  //         const Cell* const& w2 = tree.neighbours[4]->tree.neighbours[4];
  //         s11x = dh_inv * -( 
  //           -3.0*vc.s11 + 4.0*w1->vc.s11 - w2->vc.s11 ) / 2.0;
  //         s12x = dh_inv * -( 
  //           -3.0*vc.s12 + 4.0*w1->vc.s12 - w2->vc.s12 ) / 2.0;
  //         }
  //         break;
  //       case 1: // backward 4th order, 5-point stencil
  //         {
  //         const Cell* const& e1 = tree.neighbours[0];
  //         const Cell* const& w1 = tree.neighbours[4];
  //         const Cell* const& w2 = tree.neighbours[4]->tree.neighbours[4];
  //         const Cell* const& w3 = 
  //           tree.neighbours[4]->tree.neighbours[4]->tree.neighbours[4];
  //         s11x = dh_inv * -( 
  //           - 3.0*e1->vc.s11 - 10.0*vc.s11 + 18.0*w1->vc.s11 
  //           - 6.0*w2->vc.s11 + w3->vc.s11 ) / 12.0;
  //         s12x = dh_inv * -( 
  //           - 3.0*e1->vc.s12 - 10.0*vc.s12 + 18.0*w1->vc.s12 
  //           - 6.0*w2->vc.s12 + w3->vc.s12 ) / 12.0;
  //         }
  //         break;
  //       default: // now check other direction. 
  //         switch( tree.nn[2] )
  //         {
  //           case 0: // forward 2nd order, 3-point stencil
  //             {
  //             const Cell* const& e1 = tree.neighbours[0];
  //             const Cell* const& e2 = tree.neighbours[0]->tree.neighbours[0];
  //             s11x = dh_inv * (
  //               -3.0*vc.s11 + 4.0*e1->vc.s11 - e2->vc.s11 ) / 2.0;
  //             s12x = dh_inv * ( 
  //               -3.0*vc.s12 + 4.0*e1->vc.s12 - e2->vc.s12 ) / 2.0;
  //             }
  //             break;
  //           case 1: // forward 4th order, 5-point stencil
  //             {
  //             const Cell* const& w1 = tree.neighbours[4];
  //             const Cell* const& e1 = tree.neighbours[0];
  //             const Cell* const& e2 = tree.neighbours[0]->tree.neighbours[0];
  //             const Cell* const& e3 = 
  //               tree.neighbours[0]->tree.neighbours[0]->tree.neighbours[0];
  //             s11x = dh_inv * ( 
  //               - 3.0*w1->vc.s11 - 10.0*vc.s11 + 18.0*e1->vc.s11 
  //               - 6.0*e2->vc.s11 + e3->vc.s11 ) / 12.0;
  //             s12x = dh_inv * ( 
  //               - 3.0*w1->vc.s12 - 10.0*vc.s12 + 18.0*e1->vc.s12 
  //               - 6.0*e2->vc.s12 + e3->vc.s12 ) / 12.0;
  //             }
  //             break;
  //           default: // centered 4th order, 4-point stencil
  //             {
  //             const Cell* const& e1 = tree.neighbours[0];
  //             const Cell* const& e2 = tree.neighbours[0]->tree.neighbours[0];
  //             const Cell* const& w1 = tree.neighbours[4];
  //             const Cell* const& w2 = tree.neighbours[4]->tree.neighbours[4];
  //             s11x = dh_inv * ( 
  //               -e2->vc.s11 + 8.0*e1->vc.s11 - 8.0*w1->vc.s11 + w2->vc.s11 ) / 12.0;
  //             s12x = dh_inv * ( 
  //               -e2->vc.s12 + 8.0*e1->vc.s12 - 8.0*w1->vc.s12 + w2->vc.s12 ) / 12.0;
  //             }
  //             break;
  //         }
  //       break;
  //     }
  //     // y-direction
  //     switch( tree.nn[1] )
  //     {
  //       case 0: // backward 2nd order, 3-point stencil
  //         {
  //         const Cell* const& s1 = tree.neighbours[6];
  //         const Cell* const& s2 = tree.neighbours[6]->tree.neighbours[6];
  //         s11x = dh_inv * -(
  //           -3.0*vc.s11 + 4.0*s1->vc.s11 - s2->vc.s11 ) / 2.0;
  //         s12x = dh_inv * -(
  //           -3.0*vc.s12 + 4.0*s1->vc.s12 - s2->vc.s12 ) / 2.0;
  //         }
  //         break;
  //       case 1: // backward 4th order, 5-point stencil
  //         {
  //         const Cell* const& n1 = tree.neighbours[2];
  //         const Cell* const& s1 = tree.neighbours[6];
  //         const Cell* const& s2 = tree.neighbours[6]->tree.neighbours[6];
  //         const Cell* const& s3 = 
  //           tree.neighbours[6]->tree.neighbours[6]->tree.neighbours[6];
  //         s11x = dh_inv * -( 
  //           - 3.0*n1->vc.s11 - 10.0*vc.s11 + 18.0*s1->vc.s11 
  //           - 6.0*s2->vc.s11 + s3->vc.s11 ) / 12.0;
  //         s12x = dh_inv * -( 
  //           - 3.0*n1->vc.s12 - 10.0*vc.s12 + 18.0*s1->vc.s12 
  //           - 6.0*s2->vc.s12 + s3->vc.s12 ) / 12.0;
  //         }
  //         break;
  //       default: // now check other direction. 
  //         switch( tree.nn[3] )
  //         {
  //           case 0: // forward 2nd order, 3-point stencil
  //             {
  //             const Cell* const& n1 = tree.neighbours[2];
  //             const Cell* const& n2 = tree.neighbours[2]->tree.neighbours[2];
  //             s11x = dh_inv * (
  //               -3.0*vc.s11 + 4.0*n1->vc.s11 - n2->vc.s11 ) / 2.0;
  //             s12x = dh_inv * ( 
  //               -3.0*vc.s12 + 4.0*n1->vc.s12 - n2->vc.s12 ) / 2.0;
  //             }
  //             break;
  //           case 1: // forward 4th order, 5-point stencil
  //             {
  //             const Cell* const& s1 = tree.neighbours[6];
  //             const Cell* const& n1 = tree.neighbours[2];
  //             const Cell* const& n2 = tree.neighbours[2]->tree.neighbours[2];
  //             const Cell* const& n3 = 
  //               tree.neighbours[2]->tree.neighbours[2]->tree.neighbours[2];
  //             s11x = dh_inv * ( 
  //               - 3.0*s1->vc.s11 - 10.0*vc.s11 + 18.0*n1->vc.s11 
  //               - 6.0*n2->vc.s11 + n3->vc.s11 ) / 12.0;
  //             s12x = dh_inv * ( 
  //               - 3.0*s1->vc.s12 - 10.0*vc.s12 + 18.0*n1->vc.s12 
  //               - 6.0*n2->vc.s12 + n3->vc.s12 ) / 12.0;
  //             }
  //             break;
  //           default: // centered 4th order, 4-point stencil
  //             {
  //             const Cell* const& n1 = tree.neighbours[2];
  //             const Cell* const& n2 = tree.neighbours[2]->tree.neighbours[2];
  //             const Cell* const& s1 = tree.neighbours[6];
  //             const Cell* const& s2 = tree.neighbours[6]->tree.neighbours[6];
  //             s12y = dh_inv * ( 
  //               -n2->vc.s12 + 8.0*n1->vc.s12 - 8.0*s1->vc.s12 + s2->vc.s12 ) / 12.0;
  //             s22y = dh_inv * ( 
  //               -n2->vc.s22 + 8.0*n1->vc.s22 - 8.0*s1->vc.s22 + s2->vc.s22 ) / 12.0;
  //             }
  //             break;
  //         }
  //       break;
  //     }
  //     }
  //     break;
  //   default:
  //     break;
  // }
}
inline void Cell::compute_vc_body_force( double g[9], double nuc ) const
{
  double Fx = -nuc * ( vc.s11x + vc.s12x );
  double Fy = -nuc * ( vc.s22y + vc.s12y );
  const double& u = state.u;
  const double& v = state.v;
  g[0] = 4.0 / 3.0 * ( Fx*(-u) + Fy*(-v) );
  g[1] = 1.0 / 3.0 * ( Fx*(1+2*u) + Fy*(-v) );
  g[2] = 1.0 / 12.0 * ( Fx*(1+2*u+3*v) + Fy*(1+3*u+2*v) );
  g[3] = 1.0 / 3.0 * ( Fx*(-u) + Fy*(1+2*v) );
  g[4] = 1.0 / 12.0 * ( Fx*(-1+2*u-3*v) + Fy*(1-3*u+2*v) );
  g[5] = 1.0 / 3.0 * ( Fx*(-1+2*u) + Fy*(-v) );
  g[6] = 1.0 / 12.0 * ( Fx*(-1+2*u+3*v) + Fy*(-1+3*u+2*v) );
  g[7] = 1.0 / 3.0 * ( Fx*(-u) + Fy*(-1+2*v) );
  g[8] = 1.0 / 12.0 * ( Fx*(1+2*u-3*v) + Fy*(-1-3*u+2*v) );
}
inline void Cell::apply_steady_vc_body_force( double omega, double dt, 
  double dh_inv, double nuc )
{
  double g[9] = {0,0,0,0,0,0,0,0,0};
  fill_strain_terms( omega );
  compute_strain_differences( vc.s11x, vc.s12x, vc.s12y, vc.s22y, dh_inv );
  compute_vc_body_force( g, nuc );
  state.fc += dt*g[0];
  for(size_t i = 0; i < 8; ++i) state.f[i] += dt*g[i+1];
}

// Single-relaxation time advance of fc.
inline double Cell::next_fc_srt( double msq, double omega ) const
{
  return ( omega * FEQ( 4.0/9.0, state.rho, 0, msq ) 
    + ( 1 - omega ) * state.fc );
}
// Single-relaxation time advance of f[i].
inline double Cell::next_fi_srt( size_t i, double msq, double omega ) const
{
  return ( omega 
    * FEQ(WEIGHT(i), state.rho, CX[i]*state.u + CY[i]*state.v, msq) 
    + ( 1 - omega ) * state.f[i] );
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
inline double Cell::premultiply_MinvS( size_t i, 
  const double m[9], double omega ) const
{
  switch( i )
  {
    case 0:
      return ( ( 4.0*m[0] - 4.8*m[1] + 4.0*m[2] ) / 36.0 );
      break;
    case 1:
      return ( 
        ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] + 6.0*m[3] - 7.2*m[4] 
          + 9*omega*m[7] ) 
        / 36.0 );
      break;
    case 2:
      return ( 
        ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] + 6.0*m[5] - 7.2*m[6] 
          - 9*omega*m[7] ) 
        / 36.0 );
      break;
    case 3:
      return ( 
        ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] - 6.0*m[3] + 7.2*m[4] 
          + 9*omega*m[7] ) 
        / 36.0 );
      break;
    case 4:
      return ( 
        ( 4.0*m[0] - 1.2*m[1] - 2.0*m[2] - 6.0*m[5] + 7.2*m[6] 
          - 9*omega*m[7] ) 
        / 36.0 );
      break;
    case 5:
      return ( 
        ( 4.0*m[0] + 2.4*m[1] + 1.0*m[2] + 6.0*m[3] + 3.6*m[4] + 6.0*m[5] + 3.6*m[6] 
          + 9*omega*m[8] ) 
        / 36.0 );
      break;
    case 6:
      return ( 
        ( 4.0*m[0] + 2.4*m[1] + 1.0*m[2] - 6.0*m[3] - 3.6*m[4] + 6.0*m[5] + 3.6*m[6] 
          - 9*omega*m[8] ) 
        / 36.0 );
      break;
    case 7:
      return ( 
        ( 4.0*m[0] + 2.4*m[1] + 1.0*m[2] - 6.0*m[3] - 3.6*m[4] - 6.0*m[5] - 3.6*m[6] 
          + 9*omega*m[8] ) 
        / 36.0 );
      break;
    case 8:
      return ( 
        ( 4.0*m[0] + 2.4*m[1] + 1.0*m[2] + 6.0*m[3] + 3.6*m[4] - 6.0*m[5] - 3.6*m[6] 
          - 9*omega*m[8] ) 
        / 36.0 );
      break;
    default:
      return 0;
      break;
  }
}
// MRT advance of f[i].
void Cell::next_f_mrt( const double m[9], double omega )
{
  state.fc -= premultiply_MinvS(0, m, omega );
  state.f[0] -= premultiply_MinvS(1, m, omega );
  state.f[2] -= premultiply_MinvS(2, m, omega );
  state.f[4] -= premultiply_MinvS(3, m, omega );
  state.f[6] -= premultiply_MinvS(4, m, omega );
  state.f[1] -= premultiply_MinvS(5, m, omega );
  state.f[3] -= premultiply_MinvS(6, m, omega );
  state.f[5] -= premultiply_MinvS(7, m, omega );
  state.f[7] -= premultiply_MinvS(8, m, omega );
}

// g: grid level of this cell
void Cell::stream_parallel( vector<Cell>& g )
{
  if( state.active )
  {
    for(size_t i = 0; i < 8; ++i)
    {
      state.b[i] = ( local.neighbours[OPPOSITE(i)] > -1 ) ? 
        g[ local.neighbours[OPPOSITE(i)] ].state.f[i] : state.b[i];
    }
  }
}

void Cell::bufferize_parallel()
{
  for(size_t i = 0; i < 8; ++i) state.f[i] = state.b[i];
}

double Cell::get_mag() const
{
  return sqrt(state.u*state.u + state.v*state.v);
}

void Cell::reconstruct_macro()
{
  if ( state.active )
  {
    state.rho = state.fc;
    for (int i = 0; i < 8; ++i) state.rho += state.f[i];
    state.u = 0;
    state.v = 0;
    for (int i = 0; i < 8; ++i) state.u += state.f[i]*CX[i];
    for (int i = 0; i < 8; ++i) state.v += state.f[i]*CY[i];
    state.u /= state.rho;
    state.v /= state.rho;
  }
}

// Stationary wall boundary condition.
// Actually stores the bounced distribution in the buffer.
void Cell::bounce_back(char side)
{
  // cout << "applying bounce back" << endl;
  switch(side)
  {
    case 't': // top
      state.b[7] = state.f[3];
      state.b[6] = state.f[2];
      state.b[5] = state.f[1];
      break;
    case 'l': // left
      state.b[1] = state.f[5];
      state.b[0] = state.f[4];
      state.b[7] = state.f[3];
      break;
    case 'r': // right
      state.b[3] = state.f[7];
      state.b[4] = state.f[0];
      state.b[5] = state.f[1];
      break;
    case 'b': // bottom
      state.b[1] = state.f[5];
      state.b[2] = state.f[6];
      state.b[3] = state.f[7];
      break;
    default:
      break;
  }
}

// Moving wall bc.
// Actually stores the bounced distribution in the buffer.
void Cell::moving_wall(char side, double U)
{
  double incident = 0;
  double rho = 0;
  double A = 0;
  // cout << "applying " << U << " moving wall to " << side << endl;
  switch(side)
  {
    case 't': // top
      incident = state.f[1] + state.f[2] + state.f[3];
      rho = state.fc + state.f[0] + state.f[4] + 2*incident;
      A = ( rho * U ) / 6.0;
      state.b[7] = state.f[3] + A;
      state.b[6] = state.f[2];
      state.b[5] = state.f[1] - A;
      break;
    case 'l': // left
      incident = state.f[3] + state.f[4] + state.f[5];
      rho = state.fc + state.f[2] + state.f[6] + 2*incident;
      A = ( rho * U ) / 6.0;
      state.b[1] = state.f[5] + A;
      state.b[0] = state.f[4];
      state.b[7] = state.f[3] - A;
      break;
    case 'r': // right
      incident = state.f[1] + state.f[0] + state.f[7];
      rho = state.fc + state.f[2] + state.f[6] + 2*incident;
      A = ( rho * U ) / 6.0;
      state.b[3] = state.f[7] + A;
      state.b[4] = state.f[0];
      state.b[5] = state.f[1] - A;
      break;
    case 'b': // bottom
      incident = state.f[5] + state.f[6] + state.f[7];
      rho = state.fc + state.f[0] + state.f[4] + 2*incident;
      A = ( rho * U ) / 6.0;
      state.b[1] = state.f[5] + A;
      state.b[2] = state.f[6];
      state.b[3] = state.f[7] - A;
      break;
    default:
      break;
  }
}
