#pragma once

// A D2Q9 Lattice Boltzmann cell.

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>

#define OPPOSITE(index) (((index)<4)?((index)+4):((index)-4))
#define WEIGHT(index) (((index)%2==0)?((1.0/9.0)):((1.0/36.0)))
#define FEQ(W,RHO,UVC,MSQ) ((W)*(RHO)*( 1 + 3.0*(UVC) + 4.5*(UVC)*(UVC) - 1.5*(MSQ) ))
const static int CX[8] = {1,1,0,-1,-1,-1,0,1};
const static int CY[8] = {0,1,1,1,0,-1,-1,-1};
const static double W_orth = 1.0/9.0;
const static double W_diag = 1.0/36.0;
const static double W[8] = { 
  W_orth, W_diag, 
  W_orth, W_diag, 
  W_orth, W_diag, 
  W_orth, W_diag
};

// Indexing of distribution functions and neighbours are the same. 
// Index goes from 0 to 7 in this order (starting from E, go CCW): 
//  E, NE, N, NW, W, SW, S, SE.
// Result is that even-numbered directions are othogonal to coordinate axes, 
//  and odd-numbered indices are diagonal.

class Cell
{
public:
  // Constructors 
  // Cell( double u_, double v_, double rho_ );// Meant for coarsest cells.
  Cell( double u_, double v_, double rho_, 
    std::vector<Cell>* g, std::vector<Cell>* cg );// Meant for coarsest cells.
  // Cell( Cell* parent ); // meant for (single-level) refine operation.
  Cell( Cell* parent, 
    std::vector<Cell>* cg ); // meant for (single-level) refine operation.
  
  // state
  double u() const { return state.u; }
  double v() const { return state.v; }
  bool active() const { return state.active; }

  // Iteration
  void collide( std::size_t relax_model, std::size_t vc_model, 
    double omega, double scale_decrease, double scale_increase, double nuc );
  void stream_parallel( std::vector<Cell>& g );
  void bufferize_parallel();
  void reconstruct_macro();

  // BCs
  void bounce_back(char side);
  void moving_wall(char side, double U);
  
  // For dynamic mesh.
  void explode();
  void coalesce();
  void refine( std::vector<Cell>& next_level_cells, 
    std::vector<Cell>& grandchild_cells );
  void create_interface_children( std::vector<Cell>& next_level_cells, 
    std::vector<Cell>& grandchild_cells );
  // Local connectivity.
  bool has_neighbour(std::size_t i) const { return local.neighbours[i] > -1; }
  bool has_children() const { return local.children[0] > -1; }
  bool has_interface_children( std::vector<Cell>& next_level_cells ) const 
    { return next_level_cells[ local.children[0] ].state.interface; }
  Cell& parent() { return (*(local.pg))[local.parent]; }
  
  // VC
  void stream_body_force_parallel();
  void bufferize_body_force_parallel();
  void apply_advected_vc_body_force( double omega, double dt, 
    double dh_inv, double nuc );

  // For initialization from file.
  void set_uv( double u, double v ) { state.u = u; state.v = v; };
  void reconstruct_distribution( double u_, double v_, double rho_ );
  
  // For post-processing.
  double get_mag() const;
  double rho() const { return state.rho; }
  void link_children( std::vector<Cell>& pg, std::vector<Cell>& cg );
  // std::size_t max_active_level(std::vector<Cell>& next_level_cells)

  // Operators
  // Neighbour accessor
  Cell& operator[](std::size_t i) const { return (*(local.g))[local.neighbours[i]]; }
  // Child accessor
  Cell& operator()(std::size_t i) const { return (*(local.cg))[local.children[i]]; }


  struct
  {
    double fc; // center distribution.
    // advected distributions ( to index-correspond to the neighbours ).
    double f[8]; 
    double rho;
    double u;
    double v;
    // buffers for advected distributions (for parallel advection).
    // also useful for holding bounced-back distributions!
    double b[8];
    // 
    bool active;
    // An interface cell will not collide, only advect. 
    // A cell can be both an interface and active if it is on the coarser 
    //  layer of an interface.
    // Non-boundary real cells must ALWAYS have neighbours 
    //  (whether real or interface) but interface cells do not need neighbours.
    bool interface;
  } state;
  struct
  {
    // < 0 indices indicate invalid values.
    // The vector index for this cell in its grid level.
    int me;
    // Since Cells at each level are stored in vectors, we cannot use pointers 
    //  because vectors may resize (which invalidates pointers). So, we know 
    //  that parents are on the immediate level above and children on the 
    //  immediate level below. Therefore, we can simply store vector indices!
    int parent;
    // Children are ordered in N-order curve.
    int children[4];
    // cell neighbours for every lattice direction (have same level)
    // the number of neighbours in each orthogonal direction (at least).
    // going ccw from {1,0}.
    int neighbours[8];
    // grid levels
    std::vector<Cell>* pg; // parent grid.
    std::vector<Cell>* g; // the current grid level.
    std::vector<Cell>* cg; // child grid level.
    // Currently used for VC
    std::size_t nn[4];
    // following meaning there are at least 2 neighbours in every direction.
    bool fully_interior_cell;
  } local;
  // These variables only last during one entire-grid iteration 
  //  (persists through sub-grid iterations).
  // These need to be explicitly set during every whole-grid iteration.
  struct
  {
    // A flag to identify cells that need refinement at the end of the current 
    //  whole-grid iteration.
    bool refine;
    // Link newly-created children (not this cell itself).
    bool link_children;
  } action;
  // For viscosity-counteracting approach.
  // These are merely buffers for differencing.
  struct
  {
    // Strain rate tensor values 
    // (with 2*rho omitted, since it cancels in differencing).
    double s11; 
    double s12;
    double s22;
    // spatial differences.
    double s11x;
    double s12x;
    double s12y;
    double s22y;
    // body forces
    double gc; // center body force term
    double g[8]; // directional body force terms
    double b[8]; // buffer for transported g
  } vc;
  // Okay, this adds some memory. 
  // But hey, it makes things a lot easier!!!
  struct
  {
    int corner; // -1: not corner, 0-3: Morton-order N-curve corners.
    // Which directional components to coalesce when interface. 
    bool coalesce[8];
  } bc;
private:
  // basic
  void initialize();

  // SRT
  inline double next_fc_srt( double msq, double omega ) const;
  inline double next_fi_srt( std::size_t i, double msq, double omega ) const;
  
  // MRT
  inline double premultiply_M(size_t i) const;
  inline double premultiply_MinvS( 
    size_t i, const double m[9], double omega ) const;
  inline void compute_moment( double m[9] ) const;
  void next_f_mrt( const double m[9], double omega );
  
  // VC
  inline void compute_feq( double feq[9] ) const;
  inline void compute_strain_terms( 
    double& s11, double& s12, double& s22, 
    const double feq[9], double omega ) const;
  inline void fill_strain_terms( double omega );
  inline void compute_strain_differences( 
    double& s11x,double& s12x,double& s12y,double& s22y,double dh_inv ) const;
  inline void compute_vc_body_force( double nuc );
  inline void apply_steady_vc_body_force( 
    double omega, double scale_increase, double scale_decrease, double nuc );
  
  // Dynamic grid
  void link_new_children(size_t child);
  void activate_children( std::vector<Cell>& cg );
  void create_children( std::vector<Cell>& next_level_cells, 
    std::vector<Cell>& grandchild_cells );
};
