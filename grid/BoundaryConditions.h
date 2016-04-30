#pragma once

#include <vector>
#include <algorithm>

#include "../cell/Cell.h"

class BoundaryConditions
{
public:
  BoundaryConditions();
  void initialize( char sides[4], char bc[4], double U_, 
    BoundaryConditions* next_level_bcs_, std::vector<Cell>* g_ );
  void add_cell( int c, char side );
  void apply_bc();
  void refined_cell_bc();
private:
  double U; // value for velocity boundary conditions.
  // Auxiliary data structure.
  typedef struct Face
  {
    std::vector<int> cells; // indices of cells.
    // 't': top
    // 'b': bottom
    // 'r': right
    // 'l': left
    char side;
    // 'w': wall
    // 'm': tangentially moving wall
    // 'i': velocity inlet
    // 'o': zero-gradient outlet
    char type;
    // This defines the order of applying boundary conditions.
    // For example, a corner with one side wall, and one side moving wall, 
    //  the moving wall should be applied last and be dominant.
    int rank;
    int type_rank( char type ) const;
    Face(char s, char t) : side(s), type(t) { rank = type_rank(type); }
    bool static compare(const Face& a, const Face& b) 
      { return a.rank < b.rank; }
  } Face;
  std::vector<Face> faces;
  BoundaryConditions* next_level_bcs;
  std::vector<Cell>* g;

  void cell_bc( int c, char type, char side ) const;
  char next_side(char side) const;
  char prev_side(char side) const;
  char get_next_type(char side) const;
  char get_prev_type(char side) const;
  void facilitate_split( int c, char side );
};