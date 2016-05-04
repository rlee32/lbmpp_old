#include "BoundaryConditions.h"

using namespace std;

BoundaryConditions::BoundaryConditions() 
: U(0), next_level_bcs(nullptr), g(nullptr)
{
  
}

// We assume that the first and last cells in Face::cells are the corners!!!
void BoundaryConditions::apply_bc()
{
  // for (size_t i = 0; i < 4; ++i)
  // {
  //   Face& f = faces[i];
  //   if (f.cells.size() > 0)
  //   {
  //     // Corner cells.
  //     // Applying moving wall to corners can dramatically reduce maximum Re!
  //     // So we apply the neighbouring BCs here.
  //     if (f.type == 'm')
  //     {
  //       cell_bc( f.cells[0], get_prev_type(f.side), f.side );
  //       cell_bc( f.cells[f.cells.size()-1], get_next_type(f.side), f.side );
  //     }
  //     else
  //     {
  //       cell_bc( f.cells[0], f.type, f.side );
  //       cell_bc( f.cells[f.cells.size()-1], f.type, f.side );
  //     }
  //     // cout << "Non-corner " << f.cells.size() << endl;
  //     // Non-corner cells.
  //     for( size_t i = 1; i < f.cells.size()-1; ++i )
  //     {
  //       // cout << i << ", " << f.type << ", " << f.side << endl;
  //       cell_bc( f.cells[i], f.type, f.side );
  //     }
  //   }
  // }
  for (size_t i = 0; i < 4; ++i)
  {
    Face& f = faces[i];
    if (f.cells.size() > 0)
    {
      // Corner cells.
      // Applying moving wall to corners can dramatically reduce maximum Re!
      // So we apply the neighbouring BCs here.
      if (f.type == 'm')
      {
        for( size_t i = 0; i < f.cells.size(); ++i )
        {
          char t = ( (*g)[f.cells[i]].bc.corner > -1 ) ? 'w' : f.type;
          cell_bc(f.cells[i],t,f.side);
        }
      }
      else
      {
        for( size_t i = 0; i < f.cells.size(); ++i )
        {
          // cout << i << ", " << f.type << ", " << f.side << endl;
          cell_bc( f.cells[i], f.type, f.side );
        }
      }
      // cout << "Non-corner " << f.cells.size() << endl;
      // Non-corner cells.
    }
  }
}

// After refinement, we add the refined cells to the next level BCs.
void BoundaryConditions::refined_cell_bc()
{
  for (size_t i = 0; i < 4; ++i)
  {
    Face& f = faces[i];
    // cout << "Adding " << f.cells.size() << "to " << f.side << endl;
    for( size_t i = 0; i < f.cells.size(); ++i )
    {
      // cout << "split side: " << f.side << endl;
      facilitate_split( f.cells[i], f.side );
    }
  }
}

// After a cell has been split, add the appropriate children to the 
//  boundary condition level.
void BoundaryConditions::facilitate_split( int ci, char side )
{
  Cell& c = (*g)[ ci ];
  if ( c.action.refine or c.action.link_children )
  {
    if ( c.bc.corner > -1 ) c(c.bc.corner).bc.corner = c.bc.corner;
    switch( side )
    {
      case 'l':
        next_level_bcs->add_cell( c.local.children[0], side );
        next_level_bcs->add_cell( c.local.children[1], side );
        break;
      case 'r':
        next_level_bcs->add_cell( c.local.children[2], side );
        next_level_bcs->add_cell( c.local.children[3], side );
        break;
      case 'b':
        next_level_bcs->add_cell( c.local.children[0], side );
        next_level_bcs->add_cell( c.local.children[2], side );
        break;
      case 't':
        next_level_bcs->add_cell( c.local.children[1], side );
        next_level_bcs->add_cell( c.local.children[3], side );
        break;
      default:
        cout << "Error! Unknown side added." << endl;
        break;
    }
  }
}

// 
void BoundaryConditions::cell_bc( int ci, char type, char side ) const
{
  Cell& c = (*g)[ ci ];
  // if ( c.state.active or c.state.interface )
  if ( c.state.active )
  {
    switch( type )
    {
      case 'm':
        c.moving_wall(side, U);
        break;
      case 'w':
        c.bounce_back(side);
        break;
      default:
        break;
    }
  }
}

void BoundaryConditions::initialize( char sides[4], char bc[4], double U_, 
  BoundaryConditions* next_level_bcs_, vector<Cell>* g_ )
{
  U = U_;
  for( size_t i = 0; i < 4; ++i )
  {
    Face face( sides[i], bc[i] );
    faces.push_back( face );
  }
  sort( faces.begin(), faces.end(), Face::compare );
  next_level_bcs = next_level_bcs_;
  g = g_;
}

// Slightly inefficient search and add procedure, but this should only be used
//  when initializing the coarse cells.
void BoundaryConditions::add_cell( int c, char side)
{
  for( size_t i = 0; i < 4; ++i )
  {
    Face& f = faces[i];
    if ( f.side == side )
    {
      f.cells.push_back(c);
      break;
    }
  }
}

// We assume that the first and last cells in Face::cells are the corners.
// We assume 'next' corresponds to the side that touches the last i or j index.
// 'prev' corresponds to the side that touches the first i or j index.
char BoundaryConditions::next_side(char side) const
{
  switch( side )
  {
    case 'r':
      return 't';
    case 'l':
      return 't';
    case 't':
      return 'r';
    case 'b':
      return 'r';
    default:
      return 'x';
  }
}
char BoundaryConditions::prev_side(char side) const
{
  switch( side )
  {
    case 'r':
      return 'b';
    case 'l':
      return 'b';
    case 't':
      return 'l';
    case 'b':
      return 'l';
    default:
      return 'x';
  }
}
char BoundaryConditions::get_next_type(char side) const
{
  char next = next_side(side);
  for(size_t i = 0; i < 4; ++i)
  {
    if(faces[i].side == next)
    {
      return faces[i].type; 
    }
  }
  return 'x';
}
char BoundaryConditions::get_prev_type(char side) const
{
  char prev = prev_side(side);
  for(size_t i = 0; i < 4; ++i)
  {
    if(faces[i].side == prev)
    {
      return faces[i].type; 
    }
  }
  return 'x';
}

// Determines which boundary type should dominate when two exist (at corner).
// Higher means it should be applied last (so its effects dominate).
int BoundaryConditions::Face::type_rank( char type ) const
{
  switch( type )
  {
    case 'm':
      return 0;
    case 'o':
      return 2;
    case 'i':
      return 1;
    case 'w':
      return 3;
    default:
      return -1;
  }
}


