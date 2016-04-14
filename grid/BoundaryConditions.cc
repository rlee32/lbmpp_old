#include "BoundaryConditions.h"

using namespace std;

void BoundaryConditions::initialize( char sides[4], char bc[4], double U_ )
{
  U = U_;
  for( size_t i = 0; i < 4; ++i )
  {
    Face face( sides[i], bc[i] );
    faces.push_back( face );
  }
  sort( faces.begin(), faces.end(), Face::compare );
}

// Slightly inefficient search and add procedure, but this should only be used
//  when initializing the coarse cells.
void BoundaryConditions::add_cell(Cell* c, char side)
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

// 
void BoundaryConditions::apply_bc()
{
  for (size_t i = 0; i < 4; ++i)
  {
    Face& f = faces[i];
    // cout << "Applying bc " << f.type << ", on " << f.side << endl;
    // if (f.type=='m') cout << U << endl;
    for( size_t i = 0; i < f.cells.size(); ++i )
    {
      cell_bc( f.cells[i], f.type, f.side );
    }
  }
}

// 
void BoundaryConditions::cell_bc(Cell* c, char type, char side)
{
  switch( type )
  {
    case 'm':
      c->moving_wall(side, U);
      break;
    case 'w':
      c->bounce_back(side);
      break;
    default:
      break;
  }
}

// Determines which boundary type should dominate when two exist (at corner).
// Higher means it should be applied last (so its effects dominate).
int BoundaryConditions::Face::type_rank( char type )
{
  switch( type )
  {
    case 'm':
      return 3;
    case 'o':
      return 2;
    case 'i':
      return 1;
    case 'w':
      return 0;
    default:
      return -1;
  }
}


