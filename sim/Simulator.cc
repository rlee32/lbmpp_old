#include "Simulator.h"

using namespace std;

Simulator::Simulator(string filename) :
  refinement(false), timesteps(0), relax_model(1), vc_model(0),
  Re(0),M(0),U(0),nu(0),L(0),tau(0),omega(0),rho0(0),u0(0),v0(0),
  u0file(""),v0file(""),
  M0(0),nuc(0),nucf(0),
  display_interval(1),
  picset(0)
{
  bc[0] = 'w';
  bc[1] = 'w';
  bc[2] = 'w';
  bc[3] = 'w';
  face_order[0] = "bottom";
  face_order[1] = "right";
  face_order[2] = "top";
  face_order[3] = "left";
  face_order_char[0] = 'b';
  face_order_char[1] = 'r';
  face_order_char[2] = 't';
  face_order_char[3] = 'l';
  cell_count[0] = 0;
  cell_count[1] = 0;
  read_settings(filename);
  string grid_string = "G"+to_string( (long long)cell_count[0] );
  if( cell_count[0] != cell_count[1] )
  {
    grid_string += "x"+to_string( (long long)cell_count[1] );
  }
  string mach_string = "M"+to_string( (long long)(M*1000.0) );
  string timesteps_string = "T"+to_string( (long long)(timesteps/1000) );
  string relax_model_string = "RM"+to_string( (long long)(relax_model) );
  string vc_model_string = "VCM"+to_string( (long long)(vc_model) );
  string nucf_string = "VCF"+to_string( (long long)round((nuc/nu)*10.0) );
  string re_string = "Re"+to_string( (long long)round(Re) );
  output_suffix = grid_string+"_"+mach_string+"_"+timesteps_string+"_"
    +relax_model_string+"_"+vc_model_string+"_"+nucf_string+"_"+re_string;
}

void Simulator::read_settings(string filename)
{
  // Let's begin to read settings.
  ifstream settings_file( filename );
  string line;
  if(settings_file.is_open())
  {
    // Reading lines.
    while( getline(settings_file, line) )
    {
      istringstream iss(line);
      string parameter;
      // If content available
      if( (iss >> parameter) )
      {
        if (parameter[0] == '#') continue;
        // Assign the proper parameter entry.
        if ( not parameter.compare("coarse_cells") ) 
        {
          iss >> cell_count[0];
          iss >> cell_count[1];
        }
        if ( not parameter.compare("M") ) iss >> M;
        if ( not parameter.compare("L") ) iss >> L;
        if ( not parameter.compare("Re") ) iss >> Re;
        if ( not parameter.compare("vc_model") ) iss >> vc_model;
        if ( not parameter.compare("nucf") ) iss >> nucf;
        if ( not parameter.compare("relax_model") ) iss >> relax_model;
        if ( not parameter.compare("timesteps") ) iss >> timesteps;
        if ( not parameter.compare("refinement") ) iss >> refinement;
        if ( not parameter.compare("rho0") ) iss >> rho0;
        if ( not parameter.compare("u0") ) iss >> u0;
        if ( not parameter.compare("v0") ) iss >> v0;
        if ( not parameter.compare("u0file") ) iss >> u0file;
        if ( not parameter.compare("v0file") ) iss >> v0file;
        if ( not parameter.compare("M0") ) iss >> M0;
        if ( not parameter.compare("display_interval") ) iss >> display_interval;
        if ( not parameter.compare("picset") ) iss >> picset;
        for ( size_t i = 0; i < 4; ++i )
        {
          if ( not parameter.compare(face_order[i]) ) { iss >> bc[i]; }
        }
        if ( not parameter.compare("experimental") ) iss >> experimental;
      }
    }
    U = M / sqrt(3);
    nu = U * L / Re;
    nuc = ( vc_model != 0 ) ? (nucf*nu) : 0; // counteracting viscosity.
    tau = 3 * (nu + nuc) + 0.5;
    omega = 1.0 / tau;
    grid.initialize( 
      cell_count[0], cell_count[1], 
      rho0, u0, v0, 
      nu, nuc, 
      face_order_char, bc, U,
      relax_model, vc_model, experimental );
    if ( u0file != "" and v0file != "" ) read_coarse_solution();
  }
}

void Simulator::iteration()
{
  grid.iteration(0);
}

void Simulator::output_solution()
{
  output_coarse_field();
  output_centerlines();
}

void Simulator::output_coarse_field()
{
  ofstream u, v;
  const vector<Cell>& g = grid.get_cells(0);

  u.open("u_"+output_suffix+".dat");
  v.open("v_"+output_suffix+".dat");
  u.precision(std::numeric_limits< double >::digits10);
  v.precision(std::numeric_limits< double >::digits10);
  for (size_t j = 0; j < cell_count[1]; ++j)
  {
    for (size_t i = 0; i < cell_count[0]; ++i)
    {
      size_t ii = i + j * cell_count[0];
      u << g[ii].state.u << "\t";
      v << g[ii].state.v << "\t";
    }
    u << endl;
    v << endl;
  }
  u.close();
  v.close();
}

void Simulator::output_picset_field(size_t k)
{
  ofstream u, v;
  const vector<Cell>& g = grid.get_cells(0);
  u.open("picset/u_"+output_suffix+"_"+to_string((long long)k)+".dat");
  v.open("picset/v_"+output_suffix+"_"+to_string((long long)k)+".dat");
  u.precision(std::numeric_limits< double >::digits10);
  v.precision(std::numeric_limits< double >::digits10);
  for (size_t j = 0; j < cell_count[1]; ++j)
  {
    for (size_t i = 0; i < cell_count[0]; ++i)
    {
      size_t ii = i + j * cell_count[0];
      u << g[ii].state.u << "\t";
      v << g[ii].state.v << "\t";
    }
    u << endl;
    v << endl;
  }
  u.close();
  v.close();
}

// Read in the components of a row-major coarse grid 
//  and interpolate values to current grid.
// filename: name of file where each line is a row of cells.
//  Example: 256, 128 is represented by 256 columns, 128 lines. 
// Returns x_cells.
size_t Simulator::read_coarse_field(string filename, vector<double>& phi, 
  double scale)
{
  size_t x_cells = 0;
  ifstream myfile (filename);
  if (myfile.is_open())
  {
    string line;
    bool first_line = true;
    while ( getline (myfile,line) )
    {
      istringstream iss(line);
      double value = 0;
      while( iss >> value )
      {
        phi.push_back( scale*value );
      }
      if ( first_line )
      {
        x_cells = phi.size();
        first_line = false;
      }
    }
    myfile.close();
    if ( phi.size() % x_cells != 0 )
    {
      cout << "Error! Read " << phi.size() << " total cells, but x cells = " 
        << x_cells << endl;
    }
    else
    {
      cout << "Read in " << filename 
        << " (" << x_cells << "x" << phi.size()/x_cells << ")" << endl;
    }
  }
  else
  {
    cout << "Unable to open " << filename << endl; 
  }
  return x_cells;
}

void Simulator::interpolate_field( size_t source_x_cells, vector<double>& source, 
  size_t target_x_cells, size_t target_y_cells, vector<double>& target )
{
  // Assume cell-centered fields.
  // Assume uniformly-sized square cells.
  // Target indices are converted to source indices, then bilinear interpolation
  //  is applied (as on Wikipedia page for bilinear interpolation).

  // Converting target indices to source indices for interpolation.
  double x_convert = ((double)source_x_cells) / ((double)target_x_cells);
  size_t source_y_cells = source.size() / source_x_cells;
  double y_convert = ((double)source_y_cells) / ((double)target_y_cells);
  for( size_t j = 0; j < target_y_cells; ++j )
  {
    double y = ( 0.5 + j ) * y_convert;
    size_t jj = ceil( y - 0.5 );
    double dy = y - floor(y);
    if ( jj >= source_y_cells ) 
    {
      jj = source_y_cells - 1;
      dy = 1.0;
    }
    if ( jj <= 0 ) 
    {
      jj = 1;
      dy = 0;
    }
    for( size_t i = 0; i < target_x_cells; ++i )
    {
      double x = ( 0.5 + i ) * x_convert;
      size_t ii = ceil( x - 0.5 );
      double dx = x - floor(x);
      if ( ii >= source_x_cells ) 
      {
        ii = source_x_cells - 1;
        dx = 1.0;
      }
      if ( ii <= 0 ) 
      {
        ii = 1;
        dx = 0;
      }
      // serial index for upper-right-most point
      size_t si = ii + jj*source_x_cells; 
      double f11 = source[si - 1 - source_x_cells];
      double f12 = source[si - 1];
      double f21 = source[si - source_x_cells];
      double f22 = source[si];
      target.push_back( 
        f11*(1-dx)*(1-dy) 
        + f21*(dx)*(1-dy) 
        + f12*(1-dx)*(dy) 
        + f22*(dx)*(dy) );
    }
  }
}

// Read in the complete macroscopic state rho, u, and v velocity components 
//  of a row-major coarse grid and interpolate values to current grid.
void Simulator::read_coarse_solution()
{
  vector<double> u_source;
  vector<double> v_source;
  vector<size_t> x_cells_candidates;
  vector<size_t> total_cells_candidates;
  x_cells_candidates.push_back( read_coarse_field(u0file, u_source, M / M0 ) );
  x_cells_candidates.push_back( read_coarse_field(v0file, v_source, M / M0 ) );
  size_t x_cells = *( max_element(
    x_cells_candidates.begin(), x_cells_candidates.end()) );
  total_cells_candidates.push_back( u_source.size() );
  total_cells_candidates.push_back( v_source.size() );
  size_t total_cells = *( max_element(
    total_cells_candidates.begin(), total_cells_candidates.end()) );
  size_t y_cells = total_cells / x_cells;
  if ( y_cells*x_cells != total_cells )
  {
    cout << "Error! Cell dimensions not consistent with total cells: "
      << x_cells << "x" << y_cells << " != " << total_cells << endl;
    return;
  }
  // Now kith, I mean interpolate.
  interpolate_field( x_cells, u_source, 
    cell_count[0], cell_count[1], u0field );
  interpolate_field( x_cells, v_source, 
    cell_count[0], cell_count[1], v0field );
  // cout << u_target.size() << ", " << v_target.size() << endl;
  // Finally, add to the grid.
  grid.set_coarse_solution( rho0, u0field, v0field );
  // output_coarse_field("test.tsv");
}


// the public method
void Simulator::output_centerlines()
{
  vector<CellData> top;
  vector<CellData> bottom;
  vector<CellData> right;
  vector<CellData> left;
  centerline_x( left, right );
  centerline_y( top, bottom );
  vector<CellData> centerx;
  vector<CellData> centery;
  produce_centerline_y( top, bottom, centery );
  produce_centerline_x( left, right, centerx );
  print_centerlines( centerx, centery );
}


// Inherently sequential.
// Traverses children to get the data points closest to a particular edge of 
//  this cell. 
// cell: cell under current investigation.
// xstart, ystart: the x,y coordinates of the lower left corner of this cell.  
// dim: dimension of current cell.
void Simulator::get_data( Cell& cell, 
  vector<CellData>& data, char side, 
  double xstart, double ystart, double dim )
{
  if ( cell.active() )
  {
    CellData cd;
    cd.u = cell.u();
    cd.v = cell.v();
    cd.x = xstart + 0.5*dim;
    cd.y = ystart + 0.5*dim;
    data.push_back( cd );
  }
  else
  {
    dim *= 0.5;
    double xnew = xstart;
    double ynew = ystart;
    switch ( side )
    {
      case 't': // top
        ynew += dim;
        get_data( cell(1), data, side, xnew, ynew, dim );
        xnew += dim;
        get_data( cell(3), data, side, xnew, ynew, dim );
        break;
      case 'b': // bottom
        get_data( cell(0), data, side, xnew, ynew, dim );
        xnew += dim;
        get_data( cell(2), data, side, xnew, ynew, dim );
        break;
      case 'l': // left
        get_data( cell(0), data, side, xnew, ynew, dim );
        ynew += dim;
        get_data( cell(1), data, side, xnew, ynew, dim );
        break;
      case 'r': // right
        xnew += dim;
        get_data( cell(2), data, side, xnew, ynew, dim );
        ynew += dim;
        get_data( cell(3), data, side, xnew, ynew, dim );
        break;
      case 'm': // middle should be taken care of outside
        break;
      default: break;
    }

  }
}

// Usable on our tree dynamic grid.
// Assumes lid-driven cavity, with top surface as moving lid.
// Assumes square domain with square cells.
void Simulator::centerline_y(
  vector<CellData>& top_values, vector<CellData>& bottom_values)
{
  CellData left_wall;
  left_wall.u = 0;
  left_wall.v = 0;
  left_wall.x = 0;
  left_wall.y = 0.5;
  top_values.push_back(left_wall);
  bottom_values.push_back(left_wall);
  double dim = 1.0 / cell_count[0];
  // odd
  if ( cell_count[1] & 1 )
  {
    size_t start_j = cell_count[1] / 2;
    for( size_t i = 0; i < cell_count[0]; ++i )
    {
      size_t ii = start_j*cell_count[0] + i;
      Cell& cell = grid[0][ii];
      if ( cell.active() )
      {
        // right on centerline
        CellData cd;
        cd.u = cell.u();
        cd.v = cell.v();
        cd.x = (i+0.5)*dim;
        cd.y = 0.5;
        top_values.push_back(cd);
        bottom_values.push_back(cd);
      }
      else
      {
        // below centerline
        get_data( cell(0), bottom_values, 't', i*dim, 0.5-0.5*dim, 0.5*dim );
        get_data( cell(2), bottom_values, 't', (i+0.5)*dim, 0.5-0.5*dim, 0.5*dim );
        // above centerline
        get_data( cell(1), top_values, 'b', i*dim, 0.5, 0.5*dim );
        get_data( cell(3), top_values, 'b', (i+0.5)*dim, 0.5, 0.5*dim );
      }
    }
  }
  // even
  else
  {
    size_t top_j = cell_count[1] >> 1;
    size_t bottom_j = top_j - 1;
    // below centerline
    for( size_t i = 0; i < cell_count[0]; ++i )
    {
      size_t ii = bottom_j * cell_count[0] + i;
      double xstart = i*dim;
      double ystart = bottom_j*dim;
      get_data( grid[0][ii], bottom_values, 't', xstart, ystart, dim );
    }
    // above centerline
    for( size_t i = 0; i < cell_count[0]; ++i )
    {
      size_t ii = top_j * cell_count[0] + i;
      double xstart = i*dim;
      double ystart = top_j*dim;
      get_data( grid[0][ii], top_values, 'b', xstart, ystart, dim );
    }
  }
  CellData right_wall;
  right_wall.u = 0;
  right_wall.v = 0;
  right_wall.x = 1.0;
  right_wall.y = 0.5;
  top_values.push_back(right_wall);
  bottom_values.push_back(right_wall);
}
// Usable on our tree dynamic grid.
// Assumes lid-driven cavity, with top surface as moving lid.
// Assumes square domain with square cells.
void Simulator::centerline_x(
  vector<CellData>& left_values, vector<CellData>& right_values)
{
  CellData bottom_wall;
  bottom_wall.u = 0;
  bottom_wall.v = 0;
  bottom_wall.x = 0.5;
  bottom_wall.y = 0;
  left_values.push_back(bottom_wall);
  right_values.push_back(bottom_wall);
  double dim = 1.0 / cell_count[1];
  // odd
  if ( cell_count[0] & 1 )
  {
    size_t start_i = cell_count[0] >> 1;
    for( size_t j = 0; j < cell_count[1]; ++j )
    {
      size_t ii = j*cell_count[0] + start_i;
      Cell& cell = grid[0][ii];
      if ( cell.active() )
      {
        // right on centerline
        CellData cd;
        cd.u = cell.u();
        cd.v = cell.v();
        cd.x = 0.5;
        cd.y = (j+0.5)*dim;
        left_values.push_back(cd);
        right_values.push_back(cd);
      }
      else
      {
        // left of centerline
        get_data( cell(0), left_values, 'r', 0.5-0.5*dim, j*dim, 0.5*dim );
        get_data( cell(1), left_values, 'r', 0.5-0.5*dim, (j+0.5)*dim, 0.5*dim );
        // right of centerline
        get_data( cell(2), right_values, 'l', 0.5, j*dim, 0.5*dim );
        get_data( cell(3), right_values, 'l', 0.5, (j+0.5)*dim, 0.5*dim );
      }
    }
  }
  // even
  else
  {
    size_t right_i = cell_count[0] >> 1;
    size_t left_i = right_i - 1;
    // left of centerline
    for( size_t j = 0; j < cell_count[1]; ++j )
    {
      size_t ii = j * cell_count[0] + left_i;
      double xstart = left_i*dim;
      double ystart = j*dim;
      get_data( grid[0][ii], left_values, 'r', xstart, ystart, dim );
    }
    // left of centerline
    for( size_t j = 0; j < cell_count[1]; ++j )
    {
      size_t ii = j * cell_count[0] + right_i;
      double xstart = right_i*dim;
      double ystart = j*dim;
      get_data( grid[0][ii], right_values, 'l', xstart, ystart, dim );
    }
  }
  CellData top_wall;
  top_wall.u = U;
  top_wall.v = 0;
  top_wall.x = 0.5;
  top_wall.y = 1.0;
  left_values.push_back(top_wall);
  right_values.push_back(top_wall);
}

// From side1 and side2 vector data, produces centerline data
// we assume 0.5 is the centerline for x or y
void Simulator::produce_centerline_y(vector<CellData>& side1, 
  vector<CellData>& side2, vector<CellData>& center)
{
  CellData cd;
  cd.u = 0;
  cd.v = 0;
  cd.x = 0;
  cd.y = 0.5;
  while( not side1.empty() and not side2.empty() )
  {
    double x1 = side1.back().x;
    double x2 = side2.back().x;
    double y1 = side1.back().y;
    double y2 = side2.back().y;
    if( x1 == x2 )
    {
      cd.x = x1;
      if(y1 == y2)
      {
        cd.u = side1.back().u;
        cd.v = side1.back().v;
      }
      else
      {
        cd.u = ( side1.back().u + side2.back().u ) * 0.5;
        cd.v = ( side1.back().v + side2.back().v ) * 0.5;
      }
      center.push_back(cd);
      side1.pop_back();
      side2.pop_back();
    }
    else
    {
      double dx = x2 - x1;
      double dy = y2 - y1;
      double mag = sqrt( dx*dx + dy*dy );
      double dy1 = 0.5 - y1;
      double dx1 = dy1 / dy * dx;
      double mag1 = sqrt( dx1*dx1 + dy1*dy1 );
      double r2 = mag1 / mag;

      cd.x = ( 1.0 - r2 ) * x1 + r2 * x2;
      cd.u = ( 1.0 - r2 ) * side1.back().u + r2 * side2.back().u;
      cd.v = ( 1.0 - r2 ) * side1.back().v + r2 * side2.back().v;

      center.push_back(cd);
      if( x1 < x2 ) side1.pop_back();
      if( x2 < x1 ) side2.pop_back();
    }
  }
}
// From side1 and side2 vector data, produces centerline data
// we assume 0.5 is the centerline for x or y
void Simulator::produce_centerline_x(vector<CellData>& side1, 
  vector<CellData>& side2, vector<CellData>& center)
{
  CellData cd;
  cd.u = 0;
  cd.v = 0;
  cd.x = 0.5;
  cd.y = 0;
  while( not side1.empty() and not side2.empty() )
  {
    double x1 = side1.back().x;
    double x2 = side2.back().x;
    double y1 = side1.back().y;
    double y2 = side2.back().y;
    if( y1 == y2 )
    {
      cd.y = y1;
      if( x1 == x2 )
      {
        cd.u = side1.back().u;
        cd.v = side1.back().v;
      }
      else
      {
        cd.u = ( side1.back().u + side2.back().u ) * 0.5;
        cd.v = ( side1.back().v + side2.back().v ) * 0.5;
      }
      center.push_back(cd);
      side1.pop_back();
      side2.pop_back();
    }
    else
    {
      double dx = x2 - x1;
      double dy = y2 - y1;
      double mag = sqrt( dx*dx + dy*dy );
      double dx1 = 0.5 - x1;
      double dy1 = dx1 / dx * dy;
      double mag1 = sqrt( dx1*dx1 + dy1*dy1 );
      double r2 = mag1 / mag;

      cd.y = ( 1.0 - r2 ) * y1 + r2 * y2;
      cd.u = ( 1.0 - r2 ) * side1.back().u + r2 * side2.back().u;
      cd.v = ( 1.0 - r2 ) * side1.back().v + r2 * side2.back().v;

      center.push_back(cd);
      if( y1 > y2 ) side1.pop_back();
      if( y2 > y1 ) side2.pop_back();
    }
  }
}

// outputs a file with centerline info. these are the line descriptions:
// x positions of y-centerline
// v values of y-centerline
// y positions of x-centerline
// u values of x-centerline
void Simulator::print_centerlines( 
  vector<CellData>& centerx, vector<CellData>& centery )
{
  ofstream out;
  out.open("centerlines_"+output_suffix+".tsv");
  out.precision(std::numeric_limits< double >::digits10);
  vector<CellData>::iterator it = centery.begin();
  for ( ; it != centery.end(); ++it ) out << it->x << "\t";
  out << endl;
  it = centery.begin();
  for ( ; it != centery.end(); ++it ) out << it->v << "\t";
  out << endl;
  it = centerx.begin();
  for ( ; it != centerx.end(); ++it ) out << it->y << "\t";
  out << endl;
  it = centerx.begin();
  for ( ; it != centerx.end(); ++it ) out << it->u << "\t";
  out << endl;
  out.close();
}


