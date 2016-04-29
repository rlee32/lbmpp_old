#include "Simulator.h"

using namespace std;

Simulator::Simulator(string filename)
{
  read_settings(filename); 
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
        for ( size_t i = 0; i < 4; ++i )
        {
          if ( not parameter.compare(face_order[i]) ) { iss >> bc[i]; }
        }
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
      relax_model, vc_model );
    if ( u0file != "" and v0file != "" ) read_coarse_solution();
  }
}

void Simulator::iteration()
{
  grid.iteration(0);
}

void Simulator::output_coarse_field(string output_suffix)
{
  ofstream u, v;
  const vector<Cell>& g = grid.get_cells(0);

  u.open("u_"+output_suffix);
  v.open("v_"+output_suffix);
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
  output_coarse_field("test.tsv");
}

// Traverses children to get the data points closest to the top side of 
//  the surface of this cell. 
// x: x-distance from the leftmost side of the root cell.
// y: y-distance from the bottommost side of the root cell.
// u: u-velocity
// v: v-velocity
// vector<> get_topmost_data(double& x, double& y, double& u, double &v,)
// {
//   if ()
//   {

//   }
// }


// // Usable on our tree dynamic grid.
// vector<double> centerline_x()
// {
//   vector<double> values;
//   size_t max_level = 0;
//   for( size_t i = 0; i < cell_count[0]; ++i )
//   {
//     size_t ii = i + 
//     size_t test_level = 
//   }
// }