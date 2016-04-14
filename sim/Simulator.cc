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
        for ( size_t i = 0; i < 4; ++i )
        {
          if ( not parameter.compare(face_order[i]) ) { iss >> bc[i]; }
        }
      }
    }
    U = M / sqrt(3);
    nu = U * L / Re;
    nuc = nucf*nu; // counteracting viscosity.
    tau = 3 * (nu + nuc) + 0.5;
    omega = 1.0 / tau;
    grid.initialize( 
      cell_count[0], cell_count[1], 
      rho0, u0, v0, 
      nu, nuc, 
      face_order_char, bc, U,
      relax_model, vc_model );
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

