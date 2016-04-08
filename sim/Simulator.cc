#include "Simulator.h"

using namespace std;

Simulator::Simulator(string filename)
{
  read_settings(filename); 
}

void Simulator::read_settings(string filename)
{
  double nucf = 0;
  double M = 0;
  double L = 0;

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
        if ( not parameter.compare("nucf") ) iss >> nucf;
        if ( not parameter.compare("timesteps") ) iss >> timesteps;
        if ( not parameter.compare("refinement") ) iss >> refinement;
        if ( not parameter.compare("rho0") ) iss >> rho0;
        if ( not parameter.compare("u0") ) iss >> u0;
        if ( not parameter.compare("v0") ) iss >> v0;
        if ( not parameter.compare("bottom") ) 
        {
          iss >> bc[0];
          // iss >> bcv[0];
        }
        if ( not parameter.compare("right") )
        {
          iss >> bc[1];
          // iss >> bcv[1];
        }
        if ( not parameter.compare("top") )
        {
          iss >> bc[2];
          // iss >> bcv[2];
        }
        if ( not parameter.compare("left") )
        {
          iss >> bc[3];
          // iss >> bcv[3];
        }
      }
    }
    U = M / sqrt(3);
    nu = U * L / Re;
    nuc = nucf*nu; // counteracting viscosity.
    // cs2 = uc*uc / 3; // lattice speed squared.
    // tau = (nu + nuc) / cs2 + dt / 2.0;
    tau = 3 * (nu + nuc) + 0.5;
    omega = 1.0 / tau;
    // grid.initialize(cell_count[0], cell_count[1], rho0, u0, v0, 
    // tau, omega, nu, nuc, bc, bcv, uc );
    grid.initialize( cell_count[0], cell_count[1], rho0, u0, v0, 
    tau, omega, nu, nuc, bc, U );
  }
}

void Simulator::iterate()
{
  grid.iterate(0);
}

void Simulator::output_coarse_field(string output_suffix)
{
  ofstream u, v;
  vector<Cell>& g = grid.grid_levels[0];

  u.open("u_"+output_suffix);
  v.open("v_"+output_suffix);
  for (size_t j = 0; j < grid.cell_count[1]; ++j)
  {
    for (size_t i = 0; i < grid.cell_count[0]; ++i)
    {
      u << g[i+j*grid.cell_count[0]].state.u << "\t";
      v << g[i+j*grid.cell_count[0]].state.v << "\t";
    }
    u << endl;
    v << endl;
  }
  u.close();
  v.close();
}
