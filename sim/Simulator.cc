#include "Simulator.h"

using namespace std;

Simulator::Simulator(string filename)
{
  read_settings(filename); 
  process_settings();
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
        if ( not parameter.compare("coarse_cell_size") ) iss >> coarse_cell_size;
        if ( not parameter.compare("coarse_cells_x") ) iss >> cell_count[0];
        if ( not parameter.compare("coarse_cells_y") ) iss >> cell_count[1];
        if ( not parameter.compare("viscosity_physical") ) iss >> viscosity_physical;
        if ( not parameter.compare("velocity_physical") ) iss >> velocity_physical;
        if ( not parameter.compare("length_physical") ) iss >> length_physical;
        if ( not parameter.compare("dt_lattice") ) iss >> dt_lattice;
        if ( not parameter.compare("buffer_viscosity_factor") ) iss >> buffer_viscosity_factor;
        if ( not parameter.compare("timesteps") ) iss >> timesteps;
        if ( not parameter.compare("refinement") ) iss >> refinement;
        if ( not parameter.compare("rho0") ) iss >> rho0;
        if ( not parameter.compare("u0") ) iss >> u0;
        if ( not parameter.compare("v0") ) iss >> v0;
        if ( not parameter.compare("bottom") ) 
        {
          iss >> bc[0];
          iss >> bcv[0];
        }
        if ( not parameter.compare("right") )
        {
          iss >> bc[1];
          iss >> bcv[1];
        }
        if ( not parameter.compare("top") )
        {
          iss >> bc[2];
          iss >> bcv[2];
        }
        if ( not parameter.compare("left") )
        {
          iss >> bc[3];
          iss >> bcv[3];
        }
        // cout << parameter << endl;
      }
    }
    // cout << coarse_cell_size << endl;
    // cout << cell_count[0] << endl;
    // cout << cell_count[1] << endl;
    // cout << viscosity_physical << endl;
    // cout << velocity_physical << endl;
    // cout << length_physical << endl;
    // cout << dt_lattice << endl;
    // cout << buffer_viscosity_factor << endl;
    // cout << timesteps << endl;
    // cout << refinement << endl;
  }
}

void Simulator::process_settings()
{
  Re = length_physical * velocity_physical / viscosity_physical;
  cout << "Reynolds number: " << Re << endl;
  dt_physical = length_physical / velocity_physical;
  viscosity_lattice = dt_lattice / coarse_cell_size / coarse_cell_size / Re;
  cout << "Lattice viscosity: " << viscosity_lattice << endl;
  double nuc = buffer_viscosity_factor*viscosity_lattice; // counteracting viscosity.
  tau = 3.0 * (viscosity_lattice + nuc) + 0.5;
  velocity_lattice = dt_lattice / coarse_cell_size;
  omega = 1.0 / tau;
  cout << "Let there be grid! (creating coarse grid)" << endl;
  grid.initialize(cell_count[0], cell_count[1], rho0, u0, v0, 
  tau, omega, viscosity_lattice, nuc, bc, bcv, velocity_lattice );
}

void Simulator::iterate()
{
  grid.iterate(0);
}

// void Simulator::run()
// {
// }
