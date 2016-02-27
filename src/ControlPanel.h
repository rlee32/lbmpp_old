#pragma once

// Contains settings for the simulation.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class ControlPanel
{ 
  public:
    ControlPanel();
    void read_settings(std::string filename);
  private:
    double coarse_cell_dim; // dimension of the square coarsest cell.
    int coarse_cells_x, coarse_cells_y; // coarse cells in the x and y direction.
    double lattice_dt; // discrete, Lattice Boltzmann time step.
    int timesteps;

    double physical_viscosity;
    double physical_velocity;
    double physical_length;

    bool refinement;
    double buffer_viscosity;

    double Re; // Reynolds number.
    double rho;
} ControlPanel;