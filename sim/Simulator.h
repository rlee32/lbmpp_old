#pragma once

// Contains settings for the simulation.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>

#include "../grid/Grid.h"
#include "../grid/GridLevel.h"

// For post-processing purposes.
typedef struct CellData
{
  double u,v;
  double x,y;
  CellData():u(0),v(0),x(0),y(0){}
} CellData;

class Simulator
{ 
public:
  Simulator(std::string filename);
  void iteration(); // performs one time step advancement.
  Grid grid;
  std::size_t get_cell_count_0() const { return cell_count[0]; }
  std::size_t get_cell_count_1() const { return cell_count[1]; }
  double get_Re() const { return Re; }
  double get_nu() const { return nu; }
  double get_tau() const { return tau; }
  double get_M() const { return M; }
  double get_U() const { return U; }
  double get_timesteps() const { return timesteps; }
  double get_relax_model() const { return relax_model; }
  double get_vc_model() const { return vc_model; }
  double get_nucf() const { return nuc / nu; }
  std::size_t get_display_interval() const { return display_interval; }
  void output_solution();
  bool do_picset() { return picset; }
  void output_picset_field(std::size_t k);
  std::string get_output_suffix() { return output_suffix; }
private:
  // Run time control.
  bool refinement; // If true, enables solution-adaptive refinement.
  std::size_t timesteps;
  // Solver parameters.
  std::size_t relax_model;
  std::size_t vc_model;
  // Physical parameters.
  double Re; // Reynolds number.
  double M; // Mach number for characteristic velocity.
  double U; // characteristic velocity, used for velocity BC values.
  double nu; // lattice viscosity of the coarsest cells.
  double L;
  double tau; // relaxation time of the coarsest cells.
  double omega; // relaxation frequency of the coarsest cells.
  // Initial values.
  double rho0;
  double u0;
  double v0;
  std::string u0file;
  std::string v0file;
  std::vector<double> u0field;
  std::vector<double> v0field;
  double M0;
  // For the viscosity-counteracting approach. The buffer viscosity.
  double nuc;
  double nucf;
  // Boundary conditions
  char bc[4]; // bottom, right, top, left
  std::string face_order[4];
  char face_order_char[4];
  // coarse grid dimension
  std::size_t cell_count[2];
  // output control
  std::size_t display_interval;
  bool picset;
  std::string output_suffix;
  static const size_t output_precision = 16;

  bool experimental;

  // Preprocessing 
  void read_settings(std::string filename);
  void process_settings();
  // coarse input solution.
  std::size_t read_coarse_field(std::string filename, std::vector<double>& phi, 
    double scale);
  void interpolate_field( std::size_t source_x_cells, std::vector<double>& source, 
    std::size_t target_x_cells, std::size_t target_y_cells, 
    std::vector<double>& target );
  void read_coarse_solution();

  // Postprocessing
  void output_coarse_field();
  void get_data( Cell& cell, 
    std::vector<CellData>& data, char side, 
    double xstart, double ystart, double dim );
  void produce_centerline_y( std::vector<CellData>& side1, 
    std::vector<CellData>& side2, std::vector<CellData>& center);
  void produce_centerline_x( std::vector<CellData>& side1, 
    std::vector<CellData>& side2, std::vector<CellData>& center);
  void centerline2file( std::vector<CellData>& center );
  void print_centerlines( 
    std::vector<CellData>& centerx, std::vector<CellData>& centery );
  void centerline_x( 
    std::vector<CellData>& left_values, std::vector<CellData>& right_values );
  void centerline_y(
    std::vector<CellData>& top_values, std::vector<CellData>& bottom_values );
  void output_centerlines();
  // Initialization

};
 
