// Author: Robert Lee
// Email: rlee32@gatech.edu

// This version simply runs to the specified timesteps and writes out results.
// No visualization.

#include <cstdlib>
#include <iostream>

#include "sim/Simulator.h"
#include "Timer.h"

using namespace std;

int main(int argc, char ** argv)
{
  // First read inputs and instantiate control panel.
  Simulator sim("settings");

  // timer
  Timer timer;

  // Solution loop.
  cout << "Running " << sim.get_timesteps() << " timesteps." << endl;
  // bool done = false;
  int k = 0;
  timer.start();
  while ( k <= sim.get_timesteps() )
  {
    sim.iteration();
    if ( k % sim.get_display_interval() == 0 )
    {
      cout << "Average speed: " << k / timer.stop()
        << " iterations / second" << endl;
    }
    ++k;
  }
  
  std::size_t cx = sim.get_cell_count_0();
  std::size_t cy = sim.get_cell_count_1();
  string grid_string = "G"+to_string( cx );
  if( cx != cy )
  {
    grid_string += "x"+to_string( cy );
  }
  string mach_string = "M"+to_string( (size_t)(sim.get_M()*1000.0) );
  string timesteps_string = "T"+to_string( (size_t)(sim.get_timesteps()/1000) );
  string relax_model_string = "RM"+to_string( (size_t)(sim.get_relax_model()) );
  string vc_model_string = "VCM"+to_string( (size_t)(sim.get_vc_model()) );
  string nucf_string = "VCF"+to_string( (size_t)round(sim.get_nucf()*10.0) );
  string re_string = "Re"+to_string( (size_t)round(sim.get_Re()) );
  string case_name = grid_string+"_"+mach_string+"_"+timesteps_string+"_"
    +relax_model_string+"_"+vc_model_string+"_"+nucf_string+"_"+re_string;
  
  if ( k >= sim.get_timesteps() )
  {
    sim.output_solution( case_name );
  }

  cout << "Simulation finished!" << endl;

  return EXIT_SUCCESS;
}
