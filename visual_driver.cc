// Author: Robert Lee
// Email: rlee32@gatech.edu

#include <cstdlib>
#include <iostream>

#include "viz/CImg.h"
#include "viz/SolutionViewer.h"
#include "sim/Simulator.h"
#include "Timer.h"

using namespace std;
using namespace cimg_library;

int main(int argc, char ** argv)
{
  // First read inputs and instantiate control panel.
  Simulator sim("settings");

  // Instantiate the solution viewer.
  int max_pixel_dim = (argc < 2) ? 800 : atoi(argv[1]);
  SolutionViewer sv(sim, max_pixel_dim);

  // timer
  Timer timer;

  // Solution loop.
  // cout << "Running " << sim.timesteps << " timesteps." << endl;
  // bool done = false;
  int k = 0;
  timer.start();
  while ( not sv.window.is_closed() and k <= sim.get_timesteps() )
  {
    //if (not done)
    {
      // for (int k = 0; k < sim.timesteps; ++k)
      {
        // sv.window.wait();
        // if (sv.window.button() && sv.window.mouse_y()>=0)
        {
          // cout << "Performing iteration " << k << endl;
          sim.iteration();
          if ( k % sim.get_display_interval() == 0 )
          {
            sv.draw_velocity_magnitude( sim.grid );
            sv.draw_status( k, sim, timer.stop() );
            sv.display();
          }
          ++k;
        }
        // else
        //{
        //  --k;
        //}
      }
      //done = true;
    }
    // cin.ignore();
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
    sv.save_image( "mag_"+case_name+".png" );
  }

  cout << "Simulation finished! Press enter to continue." << endl;
  std::cin.ignore();

  return EXIT_SUCCESS;
}
