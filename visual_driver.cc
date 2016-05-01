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
            if ( sim.do_picset() )
            {
              sim.output_picset_field( k/sim.get_display_interval() );
            }
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
  
  if ( k >= sim.get_timesteps() )
  {
    sim.output_solution();
    sv.save_image( "mag_"+sim.get_output_suffix()+".png" );
  }

  cout << "Simulation finished! Press enter to continue." << endl;
  std::cin.ignore();

  return EXIT_SUCCESS;
}
