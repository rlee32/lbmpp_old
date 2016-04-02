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
  SolutionViewer sv(sim.grid, max_pixel_dim);

  // timer
  Timer timer;

  // Solution loop.
  // cout << "Running " << sim.timesteps << " timesteps." << endl;
  // bool done = false;
  int k = 0;
  int display_interval = 1000;
  timer.start();
  while ( not sv.window.is_closed() )
  {
    //if (not done)
    {
      // for (int k = 0; k < sim.timesteps; ++k)
      {
        // sv.window.wait();
        // if (sv.window.button() && sv.window.mouse_y()>=0)
        {
          // cout << "Performing iteration " << k << endl;
          sim.iterate();
          if ( k % display_interval == 0 )
          {
            sv.draw_velocity_magnitude(sim.grid);
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
  }
  sim.output_coarse_field("velocity_field.dat");

  return EXIT_SUCCESS;
}
