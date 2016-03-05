// Author: Robert Lee
// Email: rlee32@gatech.edu

#include <cstdlib>
#include <iostream>

#include "viz/CImg.h"
#include "viz/SolutionViewer.h"
#include "sim/Simulator.h"

using namespace std;
using namespace cimg_library;

int main(int argc, char ** argv)
{
  // First read inputs and instantiate control panel.
  Simulator sim("settings");

  // Instantiate the solution viewer.
  int max_pixel_dim = (argc < 2) ? 800 : atoi(argv[1]);
  SolutionViewer sv(sim.grid, max_pixel_dim);

  // Solution loop.
  while ( not sv.window.is_closed() )
  {
    sv.window.wait();
    if (sv.window.button() && sv.window.mouse_y()>=0)
    {
      sv.draw_velocity_magnitude(sim.grid);
      sv.display();
    }
  }
  return EXIT_SUCCESS;
}