// Author: Robert Lee
// Email: rlee32@gatech.edu

#include <cstdlib>
#include <iostream>

#include "viz/CImg.h"
#include "viz/SolutionViewer.h"
#include "sim/Simulator.h"

using namespace std;
using namespace cimg_library;

#define MAXRESDIM 500

int main(int argc, char ** argv)
{

  // First read inputs and instantiate control panel.
  Simulator sim("settings");

  // Instantiate the solution viewer.
  SolutionViewer sv(sim.grid, MAXRESDIM);

  while ( not sv.window.is_closed() )
  {
    sv.window.wait();
    if (sv.window.button() && sv.window.mouse_y()>=0)
    {
      sv.test_draw();
      sv.display();
    }
  }
  return EXIT_SUCCESS;
}