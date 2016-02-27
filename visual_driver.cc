// Author: Robert Lee
// Email: rlee32@gatech.edu

#include <cstdlib>
#include <iostream>

#include "src/CImg.h"

using namespace std;
using namespace cimg_library;

#define MAXRESDIM 500

int main(int argc, char ** argv)
{
  CImg<unsigned char> solution_image(MAXRESDIM,MAXRESDIM,
    1, // 1 z-layer (2D)
    3, // 3 channels (RGB)
    255); // default background color.
  const float black[] = {0,0,0};
  solution_image.draw_line(1,1,50,50,black);

  CImgDisplay solution_viewer(solution_image, "Flow solution");
  while ( !solution_viewer.is_closed() )
  {
    solution_viewer.wait();
    if (solution_viewer.button() && solution_viewer.mouse_y()>=0)
    {
      solution_image.draw_line(100,100,50,50,black);
      solution_image.draw_fill(100,100,black);
      solution_image.draw_rectangle(100,100,150,150,black);
      solution_image.draw_grid(10,10,0,0,false,false,black);
      solution_image.display(solution_viewer);
    }
  }
  return EXIT_SUCCESS;
}