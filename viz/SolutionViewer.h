#pragma once

#include <iostream>

#include "CImg.h"
#include "../grid/Grid.h"

using namespace cimg_library;

class SolutionViewer
{
public:
  SolutionViewer(Grid& grid, int max_pixel_dim);
  void test_draw();
  void display();
  void draw_velocity_magnitude(Grid& grid);
  CImgDisplay window;

private:
  int pixels[2]; // The pixel dimension of the windows.
  int pixels_per_cell; // The number of pixels per coarsest cell.
  CImg<unsigned char> image;
  void scalar2rgb(double min, double max, double value, float rgb[3]);
};