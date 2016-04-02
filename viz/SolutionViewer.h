#pragma once

#include <string>
#include <iostream>

#include "CImg.h"
#include "../sim/Simulator.h"
#include "../grid/Grid.h"

using namespace cimg_library;

class SolutionViewer
{
public:
  SolutionViewer(Grid& grid, int max_pixel_dim);
  void test_draw();
  void display();
  void draw_velocity_magnitude(Grid& grid);
  void draw_status( int iteration, Simulator& sim, double elapsed_time );
  CImgDisplay window;
private:
  static const uint TextDisplayDim = 70;
  static const uint TextHeight = 10;
  static const uint TextPadding = 3;
  int pixels[2]; // The pixel dimension of the windows.
  int pixels_per_cell; // The number of pixels per coarsest cell.
  CImg<unsigned char> image;
  void scalar2rgb(double min, double max, double value, float rgb[3]);
  void draw_text_line( std::string text, int line );
};