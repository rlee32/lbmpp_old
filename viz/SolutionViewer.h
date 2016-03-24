#pragma once

#include <string>
#include <iostream>

#include "CImg.h"
#include "../grid/Grid.h"

using namespace cimg_library;

// Text display settings.
#define TEXT_DISPLAY_DIM 50
#define TEXT_HEIGHT 10
#define TEXT_PADDING 3

class SolutionViewer
{
public:
  SolutionViewer(Grid& grid, int max_pixel_dim);
  void test_draw();
  void display();
  void draw_velocity_magnitude(Grid& grid);
  void draw_status( int iteration, double Re, double lattice_viscosity );
  CImgDisplay window;

private:
  int pixels[2]; // The pixel dimension of the windows.
  int pixels_per_cell; // The number of pixels per coarsest cell.
  CImg<unsigned char> image;
  void scalar2rgb(double min, double max, double value, float rgb[3]);
  void draw_text_line( std::string text, int line );
};