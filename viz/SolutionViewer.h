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
  SolutionViewer(Simulator& sim, int max_pixel_dim);
  void test_draw();
  void display();
  void draw_velocity_magnitude(Grid& grid);
  void draw_status( 
    int iteration, Simulator& sim, double elapsed_time );
  void save_image(std::string filename) { image.save( filename.c_str() ); }
  CImgDisplay window;
private:
  static const uint TextHeight = 10;
  static const uint TextPadding = 3;
  static const uint TextDisplayDim = 8*( TextHeight + TextPadding ) 
    + TextPadding;
  std::size_t pixels[2]; // The pixel dimension of the windows.
  std::size_t pixels_per_cell; // The number of pixels per coarsest cell.
  CImg<unsigned char> image;
  void scalar2rgb(double min, double max, double value, float rgb[3]);
  void draw_text_line( std::string text, int line );
  double last_elapsed_time;
  double last_iteration;
  // Temp values used during drawing field solution.
  struct
  {
    double min;
    double max;
  } temp;
  void draw_mag_tree( GridLevel* cg, 
    Cell& cell, std::size_t i, std::size_t j, std::size_t p );
};