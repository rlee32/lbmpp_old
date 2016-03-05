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
  void draw_solution(Grid& grid);
  CImgDisplay window;
private:
  CImg<unsigned char> image;
  void scalar2rgb(double min, double max, double value, int* rgb);
};