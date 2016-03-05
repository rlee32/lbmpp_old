#include "SolutionViewer.h"

using namespace std;

SolutionViewer::SolutionViewer(Grid& grid, int max_pixel_dim)
{
  int max_cell_dim = (grid.cell_count[0] > grid.cell_count[1]) ? grid.cell_count[0] : grid.cell_count[1];
  double pixels_per_cell = (double) max_pixel_dim / (double) max_cell_dim;
  pixels_per_cell = ( pixels_per_cell < 1 ) ? 1 : pixels_per_cell;
  // double aspect_ratio = grid.cell_count[0] / grid.cell_count[1];

  cout << pixels_per_cell << endl;
  cout << grid.cell_count[0] << endl;
  cout << grid.cell_count[1] << endl;
  cout << pixels_per_cell*grid.cell_count[0] << endl;
  cout << pixels_per_cell*grid.cell_count[1] << endl;

  CImg<unsigned char> image_(
    pixels_per_cell*grid.cell_count[0], // x-resolution
    pixels_per_cell*grid.cell_count[1], // y-resolution
    1, // 1 z-layer (2D)
    3, // 3 channels (RGB)
    255); // default background color.
  image = image_;
  CImgDisplay window_(image, "Flow solution");
  window = window_;
}

void SolutionViewer::test_draw()
{
  const float black[] = {0,0,0};
  image.draw_line(1,1,50,50,black);
  image.draw_line(100,100,50,50,black);
  image.draw_fill(100,100,black);
  image.draw_rectangle(100,100,150,150,black);
  image.draw_grid(10,10,0,0,false,false,black);
}

void SolutionViewer::draw_solution(Grid& grid)
{
  // image.draw_grid(10,10,0,0,false,false,black);
}

void SolutionViewer::display()
{
  image.display(window);
}

void SolutionViewer::scalar2rgb(double min, double max, double value, int* rgb)
{
  // Derived from: https://www.particleincell.com/2014/colormap/
  // Normalize
  double f = (value - min) / (max - min);
  // Invert the value to get blue to red as 0 to 1.
  double a = 4 * ( 1 - f ); //invert and group
  int X = (int) a;  //this is the integer part
  int Y = (int) ( 255 * ( a - X ) ); //fractional part from 0 to 255
  rgb[0] = 0;
  rgb[1] = 0;
  rgb[2] = 0;
  switch(X)
  {
      case 0: 
        rgb[0]=255; 
        rgb[1]=Y; 
        break;
      case 1: 
        rgb[0]=255-Y; 
        rgb[1]=255; 
        break;
      case 2: 
        rgb[1]=255; 
        rgb[2]=Y; 
        break;
      case 3: 
        rgb[1]=255-Y; 
        rgb[2]=255; 
        break;
      case 4: 
        rgb[2]=255; 
        break;
  }
}
