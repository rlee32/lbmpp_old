#include "SolutionViewer.h"

using namespace std;

SolutionViewer::SolutionViewer(Grid& grid, int max_pixel_dim)
{
  int max_cell_dim = (grid.cell_count[0] > grid.cell_count[1]) ? grid.cell_count[0] : grid.cell_count[1];
  pixels_per_cell = (double) max_pixel_dim / (double) max_cell_dim;
  pixels_per_cell = ( pixels_per_cell < 1 ) ? 1 : pixels_per_cell;
  pixels_per_cell = (int) pixels_per_cell;
  pixels[0] = pixels_per_cell*grid.cell_count[0];
  pixels[1] = pixels_per_cell*grid.cell_count[1];
  CImg<unsigned char> image_(
    pixels[0], // x-resolution
    pixels[1], // y-resolution
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

void SolutionViewer::draw_velocity_magnitude(Grid& grid)
{
  const float black[] = {0,0,0};
  double min = grid.get_min_velocity_magnitude();
  double max = grid.get_max_velocity_magnitude();
  cout << "Velocity magnitude min, max: " << min << ", " << max << endl;
  for (int i = 0; i < grid.cell_count[0]; ++i)
  {
    for (int j = 0; j < grid.cell_count[1]; ++j)
    {    
      // Get color
      float rgb[3];
      scalar2rgb(min, max, grid.grid_levels[0][i+j*grid.cell_count[0]].get_velocity_magnitude(), rgb);
      // cout << rgb[0] << ", " << rgb[1] << ", " << rgb[2] << endl;
      // Draw
      image.draw_rectangle(
        i*pixels_per_cell, // upper left corner
        j*pixels_per_cell,
        (i+1)*pixels_per_cell, // lower right corner
        (j+1)*pixels_per_cell,
        rgb );
    }
  }
  if (pixels_per_cell > 2)
  {
    image.draw_grid(pixels_per_cell,pixels_per_cell,0,0,false,false,black);
  }
}

void SolutionViewer::display()
{
  image.display(window);
}

void SolutionViewer::scalar2rgb(double min, double max, double value, float rgb[3])
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
  if (max == min) X = 4;
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
