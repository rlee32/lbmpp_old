#include "SolutionViewer.h"

using namespace std;

SolutionViewer::SolutionViewer(Simulator& sim, int max_pixel_dim) : 
  last_elapsed_time(0), last_iteration(0)
{
  temp.min = 0;
  temp.max = 0;
  int max_cell_dim = (sim.get_cell_count_0() > sim.get_cell_count_1()) ? 
    sim.get_cell_count_0() : sim.get_cell_count_1();
  pixels_per_cell = (double) max_pixel_dim / (double) max_cell_dim;
  pixels_per_cell = ( pixels_per_cell < 1 ) ? 1 : pixels_per_cell;
  pixels_per_cell = (int) pixels_per_cell;
  pixels[0] = pixels_per_cell*sim.get_cell_count_0();
  pixels[1] = pixels_per_cell*sim.get_cell_count_1();
  CImg<unsigned char> image_(
    pixels[0], // x-resolution
    pixels[1]+TextDisplayDim, // y-resolution
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
  size_t cell_count_x = grid.cell_count_x();
  size_t cell_count_y = grid.cell_count_y();
  grid.reconstruct_macro();
  temp.min = grid.min_mag();
  temp.max = grid.max_mag();
  for (size_t i = 0; i < cell_count_x; ++i)
  {
    for (size_t j = 0; j < cell_count_y; ++j)
    {
      size_t ii = i + j * cell_count_x;
      size_t jp = cell_count_y - 1 - j;
      draw_mag_tree( &(grid.get_level(1)), (grid.get_cells(0))[ii],
        i*pixels_per_cell, jp*pixels_per_cell, pixels_per_cell );
    }
  }
  if (pixels_per_cell > 2)
  {
    const float black[] = {0,0,0};
    image.draw_grid(pixels_per_cell,pixels_per_cell,0,0,false,false,black);
  }
}

// cg: grid level of the children of this cell.
// cell: current cell whose children will be drawn if not active.
// i: the upper-left corner horizontal position in pixels.
// j: the upper-left corner vertical position in pixels.
// p: pixel dimension of this cell.
void SolutionViewer::draw_mag_tree( GridLevel* cg, 
  Cell& cell, size_t i, size_t j, size_t p )
{
  if (not cell.state.active and p > 2 and cg != nullptr )
  {
    // draw children
    size_t cdim = p / 2;
    vector<Cell>& g = cg->get_cells();
    draw_mag_tree( cg->get_next_grid_level(), 
      g[ cell.local.children[0] ], i, j+cdim, cdim );
    draw_mag_tree( cg->get_next_grid_level(), 
       g[ cell.local.children[1] ], i, j, cdim );
    draw_mag_tree( cg->get_next_grid_level(), 
       g[ cell.local.children[2] ], i+cdim, j+cdim, cdim );
    draw_mag_tree( cg->get_next_grid_level(), 
       g[ cell.local.children[3] ], i+cdim, j, cdim );
    if (pixels_per_cell > 4)
    {
      const float black[] = {0,0,0};
      // image.draw_line(i,j,i+p,j,black);
      // image.draw_line(i+p,j,i+p,j+p,black);
      // image.draw_line(i+p,j+p,i,j+p,black);
      // image.draw_line(i,j+p,i,j,black);
      image.draw_line(i,j+p/2,i+p,j+p/2,black);
      image.draw_line(i+p/2,j,i+p/2,j+p,black);
    }
  }
  else
  {
    // cell.reconstruct_macro();
    // Draw this cell
    float rgb[3];
    scalar2rgb( temp.min, temp.max, cell.get_mag(), rgb );
    image.draw_rectangle(
      i, // upper left corner
      j,
      i+p, // lower right corner
      j+p,
      rgb );
  }
}

void SolutionViewer::draw_status( int iteration, Simulator& sim,  
  double elapsed_time )
{
  const float white[] = { 255, 255, 255 };
  image.draw_rectangle(
    1, // upper left corner
    pixels[1]+1,
    pixels[0], // lower right corner
    pixels[1]+TextDisplayDim,
    white );
  double dt = elapsed_time - last_elapsed_time;
  double di = iteration - last_iteration;
  double speed = (dt > 0) ? di / dt : 0;
  double elapsed_time_seconds = fmod(elapsed_time,60);
  size_t elapsed_time_minutes = (size_t)( elapsed_time / 60 );
  size_t remaining_iter = sim.get_timesteps() - iteration;
  double avg_speed = (double)iteration / elapsed_time;
  double remaining_seconds_total = remaining_iter / avg_speed;
  double remaining_seconds = fmod(remaining_seconds_total,60);
  size_t remaining_minutes = (size_t)( remaining_seconds_total / 60 );
  string text = "";
  uint line = 0;
  text = "Iteration: " + to_string(iteration) + " / " 
    + to_string((size_t) sim.get_timesteps());
  draw_text_line( text, line++ );
  text = "Elapsed time: " + to_string(elapsed_time_minutes) + " min, " 
    + to_string(elapsed_time_seconds) + " sec"
    + " ( Estimated " + to_string(remaining_minutes) + " minutes, " 
    + to_string(remaining_seconds) + " sec remaining )";
  draw_text_line( text, line++ );
  text = "Total cells: " + to_string(sim.grid.active_cells());
  draw_text_line( text, line++ );
  text = "Computation speed: " + to_string(speed) + " timesteps / second"
    + " ( average speed: " + to_string(avg_speed) + " timesteps / second )";
  draw_text_line( text, line++ );
  text = "Reynolds number: " + to_string(sim.get_Re());
  draw_text_line( text, line++ );
  text = "Relaxation time: " + to_string(sim.get_tau());
  draw_text_line( text, line++ );
  text = "Min, max density: " + to_string(sim.grid.min_rho()) + ", "
    + to_string(sim.grid.max_rho());
  draw_text_line( text, line++ );
  text = "Min, max velocity magnitude (BC: " + to_string( sim.get_U() ) 
    + "): " + to_string( temp.min ) + ", "
    + to_string( temp.max );
  draw_text_line( text, line++ );
  last_elapsed_time = elapsed_time;
  last_iteration = iteration;
}

void SolutionViewer::draw_text_line( string text, int line )
{
  unsigned char black[] = { 0, 0, 0 };
  int pixel_increment = TextHeight + TextPadding;
  image.draw_text(
    1 + TextPadding, // upper left corner
    1 + TextPadding + pixels[1] + pixel_increment * line, 
    text.c_str(), // string to display
    black ); // foreground color
}

void SolutionViewer::display()
{
  image.display(window);
}

void SolutionViewer::scalar2rgb(
  double min, double max, double value, float rgb[3])
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
