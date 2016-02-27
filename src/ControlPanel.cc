#include "ControlPanel.h"

using namespace std;

ControlPanel::ControlPanel()
{
  
}

void ControlPanel::read_settings(string filename)
{
  ifstream settings_file( filename );
  string line;
  if(settings_file.is_open())
  {
    while( getline(settings_file, line) )
    {
      istringstream iss(line);
      string parameter;
      if( (iss >> parameter) )
      {
        if ( parameter == "coarse_cell_dim" ) iss >> coarse_cell_dim;
      }
    }
  }
}