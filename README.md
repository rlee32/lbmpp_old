A C++ LBM code.  

Focus is on incompressible, isothermal high-Re flows. The goal is to simulate multi-element airfoils.

## Assumptions

1. Incompressible and isothermal fluid.  
2. Coarsest grid and time steps are equal to 1.  

## Development Notes / TODO

1. Make tables for each grid level for quantities that are the same within a grid level.  



## Usage

Run 'make' to compile with visualizer. Run 'make NOVIZ=1' to compile without visualizer (it will be faster). 

Modify 'settings' to your needs.  

To run: './lbmpp'. For custom maximum window size, run: './lbmpp <maximum resolution dimension>'.

## Interesting Notes

You cannot use OpenMP on STL iterator for-loops, due to the presence of '!= v.end()'. 




