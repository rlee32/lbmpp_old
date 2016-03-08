A C++ LBM code.  

Focus is on incompressible, isothermal high-Re flows. Made with the target of simulating multi-element airfoils.  

## Development Notes / TODO

1. Make tables for each grid level for quantities that are the same within a grid level.  



## Usage

Run 'make' to compile with visualizer. Run 'make NOVIZ=1' to compile without visualizer (it will be faster). 

Modify 'settings' to your needs.  

To run: './lbmpp'. For custom maximum window size, run: './lbmpp <maximum resolution dimension>'.

