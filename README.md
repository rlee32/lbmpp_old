A C++ LBM code.  

Focus is on incompressible, isothermal high-Re flows. 
The goal is to simulate multi-element airfoils.

## Assumptions

1. Incompressible and isothermal fluid.  
2. Coarsest grid and time steps are equal to 1.  

## Development Notes / TODO

1. Make tables for each grid level for quantities that are the same within a grid level.  



## Usage

Run 'make' to compile with visualizer. 
Run 'make NOVIZ=1' to compile without visualizer (it will be slightly faster). 

Modify 'settings' to your needs.  

To run: './lbmpp'. 
For custom maximum window size, run: './lbmpp <maximum resolution dimension>'.


## Validation Cases

These validation cases are compared with data from papers found in the 
'ref' folder.  

### SRT

1. Re = 100:  

![SRT Re = 100]
(results/srt_100_m.png)


2. Re = 1000:  

### MRT

1. Re = 100:  

2. Re = 1000:  

3. Re = 10000:  

4. Re = 21000:  

### MRT + VC

1. Re = 100:  

2. Re = 1000:  

3. Re = 10000:  

4. Re = 21000:  

### MRT + DG (Dynamic Grid)

1. Re = 100:  

2. Re = 1000:  

3. Re = 10000:  

4. Re = 21000:  

## Interesting Notes

You cannot use OpenMP on STL iterator for-loops, due to the presence of '!= v.end()'. 




