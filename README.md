A C++ LBM code.  

Focus is on incompressible, isothermal high-Re flows. 
The ultimate goal is to simulate multi-element airfoils.

## Assumptions

1. Incompressible and isothermal fluid.  
2. Coarsest grid and time steps are equal to 1.  
3. D2Q9 lattice sites.  

## Main Parameters

1. Mach number: this is the ratio of the boundary condition velocity to the 
  lattice speed of sound, which is 1 / sqrt(3) for grid and time steps of 1. 
  The Mach number should not be set too high for compressibility to affect the 
  results, but increasing Mach number allows for faster convergence and more 
  stability via more viscosity.  
2. Characteristic length: This is not necessarily related to the 
  grid (though in the lid-driven cavity case it is). It is determined by the 
  physics of your problem.   
3. Reynolds number: Perhaps the most important parameter. This is determined by 
  the physics of your problem.  

## Development Notes / TODO

1. Dynamic Grid (DG) (solution-adaptive grid method)  
2. Viscosity counteraction with spatial variation terms.  


## Usage

Run 'make' to compile with visualizer. 
Run 'make NOVIZ=1' to compile without visualizer (it will be slightly faster). 

Modify 'settings' to your needs.  

To run: './lbmpp'. 
For custom maximum window size, run: './lbmpp <maximum resolution dimension>'.


## Validation Cases

These validation cases are compared with data from papers found in the 
'ref' folder.  

### Single-Relaxation Time (SRT) (a.k.a. Standard LBM)

1. Re = 100, Grid = 101x101  
<p align="center">
![SRT Re = 100]
(val/srt_100_u.bmp)
![SRT Re = 100]
(val/srt_100_v.bmp)
![SRT Re = 100]
(val/srt_100_s.bmp)
![SRT Re = 100]
(val/srt_100_m.png)
</p>
2. Re = 1000, Grid = 151x151  
<p align="center">
![SRT Re = 1000]
(val/srt_1000_u.bmp)
![SRT Re = 1000]
(val/srt_1000_v.bmp)
![SRT Re = 1000]
(val/srt_1000_s.bmp)
![SRT Re = 1000]
(val/srt_1000_m.png)
</p>

### Multiple-Relaxation Time (MRT)

1. Re = 100, Grid 101x101    
<p align="center">
![MRT Re = 100]
(val/mrt_100_u.bmp)
![MRT Re = 100]
(val/mrt_100_v.bmp)
![MRT Re = 100]
(val/mrt_100_s.bmp)
![MRT Re = 100]
(val/mrt_100_m.png)
</p>
2. Re = 1000, Grid 151x151  
<p align="center">
![MRT Re = 1000]
(val/mrt_1000_u.bmp)
![MRT Re = 1000]
(val/mrt_1000_v.bmp)
![MRT Re = 1000]
(val/mrt_1000_s.bmp)
![MRT Re = 1000]
(val/mrt_1000_m.png)
</p>
3. Re = 10000, Grid 225x225 (Results could be better by running longer.)  
<p align="center">
![MRT Re = 10000]
(val/mrt_10000_u.bmp)
![MRT Re = 10000]
(val/mrt_10000_v.bmp)
![MRT Re = 10000]
(val/mrt_10000_s.bmp)
![MRT Re = 10000]
(val/mrt_10000_m.png)
</p>
4. Re = 21000:  

### MRT + Viscosity Counteraction (VC)

1. Re = 100:  

2. Re = 1000:  

3. Re = 10000:  

4. Re = 21000:  

### MRT + Dynamic Grid (DG)

1. Re = 100:  

2. Re = 1000:  

3. Re = 10000:  

4. Re = 21000:  


## Stability Limits

### SRT

1. Re = 4000, Grid = 151x151, Unstable.  

## Implementation Details

1. The top corner lattice cells apply stationary walls to 
both of their boundaries, as opposed to one stationary wall 
and one moving wall. If one side is a moving wall, then the 
maximum Reynolds number is severely reduced due to instability.  
2. The lattice links are indexed as follows:  
<p align="center">
![Lattice Link Indices]
(ref/lattice_indices.png)
</p>

## Interesting Notes

1. You cannot use OpenMP on STL iterator for-loops, 
  due to the presence of '!= v.end()'.  




