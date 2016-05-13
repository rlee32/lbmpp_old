A C++ LBM code.  

Focus is on incompressible, isothermal high-Re flows. 
The ultimate goal is solution-adaptive rotor-stator simulations.

## Usage

Run 'make' to compile with visualizer. 
Run 'make NOVIZ=1' to compile without visualizer (it will be slightly faster). 

Modify 'settings' to your needs.  

To run: './lbmpp'.  
For custom maximum window size, run: './lbmpp \<maximum resolution dimension\>'.  

## Parameters

The following are the main parameters:  

1. Mach number: this is the ratio of the boundary condition velocity to the 
  lattice speed of sound, which is 1 / sqrt(3) for grid and time steps of 1. 
  The Mach number should not be set too high for compressibility to affect the 
  results, but increasing Mach number allows for less timesteps to convergence.  
2. Characteristic length: This is not necessarily related to the 
  grid (though in the lid-driven cavity case it is). It is determined by the 
  physics of your problem.   
3. Reynolds number: This is determined by the physics of your problem.  

The complete list of parameters can be found in 'settings'.  

## Assumptions, Implementation Details, and Validation

See 'ref/report.pdf'.  

## Development Notes / TODO

1. Characterize spatial order of accuracy of grid interfaces.  
2. Solution-adaptive refinement.  

<!--
## Assumptions

1. Incompressible and isothermal fluid.  
2. Coarsest grid and time steps are equal to 1.  
3. D2Q9 lattice sites.  
-->

<!--
### Single-Relaxation Time (SRT) (a.k.a. Standard LBM)

1. Re = 100, Grid = 101x101  
<p align="center">
![SRT Re = 100]
(val/pics/srt_100_u.bmp)
![SRT Re = 100]
(val/pics/srt_100_v.bmp)
![SRT Re = 100]
(val/pics/srt_100_s.bmp)
![SRT Re = 100]
(val/pics/srt_100_m.png)
</p>
2. Re = 1000, Grid = 151x151  
<p align="center">
![SRT Re = 1000]
(val/pics/srt_1000_u.bmp)
![SRT Re = 1000]
(val/pics/srt_1000_v.bmp)
![SRT Re = 1000]
(val/pics/srt_1000_s.bmp)
![SRT Re = 1000]
(val/pics/srt_1000_m.png)
</p>

### Multiple-Relaxation Time (MRT)

1. Re = 100, Grid 101x101    
<p align="center">
![MRT Re = 100]
(val/pics/mrt_100_u.bmp)
![MRT Re = 100]
(val/pics/mrt_100_v.bmp)
![MRT Re = 100]
(val/pics/mrt_100_s.bmp)
![MRT Re = 100]
(val/pics/mrt_100_m.png)
</p>
2. Re = 1000, Grid 151x151  
<p align="center">
![MRT Re = 1000]
(val/pics/mrt_1000_u.bmp)
![MRT Re = 1000]
(val/pics/mrt_1000_v.bmp)
![MRT Re = 1000]
(val/pics/mrt_1000_s.bmp)
![MRT Re = 1000]
(val/pics/mrt_1000_m.png)
</p>
3. Re = 10000, Grid 257x257  
<p align="center">
![MRT Re = 10000]
(val/pics/mrt_10000_u.bmp)
![MRT Re = 10000]
(val/pics/mrt_10000_v.bmp)
![MRT Re = 10000]
(val/pics/mrt_10000_s.bmp)
![MRT Re = 10000]
(val/pics/mrt_10000_m.png)
</p>
-->

<!--
### MRT + Viscosity Counteraction (VC)

1. Re = 100:  

2. Re = 1000:  

3. Re = 10000:  

### MRT + Dynamic Grid (DG)

1. Re = 100:  

2. Re = 1000:  

3. Re = 10000:   
-->

<!--
## Stability Limits

### SRT

1. Re = 4000, Grid = 151x151, Unstable.  
-->

<!--
## Implementation Details

See 'ref/report.pdf'.
-->

<!--
## Miscellaneous Interesting Notes

1. You cannot use OpenMP on STL iterator for-loops, 
  due to the presence of '!= v.end()'.  
-->




