These results are for lid-driven cavity flow.  

Each run has 2 files: u and v for the velocity components.  

The name format is:  <component>_G<grid dimension>_M<milliMach number>_T<thousands of timesteps>_RM<relaxation model>_VCM<VC model>_VCF<nucf tenths>_Re<Reynolds number>.dat  

Each component is preceded by a tag to indicate the field value.  

milliMach number is thousandths of a Mach number.  

The relaxation model signifies the following:  
1: single-relaxation time  
2: two-relaxation time  
3: multiple-relaxation time  

The VC model signifies the following:  
0: No viscosity counteraction.  
1: No spatial or temporal variation approximation.  
2: No temporal variation approximation.  
3: Include both temporal and spatial effects. Inclusion of temporal makes for an implicit method.  

nucf is the ratio of additional buffer viscosity nuc to the standard 
viscosity, as tenths. So a value of 1 would corresond to nucf = 0.1  

So the u-velocity component of a 128 x 128 grid at Mach 0.005 run for 250,000 timesteps would be:  
u_128_5_250.dat  

The v-velocity component of a 150 x 150 grid at Mach 0.01 run for 250,000 timesteps with Multiple-Relaxation time would be:  
v_150_10_250_3.dat  

The u-velocity component of a 100 x 100 grid at Mach 0.006 run for 250,000 timesteps with Multiple-Relaxation time and visosity counteraction local 
approximation would be:  
u_100_6_250_3_1.dat  
