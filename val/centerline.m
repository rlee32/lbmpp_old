clear; close all; clc;

% Change the dlmread entry according to the desired solution file.
% The solution file nomencalture descriptions can be found in 
%  the README.txt in the '../results/' folder. 

% Tall cavity
% M = 0.2;
% H = 1.5;
% Re = 100;
% u = dlmread('../results/u_G128x192_M200_T100_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G128x192_M200_T100_RM1_VCM0_VCF0_Re100.dat');

% % SRT 100 Re
% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G200_M200_T50_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G200_M200_T50_RM1_VCM0_VCF0_Re100.dat');

% SRT 100 Re
M = 0.2;
H = 1;
Re = 100;
u = dlmread('../results/u_G200_M200_T20_RM3_VCM0_VCF0_Re100.dat');
v = dlmread('../results/v_G200_M200_T20_RM3_VCM0_VCF0_Re100.dat');

% % SRT 100 Re Low-res
% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G32_M200_T50_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G32_M200_T50_RM1_VCM0_VCF0_Re100.dat');

% MRT 10,000 Re (unconverged)
% M = 0.2;
% H = 1;
% Re = 10000;
% u = dlmread('../results/u_G200_M200_T500_RM3_VCM0_VCF0_Re10000.dat');
% v = dlmread('../results/v_G200_M200_T500_RM3_VCM0_VCF0_Re10000.dat');

% M = 0.2;
% H = 1.5;
% Re = 500;
% u = dlmread('../results/u_G128x192_M200_T100_RM1_VCM0_VCF0_Re500.dat');
% v = dlmread('../results/v_G128x192_M200_T100_RM1_VCM0_VCF0_Re500.dat');

% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G128_M200_T20_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G128_M200_T20_RM1_VCM0_VCF0_Re100.dat');

% M = 0.2;
% H = 1;
% Re = 500;
% u = dlmread('../results/u_G128_M200_T75_RM1_VCM0_VCF0_Re500.dat');
% v = dlmread('../results/v_G128_M200_T75_RM1_VCM0_VCF0_Re500.dat');

% M = 0.004;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G128_M4_T100_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G128_M4_T100_RM1_VCM0_VCF0_Re100.dat');

% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G50_M200_T20_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G50_M200_T20_RM1_VCM0_VCF0_Re100.dat');

% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G25_M200_T20_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G25_M200_T20_RM1_VCM0_VCF0_Re100.dat');

% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G75_M200_T20_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G75_M200_T20_RM1_VCM0_VCF0_Re100.dat');
% 
% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G200_M200_T20_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G200_M200_T20_RM1_VCM0_VCF0_Re100.dat');
% 
% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G200_M200_T30_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G200_M200_T30_RM1_VCM0_VCF0_Re100.dat');

% M = 0.1;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G128_M100_T40_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G128_M100_T40_RM1_VCM0_VCF0_Re100.dat');

% M = 0.05;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G128_M50_T80_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G128_M50_T80_RM1_VCM0_VCF0_Re100.dat');

% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G100_M200_T20_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G100_M200_T20_RM1_VCM0_VCF0_Re100.dat');

% M = 0.2;
% H = 1;
% Re = 100;
% u = dlmread('../results/u_G150_M200_T40_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G150_M200_T40_RM1_VCM0_VCF0_Re100.dat');




[rows,cols] = size(u);
rows_even = mod(rows, 2) == 0;
cols_even = mod(cols, 2) == 0;

dx = 1 / cols;
dy = 1 / rows;
% x = ( linspace(0,1-dx,cols) ) / ( 1 - dx );
% y = ( linspace(0,1-dy,rows) ) / ( 1 - dy );
x = ( linspace(dx/2,1-dx/2,cols) );
y = ( linspace(dy/2,H-dy/2,rows) );

if rows_even
    rows = [rows/2, rows/2+1];
else
    rows = ceil(rows/2);
end
if cols_even
    cols = [cols/2, cols/2+1];
else
    cols = ceil(cols/2);
end

validation = H == 1 && (Re==100 || Re==10000);

if validation
    [x_ref, v_ref] = validation_data_v_vs_x(Re);
    [y_ref, u_ref] = validation_data_u_vs_y(Re);
end


figure;
if validation
    plot(y_ref,u_ref, 'x');
    hold on;
end
U = M / sqrt(3);
plot( y, (mean(u(:,cols),2)/U) );
title(['Centerline y vs. u-velocity at Re = ' num2str(Re)]);
xlabel('y');
ylabel('u');
if validation
    legend('Ghia et al', 'Present LBM');
end

figure;
if validation
    plot(x_ref,v_ref, 'x');
    hold on;
end
plot( x, mean(v(rows,:),1)/U );
title(['Centerline x vs. v-velocity at Re = ' num2str(Re)]);
xlabel('x');
ylabel('v');
if validation
    legend('Ghia et al', 'Present LBM');
end

figure;
[X,Y] = meshgrid(x,y);
axis equal tight;
h = streamslice( X, Y, rot90(-u,2), rot90(-v,2) );
rotate(h, [0,0,1], 180);
% streamline( stream2(X,Y,u,v,x,y) );
% streamline( X,Y,u,v,x,y );
title(['Steady-State Streamlines at Re = ' num2str(Re)]);
xlabel('x');
ylabel('y');

ee = reference_error(x,v,Re);
disp(['error ' num2str(ee)]);
