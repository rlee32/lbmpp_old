clear; clc;
close all;

addpath subroutines

% This currently only works for regular grids.
% For irregular grids, griddata might be the way to go.

% Change the dlmread entry according to the desired solution file.
% The solution file nomencalture descriptions can be found in 
%  the README.txt in the '../results/' folder. 

% % SRT 100 Re
% M = 0.1;
% Re = 100;
% u = dlmread('../results/u_G101_M100_T100_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G101_M100_T100_RM1_VCM0_VCF0_Re100.dat');

% % SRT 1000 Re
% M = 0.1;
% Re = 1000;
% u = dlmread('../results/u_G151_M100_T125_RM1_VCM0_VCF0_Re1000.dat');
% v = dlmread('../results/v_G151_M100_T125_RM1_VCM0_VCF0_Re1000.dat');

% % SRT 2500 Re
% M = 0.1;
% Re = 2500;
% u = dlmread('../results/u_G151_M100_T125_RM1_VCM0_VCF0_Re2500.dat');
% v = dlmread('../results/v_G151_M100_T125_RM1_VCM0_VCF0_Re2500.dat');

% MRT 100 Re
M = 0.1;
Re = 100;
u = dlmread('../results/fields/u_G101_M100_T50_RM3_VCM0_VCF0_Re100.dat');
v = dlmread('../results/fields/v_G101_M100_T50_RM3_VCM0_VCF0_Re100.dat');

% % MRT 1000 Re
% M = 0.1;
% Re = 1000;
% % u = dlmread('../results/u_G100_M200_T50_RM3_VCM0_VCF0_Re1000.dat');
% % v = dlmread('../results/v_G100_M200_T50_RM3_VCM0_VCF0_Re1000.dat');
% u = dlmread('../results/fields/u_G151_M100_T100_RM3_VCM0_VCF0_Re1000.dat');
% v = dlmread('../results/fields/v_G151_M100_T100_RM3_VCM0_VCF0_Re1000.dat');

% % MRT 5000 Re
% M = 0.2;
% Re = 5000;
% u = dlmread('../results/u_G101_M200_T500_RM3_VCM0_VCF0_Re5000.dat');
% v = dlmread('../results/v_G101_M200_T500_RM3_VCM0_VCF0_Re5000.dat');

% % MRT 10000 Re
% M = 0.2;
% Re = 10000;
% u = dlmread('../results/u_G257_M200_T200_RM3_VCM0_VCF0_Re10000.dat');
% v = dlmread('../results/v_G257_M200_T200_RM3_VCM0_VCF0_Re10000.dat');

% % MRT 20000 Re
% M = 0.1;
% Re = 20000;
% u = dlmread('../results/u_G125_M200_T500_RM3_VCM0_VCF0_Re20000.dat');
% v = dlmread('../results/v_G125_M200_T500_RM3_VCM0_VCF0_Re20000.dat');

% % MRT 21000 Re
% M = 0.2;
% Re = 21000;
% u = dlmread('../results/u_G125_M200_T500_RM3_VCM0_VCF0_Re21000.dat');
% v = dlmread('../results/v_G125_M200_T500_RM3_VCM0_VCF0_Re21000.dat');

[rows,cols] = size(u);
dx = 1 / cols;
dy = 1 / rows;
% x = ( linspace(0,1-dx,cols) ) / ( 1 - dx );
% y = ( linspace(0,1-dy,rows) ) / ( 1 - dy );
H = 1; % height of cavity.
x = ( linspace(dx/2,1-dx/2,cols) );
y = ( linspace(dy/2,H-dy/2,rows) );

figure;
[X,Y] = meshgrid(x,y);
axis equal tight;
h = streamslice( X, Y, rot90(-u,2), rot90(-v,2), 2 );
rotate(h, [0,0,1], 180);
% streamline( stream2(X,Y,u,v,x,y) );
% streamline( X,Y,u,v,x,y );
title(['Steady-State Streamlines at Re = ' num2str(Re) ...
    ', M = ' num2str(M)]);
xlabel('x');
ylabel('y');
