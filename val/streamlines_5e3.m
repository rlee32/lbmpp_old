clear; clc;
close all;

addpath subroutines

% This currently only works for regular grids.
% For irregular grids, griddata might be the way to go.

% Change the dlmread entry according to the desired solution file.
% The solution file nomencalture descriptions can be found in 
%  the README.txt in the '../results/' folder. 

% MRT 
M = 0.2;
Re = 5e3;
u = dlmread('../results/fields/u_G251_M200_T100_RM3_VCM0_VCF0_Re5000.dat');
v = dlmread('../results/fields/v_G251_M200_T100_RM3_VCM0_VCF0_Re5000.dat');

[rows,cols] = size(u);
dx = 1 / cols;
dy = 1 / rows;
H = 1; % height of cavity.
x = ( linspace(dx/2,1-dx/2,cols) );
y = ( linspace(dy/2,H-dy/2,rows) );

figure;
[X,Y] = meshgrid(x,y);
axis equal tight;
h = streamslice( X, Y, rot90(-u,2), rot90(-v,2), 2 );
rotate(h, [0,0,1], 180);
title(['Steady-State Streamlines at Re = ' num2str(Re) ...
    ', M = ' num2str(M)]);
xlabel('x');
ylabel('y');
