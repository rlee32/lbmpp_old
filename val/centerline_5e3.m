clear; clc; close all;

addpath subroutines

% This reads centerlines from uniform fields.

% Change the dlmread entry according to the desired solution file.
% The solution file nomencalture descriptions can be found in 
%  the README.txt in the '../results/' folder. 

[x_ref, v_ref] = validation_data_v_vs_x(5000);
[y_ref, u_ref] = validation_data_u_vs_y(5000);

figure(2);
plot(x_ref,v_ref, 'x', 'DisplayName', 'Ghia et al. 257x257');
hold on;
title('Centerline x vs. v-velocity at Re = 5000');
xlabel('x');
ylabel('v');
figure(1);
plot(y_ref,u_ref, 'x', 'DisplayName', 'Ghia et al. 257x257');
hold on;
title('Centerline y vs. u-velocity at Re = 5000');
xlabel('y');
ylabel('u');

% SRT (MRT initialized) 
M = 0.2;
u = dlmread('../results/fields/u_G251_M200_T100_RM1_VCM0_VCF0_Re5000.dat');
v = dlmread('../results/fields/v_G251_M200_T100_RM1_VCM0_VCF0_Re5000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['SRT, M=' num2str(M) ', 251x251'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['SRT, M=' num2str(M) ', 251x251'] );

% MRT 
M = 0.2;
u = dlmread('../results/fields/u_G201_M200_T100_RM3_VCM0_VCF0_Re5000.dat');
v = dlmread('../results/fields/v_G201_M200_T100_RM3_VCM0_VCF0_Re5000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 201x201'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 201x201'] );

% MRT 
M = 0.2;
u = dlmread('../results/fields/u_G251_M200_T100_RM3_VCM0_VCF0_Re5000.dat');
v = dlmread('../results/fields/v_G251_M200_T100_RM3_VCM0_VCF0_Re5000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 251x251'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 251x251'] );

figure(1);
legend('show');
figure(2);
legend('show');
