clear; clc; close all;

addpath subroutines

% This reads centerlines from uniform fields.

% Change the dlmread entry according to the desired solution file.
% The solution file nomencalture descriptions can be found in 
%  the README.txt in the '../results/' folder. 

[x_ref, v_ref] = validation_data_v_vs_x(10000);
[y_ref, u_ref] = validation_data_u_vs_y(10000);

figure(2);
plot(x_ref,v_ref, 'x', 'DisplayName', 'Ghia et al. 257x257');
hold on;
title('Centerline x vs. v-velocity at Re = 10000');
xlabel('x');
ylabel('v');
figure(1);
plot(y_ref,u_ref, 'x', 'DisplayName', 'Ghia et al. 257x257');
hold on;
title('Centerline y vs. u-velocity at Re = 10000');
xlabel('y');
ylabel('u');

% % MRT 225x225
% M = 0.2;
% u = dlmread('../results/fields/u_G225_M200_T1000_RM3_VCM0_VCF0_Re10000.dat');
% v = dlmread('../results/fields/v_G225_M200_T1000_RM3_VCM0_VCF0_Re10000.dat');
% U = M / sqrt(3);
% figure(1);
% [y, uc] = get_centerline_u(u,U);
% plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 225x225'] );
% figure(2);
% [x, vc] = get_centerline_v(v);
% plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 225x225'] );

% MRT 257x257
M = 0.2;
u = dlmread('../results/fields/u_G257_M200_T200_RM3_VCM0_VCF0_Re10000.dat');
v = dlmread('../results/fields/v_G257_M200_T200_RM3_VCM0_VCF0_Re10000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 257x257'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 257x257'] );

% % MRT 257x257
% M = 0.1;
% u = dlmread('../results/fields/u_G257_M100_T200_RM3_VCM0_VCF0_Re10000.dat');
% v = dlmread('../results/fields/v_G257_M100_T200_RM3_VCM0_VCF0_Re10000.dat');
% U = M / sqrt(3);
% figure(1);
% [y, uc] = get_centerline_u(u,U);
% plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 257x257'] );
% figure(2);
% [x, vc] = get_centerline_v(v);
% plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 257x257'] );

figure(1);
legend('show');
figure(2);
legend('show');
