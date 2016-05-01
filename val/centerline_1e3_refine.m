clear; clc; close all;

addpath subroutines

% This reads centerlines from uniform fields.

% Change the dlmread entry according to the desired solution file.
% The solution file nomencalture descriptions can be found in 
%  the README.txt in the '../results/' folder. 

[x_ref, v_ref] = validation_data_v_vs_x(1000);

[y_ref, u_ref] = validation_data_u_vs_y(1000);

figure(2);
plot(x_ref,v_ref, 'x', 'DisplayName', 'Ghia et al. 129x129');
hold on;
title('Centerline x vs. v-velocity at Re = 1000');
xlabel('x');
ylabel('v');
figure(1);
plot(y_ref,u_ref, 'x', 'DisplayName', 'Ghia et al. 129x129');
hold on;
title('Centerline y vs. u-velocity at Re = 1000');
xlabel('y');
ylabel('u');

% MRT 150x150
M = 0.2;
u = dlmread('../results/fields/u_G150_M200_T100_RM3_VCM0_VCF0_Re1000.dat');
v = dlmread('../results/fields/v_G150_M200_T100_RM3_VCM0_VCF0_Re1000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 150x150'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 150x150'] );

% MRT 75 to 150
% Refined
M = 0.2;
[x, y, uc, vc] = get_centerlines( ...
    '../results/centerlines/centerlines_G75_M200_T50_RM3_VCM0_VCF0_Re1000.tsv' );
% [x, y, uc, vc] = get_centerlines( ...
%     '../results/centerlines/centerlines_G151_M100_T100_RM1_VCM0_VCF0_Re1000.tsv' );

U = M / sqrt(3);
figure(1);
plot( y, uc/U, 'DisplayName', ...
    ['MRT, M=' num2str(M) ', 75x75 refined to 150x150'] );
figure(2);
plot( x, vc/U, 'DisplayName', ...
    ['MRT, M=' num2str(M) ', 75x75 refined to 150x150'] );

figure(1);
legend('show');
figure(2);
legend('show');
