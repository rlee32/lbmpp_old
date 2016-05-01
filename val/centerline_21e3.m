clear; clc; close all;

addpath subroutines

% This reads centerlines from uniform fields.

% Change the dlmread entry according to the desired solution file.
% The solution file nomencalture descriptions can be found in 
%  the README.txt in the '../results/' folder. 

[x_ref, v_ref] = validation_data_v_vs_x(21000);
[y_ref, u_ref] = validation_data_u_vs_y(21000);

figure(2);
plot(x_ref,v_ref, 'x', 'DisplayName', 'Eturk et al. 601x601');
hold on;
title('Centerline x vs. v-velocity at Re = 21000');
xlabel('x');
ylabel('v');
figure(1);
plot(y_ref,u_ref, 'x', 'DisplayName', 'Eturk et al. 601x601');
hold on;
title('Centerline y vs. u-velocity at Re = 21000');
xlabel('y');
ylabel('u');

% MRT 
M = 0.2;
u = dlmread('../results/fields/u_G257_M200_T500_RM3_VCM0_VCF0_Re21000.dat');
v = dlmread('../results/fields/v_G257_M200_T500_RM3_VCM0_VCF0_Re21000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 257x257'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 257x257'] );

% MRT 325x325
M = 0.2;
u = dlmread('../results/fields/u_G325_M200_T75_RM3_VCM0_VCF0_Re21000.dat');
v = dlmread('../results/fields/v_G325_M200_T75_RM3_VCM0_VCF0_Re21000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 325x325'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 325x325'] );

% MRT 325x325
M = 0.2;
u = dlmread('../results/fields/u_G325_M200_T250_RM3_VCM0_VCF0_Re21000.dat');
v = dlmread('../results/fields/v_G325_M200_T250_RM3_VCM0_VCF0_Re21000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 325x325'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 325x325'] );
% MRT 301x301
M = 0.1;
u = dlmread('../results/fields/u_G301_M100_T300_RM3_VCM0_VCF1_Re21000.dat');
v = dlmread('../results/fields/v_G301_M100_T300_RM3_VCM0_VCF1_Re21000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 301x301'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 301x301'] );
% MRT 301x301
M = 0.3;
u = dlmread('../results/fields/u_G301_M300_T400_RM3_VCM0_VCF0_Re21000.dat');
v = dlmread('../results/fields/v_G301_M300_T400_RM3_VCM0_VCF0_Re21000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 301x301'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 301x301'] );

% MRT 375x375
M = 0.2;
u = dlmread('../results/fields/u_G375_M200_T151_RM3_VCM0_VCF0_Re21000.dat');
v = dlmread('../results/fields/v_G375_M200_T151_RM3_VCM0_VCF0_Re21000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 375x375'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 375x375'] );

% MRT 401x401
M = 0.2;
u = dlmread('../results/fields/u_G401_M200_T125_RM3_VCM0_VCF1_Re21000.dat');
v = dlmread('../results/fields/v_G401_M200_T125_RM3_VCM0_VCF1_Re21000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 401x401'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 401x401'] );

% MRT 601x601
M = 0.2;
u = dlmread('../results/fields/u_G601_M200_T30_RM3_VCM0_VCF0_Re21000.dat');
v = dlmread('../results/fields/v_G601_M200_T30_RM3_VCM0_VCF0_Re21000.dat');
U = M / sqrt(3);
figure(1);
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 601x601'] );
figure(2);
[x, vc] = get_centerline_v(v);
plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 601x601'] );

figure(1);
legend('show');
figure(2);
legend('show');
