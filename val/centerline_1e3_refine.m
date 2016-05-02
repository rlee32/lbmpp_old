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

% % MRT 100x100
% M = 0.2;
% u = dlmread('../results/fields/u_G100_M200_T50_RM3_VCM0_VCF0_Re1000.dat');
% v = dlmread('../results/fields/v_G100_M200_T50_RM3_VCM0_VCF0_Re1000.dat');
% U = M / sqrt(3);
% figure(1);
% [y, uc] = get_centerline_u(u,U);
% plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 100x100'] );
% figure(2);
% [x, vc] = get_centerline_v(v);
% plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 100x100'] );

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

% MRT 50 / 100
% Refined
M = 0.2;
[x, y, uc, vc] = get_centerlines( ...
    '../results/centerlines/centerlines_G50_M200_T50_RM3_VCM0_VCF0_Re1000.tsv');
U = M / sqrt(3);
figure(1);
plot( y, uc/U, 'DisplayName', ...
    ['MRT, M=' num2str(M) ', 50x25 + 100x50'] );
figure(2);
plot( x, vc/U, 'DisplayName', ...
    ['MRT, M=' num2str(M) ', 50x25 + 100x50'] );

% MRT 50 / 100 / 200
% Refined
M = 0.2;
[x, y, uc, vc] = get_centerlines( ...
    '../results/centerlines/centerlines_G50_M200_T50_RM3_VCM0_VCF0_Re1000_refine.tsv');
U = M / sqrt(3);
figure(1);
plot( y, uc/U, 'DisplayName', ...
    ['MRT, M=' num2str(M) ', 50x25 + 100x50 + 200x100'] );
figure(2);
plot( x, vc/U, 'DisplayName', ...
    ['MRT, M=' num2str(M) ', 50x25 + 100x50 + 200x100'] );

% % MRT 30 / 60 / 120
% % Refined
% M = 0.2;
% [x, y, uc, vc] = get_centerlines( ...
%     '../results/centerlines/centerlines_G30_M200_T50_RM3_VCM0_VCF0_Re1000.tsv');
% U = M / sqrt(3);
% figure(1);
% plot( y, uc/U, 'DisplayName', ...
%     ['MRT, M=' num2str(M) ', 30x10 + 60x20 + 120x40'] );
% figure(2);
% plot( x, vc/U, 'DisplayName', ...
%     ['MRT, M=' num2str(M) ', 30x10 + 60x20 + 120x40'] );

% % MRT 45 / 90 / 180
% % Refined
% M = 0.2;
% [x, y, uc, vc] = get_centerlines( ...
%     '../results/centerlines/centerlines_G45_M200_T50_RM3_VCM0_VCF0_Re1000.tsv');
% U = M / sqrt(3);
% figure(1);
% plot( y, uc/U, 'DisplayName', ...
%     ['MRT, M=' num2str(M) ', 45x15 + 90x30 + 180x60'] );
% figure(2);
% plot( x, vc/U, 'DisplayName', ...
%     ['MRT, M=' num2str(M) ', 45x15 + 90x30 + 180x60'] );

% MRT 60 / 120 / 240
% Refined
M = 0.2;
[x, y, uc, vc] = get_centerlines( ...
    '../results/centerlines/centerlines_G60_M200_T50_RM3_VCM0_VCF0_Re1000.tsv');
U = M / sqrt(3);
figure(1);
plot( y, uc/U, 'DisplayName', ...
    ['MRT, M=' num2str(M) ', 60x20 + 120x40 + 240x80'] );
figure(2);
plot( x, vc/U, 'DisplayName', ...
    ['MRT, M=' num2str(M) ', 60x20 + 120x40 + 240x80'] );

% % MRT 36 / 72 / 144
% % Refined
% M = 0.2;
% [x, y, uc, vc] = get_centerlines( ...
%     '../results/centerlines/centerlines_G36_M200_T50_RM3_VCM0_VCF0_Re1000.tsv');
% U = M / sqrt(3);
% figure(1);
% plot( y, uc/U, 'DisplayName', ...
%     ['MRT, M=' num2str(M) ', 36x12 + 72x24 + 144x48'] );
% figure(2);
% plot( x, vc/U, 'DisplayName', ...
%     ['MRT, M=' num2str(M) ', 36x12 + 72x24 + 144x48'] );

% % MRT 39 / 78 / 156
% % Refined
% M = 0.2;
% [x, y, uc, vc] = get_centerlines( ...
%     '../results/centerlines/centerlines_G39_M200_T50_RM3_VCM0_VCF0_Re1000.tsv');
% U = M / sqrt(3);
% figure(1);
% plot( y, uc/U, 'DisplayName', ...
%     ['MRT, M=' num2str(M) ', 39x13 + 78x26 + 156x52'] );
% figure(2);
% plot( x, vc/U, 'DisplayName', ...
%     ['MRT, M=' num2str(M) ', 39x13 + 78x26 + 156x52'] );

% % MRT 40x40
% M = 0.2;
% u = dlmread('../results/fields/u_G40_M200_T50_RM3_VCM0_VCF0_Re1000.dat');
% v = dlmread('../results/fields/v_G40_M200_T50_RM3_VCM0_VCF0_Re1000.dat');
% U = M / sqrt(3);
% figure(1);
% [y, uc] = get_centerline_u(u,U);
% plot( y, uc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 40x40'] );
% figure(2);
% [x, vc] = get_centerline_v(v);
% plot( x, vc/U, 'DisplayName', ['MRT, M=' num2str(M) ', 40x40'] );

figure(1);
legend('show');
figure(2);
legend('show');
