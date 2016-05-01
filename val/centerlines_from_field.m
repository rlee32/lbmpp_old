clear; clc;
% close all;

addpath subroutines

% This reads centerlines from uniform fields.

% Change the dlmread entry according to the desired solution file.
% The solution file nomencalture descriptions can be found in 
%  the README.txt in the '../results/' folder. 

% % SRT 100 Re
% M = 0.1;
% Re = 100;
% u = dlmread('../results/u_G101_M100_T100_RM1_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/v_G101_M100_T100_RM1_VCM0_VCF0_Re100.dat');

% SRT 1000 Re
M = 0.1;
Re = 1000;
u = dlmread('../results/fields/u_G101_M100_T100_RM1_VCM0_VCF0_Re1000.dat');
v = dlmread('../results/fields/v_G101_M100_T100_RM1_VCM0_VCF0_Re1000.dat');

% SRT + SVC 1000 Re
M = 0.1;
Re = 1000;
u = dlmread('../results/fields/u_G101_M100_T100_RM1_VCM2_VCF2_Re1000.dat');
v = dlmread('../results/fields/v_G101_M100_T100_RM1_VCM2_VCF2_Re1000.dat');

% % SRT 2500 Re
% M = 0.1;
% Re = 2500;
% u = dlmread('../results/u_G151_M100_T125_RM1_VCM0_VCF0_Re2500.dat');
% v = dlmread('../results/v_G151_M100_T125_RM1_VCM0_VCF0_Re2500.dat');

% % MRT 100 Re
% M = 0.1;
% Re = 100;
% u = dlmread('../results/fields/u_G175_M100_T100_RM3_VCM0_VCF0_Re100.dat');
% v = dlmread('../results/fields/v_G175_M100_T100_RM3_VCM0_VCF0_Re100.dat');

% MRT 1000 Re
M = 0.1;
Re = 1000;
u = dlmread('../results/fields/u_G101_M100_T100_RM3_VCM0_VCF0_Re1000.dat');
v = dlmread('../results/fields/v_G101_M100_T100_RM3_VCM0_VCF0_Re1000.dat');

% MRT 1000 Re + VC
M = 0.1;
Re = 1000;
u = dlmread('../results/fields/u_G101_M100_T100_RM3_VCM2_VCF1_Re1000.dat');
v = dlmread('../results/fields/v_G101_M100_T100_RM3_VCM2_VCF1_Re1000.dat');

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
rows_even = mod(rows, 2) == 0;
cols_even = mod(cols, 2) == 0;

dx = 1 / cols;
dy = 1 / rows;
H = 1; % height of cavity.
y = [0, ( linspace(dy/2,H-dy/2,rows) ), H];

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

validation = H == 1 && (Re==100 || Re==10000 || Re==1000 || Re==21000);

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
[y, uc] = get_centerline_u(u,U);
plot( y, uc/U, 'DisplayName', 'Present LBM' );
% plot( y, [0; (mean(u(:,cols),2)/U); 1] );
title(['Centerline y vs. u-velocity at Re = ' num2str(Re) ...
    ', M = ' num2str(M)]);
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

[x, vc] = get_centerline_v(v);
plot( x, vc/U );
% plot( x, mean(v(rows,:),1)/U );
title(['Centerline x vs. v-velocity at Re = ' num2str(Re) ...
    ', M = ' num2str(M)]);
xlabel('x');
ylabel('v');
if validation
    legend('Ghia et al', 'Present LBM');
end
