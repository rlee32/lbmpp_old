clear; close all; clc;

% U = 0.006 / sqrt(3);
% H = 1.5;
% Re = 100;
% u = dlmread('../results/u_128x192_6_250_1_0_1.dat');
% v = dlmread('../results/v_128x192_6_250_1_0_1.dat');

% U = 0.006 / sqrt(3);
% H = 1;
% Re = 100;
% u = dlmread('../results/u_128_6_250_3_1_10.dat');
% % v = dlmread('../results/v_128_6_250_3_1_10.dat');

% U = 0.03 / sqrt(3);
% H = 1;
% Re = 10000;
% u = dlmread('../results/u_128_30_250_3_1_10_Re10000.dat');
% v = dlmread('../results/v_128_30_250_3_1_10_Re10000.dat');

% U = 0.015 / sqrt(3);
% H = 1;
% Re = 500;
% u = dlmread('../results/u_128_15_250_1_0_10_Re500.dat');
% v = dlmread('../results/v_128_15_250_1_0_10_Re500.dat');

U = 0.015 / sqrt(3);
H = 1.5;
Re = 500;
u = dlmread('../results/u_128x192_15_250_1_0_10_Re500.dat');
v = dlmread('../results/v_128x192_15_250_1_0_10_Re500.dat');

[rows,cols] = size(u);
rows_even = mod(rows, 2) == 0;
cols_even = mod(cols, 2) == 0;

dx = 1 / cols;
dy = 1 / rows;
% x = ( linspace(0,1-dx,cols) ) / ( 1 - dx );
% y = ( linspace(0,1-dy,rows) ) / ( 1 - dy );
x = ( linspace(0,1,cols) );
y = ( linspace(0,H,rows) );

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
plot( y, flipud(mean(u(:,cols),2)/U) );
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
streamslice(X,Y, flipud(u),flipud(v));
% streamline( stream2(X,Y,u,v,x,y) );
% streamline( X,Y,u,v,x,y );
title(['Steady-State Streamlines at Re = ' num2str(Re)]);
xlabel('x');
ylabel('y');
