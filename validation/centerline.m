clear; close all; clc;

U = 0.003 * 4;

u = dlmread('../u_velocity_field.dat');
v = dlmread('../v_velocity_field.dat');

[rows,cols] = size(u);
rows_even = mod(rows, 2) == 0;
cols_even = mod(cols, 2) == 0;

dx = 1 / cols;
dy = 1 / rows;
x = ( linspace(0,1-dx,cols) ) / ( 1 - dx );
y = ( linspace(0,1-dy,rows) ) / ( 1 - dy );

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


[x_ref, v_ref] = validation_data_v_vs_x();
[y_ref, u_ref] = validation_data_u_vs_y();

figure;
plot(y_ref,u_ref, 'x');
hold on;
plot( y, flipud(mean(u(:,cols),2)/U) );
title('Centerline y vs. u-velocity at Re = 100');
xlabel('y');
ylabel('u');
legend('Ghia et al', 'Present LBM');

figure;
plot(x_ref,v_ref, 'x');
hold on;
plot( x, mean(v(rows,:),1)/U );
title('Centerline x vs. v-velocity at Re = 100');
xlabel('x');
ylabel('v');
legend('Ghia et al', 'Present LBM');

figure;
[X,Y] = meshgrid(x,y);
axis equal tight;
streamslice(X,Y, flipud(u),flipud(v));
% streamline( stream2(X,Y,u,v,x,y) );
% streamline( X,Y,u,v,x,y );
title('Steady-State Streamlines at Re = 100');
xlabel('x');
ylabel('y');
