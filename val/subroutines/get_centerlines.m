function [x, y, u, v] = get_centerlines( file_name )

% Simply reads in the centerlines as output by the simulator.

A = dlmread(file_name);
x = A(1,:);
y = A(3,:);
u = A(4,:);
v = A(2,:);


