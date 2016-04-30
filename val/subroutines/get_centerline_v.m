function [x, v] = get_centerline_v( vfield )

% From the uniform field vfield, extracts the centerline v.
% also returns corresponding x, assuming x in [0,1].

[ rows, cols ] = size(vfield);
rows_even = mod(rows, 2) == 0;
if rows_even
    rows = [rows/2, rows/2+1];
else
    rows = ceil(rows/2);
end

dx = 1 / cols;
x = [0, ( linspace(dx/2,1-dx/2,cols) ), 1];
v = [0, mean( vfield(rows,:), 1 ), 0];