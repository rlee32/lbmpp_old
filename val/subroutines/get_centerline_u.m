function [y, u] = get_centerline_u( ufield, U )

% From the uniform field ufield, extracts the centerline u.
% also returns corresponding y, assuming y in [0,1].

[ rows, cols ] = size(ufield);
cols_even = mod(cols, 2) == 0;
if cols_even
    cols = [cols/2, cols/2+1];
else
    cols = ceil(cols/2);
end

dy = 1 / rows;
y = [0, ( linspace(dy/2,1-dy/2,rows) ), 1];
u = [0; mean( ufield(:,cols), 2 ); U];