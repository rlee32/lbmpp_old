function error = reference_error(x,v,Re)
% compares the input x,v to the reference x,v (Ghia et al)
% returns the l2 norm of the error.
% x,v are the x-location and v component of velocity along the centerline.
% interp1 is used to compare the two sets of data. 
% We assume the solution data is dense, while the reference is sparse.
% Therefore, we use the reference data to query interpolations in the
%   solution data.

[xref, vref] = validation_data_v_vs_x(Re);

% Centerline
[rows, ~, ~] = size(v);
if mod(rows, 2) == 0
    c = ( sum(v(round(rows/2):round(rows/2+1),:),1) ) / 2;
else
    c = v(ceil(rows/2),:);
end

vq = interp1(x,c,xref);
error = sqrt(sum((vq - vref).^2));



